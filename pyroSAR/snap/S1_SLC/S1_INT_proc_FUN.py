##necessary libraries
import pyroSAR
from pyroSAR.snap.auxil import parse_recipe, parse_node, gpt, execute
import os
import shutil
import glob
import datetime


def S1_HA_proc(infiles, out_dir=None, tmpdir=None, t_res=20, t_crs=4326, out_format="GeoTIFF", gpt_paras=None,
               IWs=["IW1", "IW2", "IW3"], decompFeats=["Alpha", "Entropy", "Anisotropy"], ext_DEM=False,
               ext_DEM_noDatVal=-9999, ext_Dem_file=None, msk_noDatVal=False, ext_DEM_EGM=True,
               imgResamp="BICUBIC_INTERPOLATION", demResamp="BILINEAR_INTERPOLATION", decomp_win_size=5,
               speckFilter="Box Car Filter", ml_RgLook=4, ml_AzLook=1, osvPath=None, firstBurstIndex=None,
               lastBurstIndex=None, clean_tmpdir=True):
    """

    function for processing backscatter intensities VV and VH from S-1 SLC files in SNAP

    Parameters
    ----------
        infiles: list or str
            filepaths of SLC zip files
        out_dir: str or None
            output folder if None a default folder structure is provided: "INT/decompFeat/"
        tmpdir: str
            temporary dir for intermediate processing steps, its automatically created at cwd if none is provided
        t_res: int, float
            resolution in meters of final product, default is 20
        t_crs: int
            EPSG code of target coordinate system, default is 4326
        out_format: str
            format of final output, formats supported by SNAP, default is GeoTiff
        gpt_paras: none or list
            a list of additional arguments to be passed to the gpt call
        pol: str or list or "full"
            polaristations to process, "full" processes all available polarizations, default is "full"
        IWs: str or list
            selected subswath for processing, default is all 3
        extDEM: bool
            set to true if external DEM should be used in processing
        ext_DEM_noDatVal: int or float
            dependent on external DEM, default False
        ext_DEM_file: str
            path to file of external DEM, must be a format that SNAP can handle
        msk_NoDatVal: bool
            if true No data values of DEM, especially at sea, are masked out
        ext_DEM_EGM: bool
            apply earth gravitational model to external DEM, default true
        imgResamp: str
            image resampling method, must be supported by SNAP
        demResamp: str
            DEM resampling method, must be supported by SNAP
        speckFilter: str
            type of speckle filtering approach, default is Boxcar
        filterSizeX: int
            window size of speckle filter in x, default is 5
        filterSizeY: int
            window size of speckle filter in y, default is 5
        ml_RgLook: int
            number of looks in range, default is 4
        ml_AzLook: int
            number of looks in azimuth, default is 1
        firstBurstIndex: int or None
            index of first burst for TOPSAR-split
        lastBurstIndex: int or None
            index of last burst for TOPSAR-split
        l2dB: bool
            option for conversion from linear to dB scaling of output, default true
        clean_tmpdir, bool
            delete tmpdir, default true
        osvPath: None
            specify path to locally stored OSVs, if none default OSV path of SNAP is set

        Returns
        -------
        Raster files of selected output format for selected H-alpha features

        Note
        ----
        Only set first and last burstindex if all files you are processing have the same number of bursts

        Examples
        --------
        process backscatter intensities VV and VH for given SLC file

        >>> from pyroSAR.snap import S1_INT_proc
        >>> filename= 'S1A_IW_GRDH_1SDV_20180829T170656_20180829T170721_023464_028DE0_F7BD.zip'
        >>> gpt_paras = ["-e", "-x", "-c","35G", "-q", "16", "-J-Xms25G", "-J-Xmx75G"]
        >>> pol= "full"
        >>> S1_INT_proc(infiles= filename, gtp_paras= gpt_paras, pol= "full")

    """
    ##define formatName for reading zip-files
    formatName = "SENTINEL-1"
    ##specify tmp output format
    tpm_format = "BEAM-DIMAP"
    ## create temp dir for intermediate .dim files
    if tmpdir is None:
        tmpdir = os.getcwd() + "/tmp_dim"
        if os.path.isdir(tmpdir) == False:
            os.mkdir(tmpdir)
    ##check if a single IW or consecutive IWs are selected
    if isinstance(IWs, str):
        IWs = [IWs]
    if sorted(IWs) == ["IW1", "IW3"]:
        raise RuntimeError("Please select single or consecutive IW")
    
    ##extract info about files and order them by date
    ##handle length and type of infiles: str or list
    if isinstance(infiles, str):
        info = pyroSAR.identify(infiles)
        fps_lst = [info.scene]
        info = [info]
    elif isinstance(infiles, list):
        info = pyroSAR.identify_many(infiles, verbose=False, sortkey='start')
        ##collect filepaths sorted by date
        fps_lst = []
        for fp in info:
            fp_str = fp.scene
            fps_lst.append(fp_str)
    else:
        raise RuntimeError('Please provide str or list of filepaths')
        ##query and handle polarisations, raise error if selected polarisations don't match (see Truckenbrodt et al.: pyroSAR: geocode)
    ##specify auto download DEM and handle external DEM file
    if ext_DEM == False:
        demName = 'SRTM 1Sec HGT'
        ext_DEM_file = None
    else:
        demName = "External DEM"
    ##raise error if no path to external file is provided
    if ext_DEM == True and ext_DEM_file == None:
        raise RuntimeError('No DEM file provided. Specify path to DEM-file')
    ##raise error ifwrong decomp feature
    
    ##handle SNAP problem with WGS84 (EPSG: 4326) by manually constructing crs string (see Truckenbrodt et al.: pyroSAR: geocode)
    if t_crs == 4326:
        epsg = 'GEOGCS["WGS84(DD)",''DATUM["WGS84",''SPHEROID["WGS84", 6378137.0, 298.257223563]],''PRIMEM["Greenwich", 0.0],''UNIT["degree", 0.017453292519943295],''AXIS["Geodetic longitude", EAST],'           'AXIS["Geodetic latitude", NORTH]]'
    else:
        epsg = "EPSG:{}".format(t_crs)
    ##check if correct DEM resampling methods are supplied
    reSamp_LookUp = ['NEAREST_NEIGHBOUR',
                     'BILINEAR_INTERPOLATION',
                     'CUBIC_CONVOLUTION',
                     'BISINC_5_POINT_INTERPOLATION',
                     'BISINC_11_POINT_INTERPOLATION',
                     'BISINC_21_POINT_INTERPOLATION',
                     'BICUBIC_INTERPOLATION']
    
    message = '{0} must be one of the following:\n- {1}'
    if demResamp not in reSamp_LookUp:
        raise ValueError(message.format('demResamplingMethod', '\n- '.join(reSamp_LookUp)))
    if imgResamp not in reSamp_LookUp:
        raise ValueError(message.format('imgResamplingMethod', '\n- '.join(reSamp_LookUp)))
    ##check if correct speckle filter option is supplied
    speckleFilter_options = ['Box Car Filter', 'IDAN Filter', 'Refined Lee Filter', 'Improved Lee Sigma Filter']
    if speckFilter not in speckleFilter_options:
        raise ValueError(message.format('speckleFilter', '\n- '.join(speckleFilter_options)))
    ##query unique dates of files: determine if sliceAssembly is required
    dates_info = []
    for d in info:
        di = d.start.split("T")[0]
        dates_info.append(di)
    
    unique_dates_info = list(set(dates_info))
    unique_dates_info = sorted(unique_dates_info, key=lambda x: datetime.datetime.strptime(x, '%Y%m%d'))
    
    ##check for files of the same date and put them in sublists
    pair_dates_idx = []
    
    for a in unique_dates_info:
        tmp_dates = []
        for idx, elem in enumerate(dates_info):
            if (a == elem):
                tmp_dates.append(idx)
        
        pair_dates_idx.append(tmp_dates)
    
    ##selection of paired files for sliceAssembly
    for i in range(0, len(pair_dates_idx)):
        fps_grp = list(map(fps_lst.__getitem__, pair_dates_idx[i]))
        # get relative orbit number of grouped files
        info_tmp = pyroSAR.identify(fps_grp[0])
        relOrb = info_tmp.orbitNumber_rel
        sensor = info_tmp.sensor
        orbit = info_tmp.orbit
        pol = info_tmp.polarizations
        date_str = info_tmp.start
        
        ##check availability of orbit state vector
        orbitType = "Sentinel Precise (Auto Download)"
        
        match = info_tmp.getOSV(osvType='POE', returnMatch=True, osvdir=osvPath)
        if match is None:
            info_tmp.getOSV(osvType='RES', osvdir=osvPath)
            orbitType = 'Sentinel Restituted (Auto Download)'
        ##exception handling of SNAP errors    
        try:
            
            slcAs_name = sensor + "_relOrb_" + str(relOrb) + "_HA_" + unique_dates_info[i] + "_slcAs"
            slcAs_out = os.path.join(tmpdir, slcAs_name)
            ## create workflow for sliceAssembly if more than 1 file is available per date
            if len(fps_grp) > 1:
                
                workflow_slcAs = parse_recipe("blank")
                
                read1 = parse_node('Read')
                read1.parameters['file'] = fps_grp[0]
                read1.parameters['formatName'] = formatName
                readers = [read1.id]
                
                workflow_slcAs.insert_node(read1)
                
                for r in range(1, len(fps_grp)):
                    readn = parse_node('Read')
                    readn.parameters['file'] = fps_grp[r]
                    readn.parameters['formatName'] = formatName
                    workflow_slcAs.insert_node(readn, before=read1.id, resetSuccessorSource=False)
                    readers.append(readn.id)
                
                slcAs = parse_node("SliceAssembly")
                slcAs.parameters["selectedPolarisations"] = pol
                
                workflow_slcAs.insert_node(slcAs, before=readers)
                read1 = slcAs
                
                write_slcAs = parse_node("Write")
                write_slcAs.parameters["file"] = slcAs_out
                write_slcAs.parameters["formatName"] = tpm_format
                
                workflow_slcAs.insert_node(write_slcAs, before=slcAs.id)
                
                workflow_slcAs.write("HA_slc_prep_graph")
                
                gpt('HA_slc_prep_graph.xml', gpt_args=gpt_paras, outdir=tmpdir)
                
                HA_proc_in = slcAs_out + ".dim"
            ##pass file path if no sliceAssembly required
            else:
                HA_proc_in = fps_grp[0]
            
            for iw in IWs:
                tpm_name = sensor + "_HA_relOrb_" + str(relOrb) + "_" + unique_dates_info[i] + "_" + iw + "_2TPM"
                tpm_out = os.path.join(tmpdir, tpm_name)
                ##generate workflow for IW splits 
                workflow = parse_recipe("blank")
                
                read = parse_node("Read")
                read.parameters["file"] = HA_proc_in
                if len(fps_grp) == 1:
                    read.parameters["formatName"] = formatName
                workflow.insert_node(read)
                
                aof = parse_node("Apply-Orbit-File")
                aof.parameters["orbitType"] = orbitType
                aof.parameters["polyDegree"] = 3
                aof.parameters["continueOnFail"] = False
                workflow.insert_node(aof, before=read.id)
                ##TOPSAR split node
                ts = parse_node("TOPSAR-Split")
                ts.parameters["subswath"] = iw
                if firstBurstIndex is not None and lastBurstIndex is not None:
                    ts.parameters["firstBurstIndex"] = firstBurstIndex
                    ts.parameters["lastBurstIndex"] = lastBurstIndex
                workflow.insert_node(ts, before=aof.id)
                
                cal = parse_node("Calibration")
                cal.parameters["selectedPolarisations"] = pol
                cal.parameters["createBetaBand"] = False
                cal.parameters["outputBetaBand"] = False
                cal.parameters["outputSigmaBand"] = True
                cal.parameters["outputImageInComplex"] = True
                workflow.insert_node(cal, before=ts.id)
                
                tpd = parse_node("TOPSAR-Deburst")
                tpd.parameters["selectedPolarisations"] = pol
                workflow.insert_node(tpd, before=cal.id)
                
                write_tmp = parse_node("Write")
                write_tmp.parameters["file"] = tpm_out
                write_tmp.parameters["formatName"] = tpm_format
                workflow.insert_node(write_tmp, before=tpd.id)
                
                workflow.write("HA_proc_IW_graph")
                
                execute('HA_proc_IW_graph.xml', gpt_args=gpt_paras)
            
            for dc in decompFeats:
                dc_label = dc.upper()[0:3]
                ##load temporary files
                tpm_in = glob.glob(
                    tmpdir + "/" + sensor + "_HA_relOrb_" + str(relOrb) + "_" + unique_dates_info[i] + "*_2TPM.dim")
                ## parse_workflow of INT processing
                workflow_tpm = parse_recipe("blank")
                
                read1 = parse_node('Read')
                read1.parameters['file'] = tpm_in[0]
                workflow_tpm.insert_node(read1)
                last_node = read1.id
                ##merge IWs if multiple IWs were selected
                if len(tpm_in) > 1:
                    readers = [read1.id]
                    
                    for t in range(1, len(tpm_in)):
                        readn = parse_node('Read')
                        readn.parameters['file'] = tpm_in[t]
                        workflow_tpm.insert_node(readn, before=last_node, resetSuccessorSource=False)
                        readers.append(readn.id)
                    ##TOPSAR merge     
                    tpm = parse_node("TOPSAR-Merge")
                    tpm.parameters["selectedPolarisations"] = pol
                    workflow_tpm.insert_node(tpm, before=readers)
                    last_node = tpm.id
                ##create C2 covariance matrix
                polMat = parse_node("Polarimetric-Matrices")
                polMat.parameters["matrix"] = "C2"
                workflow_tpm.insert_node(polMat, before=last_node)
                last_node = polMat.id
                ##multi looking
                ml = parse_node("Multilook")
                ml.parameters["sourceBands"] = ["C11", "C12_real", "C12_imag", "C22"]
                ml.parameters["nRgLooks"] = ml_RgLook
                ml.parameters["nAzLooks"] = ml_AzLook
                ml.parameters["grSquarePixel"] = True
                ml.parameters["outputIntensity"] = False
                workflow_tpm.insert_node(ml, before=last_node)
                last_node = ml.id
                
                ##polaricmetric speckle filtering
                polSpec = parse_node("Polarimetric-Speckle-Filter")
                polSpec.parameters["filter"] = speckFilter
                workflow_tpm.insert_node(polSpec, before=last_node)
                last_node = polSpec.id
                
                ##dual-pol H/a decomposition
                polDecp = parse_node("Polarimetric-Decomposition")
                polDecp.parameters["decomposition"] = "H-Alpha Dual Pol Decomposition"
                polDecp.parameters["windowSize"] = decomp_win_size
                polDecp.parameters["outputHAAlpha"] = True
                
                workflow_tpm.insert_node(polDecp, before=last_node)
                last_node = polDecp.id
                
                # terrain correction
                tc = parse_node("Terrain-Correction")
                tc.parameters["sourceBands"] = [dc]
                tc.parameters["demName"] = demName
                tc.parameters["externalDEMFile"] = ext_Dem_file
                tc.parameters["externalDEMNoDataValue"] = ext_DEM_noDatVal
                tc.parameters["externalDEMApplyEGM"] = ext_DEM_EGM
                tc.parameters["imgResamplingMethod"] = imgResamp
                tc.parameters["demResamplingMethod"] = demResamp
                tc.parameters["pixelSpacingInMeter"] = t_res
                tc.parameters["mapProjection"] = epsg
                tc.parameters["saveSelectedSourceBand"] = True
                tc.parameters["outputComplex"] = False
                tc.parameters["nodataValueAtSea"] = msk_noDatVal
                
                workflow_tpm.insert_node(tc, before=last_node)
                last_node = tc.id
                
                ## generate str for final output file based on selected IWs
                if len(IWs) == 1:
                    out_name = sensor + "_" + orbit + "_relOrb_" + str(relOrb) + "_HA_" + dc_label + "_" + IWs[
                        0] + "_" + date_str + "_Orb_Cal_Deb_ML_TF_Spk_TC"
                elif len(IWs) == 2:
                    separator = "_"
                    iw_str = separator.join(IWs)
                    out_name = sensor + "_" + orbit + "_relOrb_" + str(
                        relOrb) + "_HA_" + dc_label + "_" + iw_str + "_" + date_str + "_Orb_Cal_Deb_ML_Spk_TC"
                else:
                    out_name = sensor + "_" + orbit + "_relOrb_" + str(
                        relOrb) + "_HA_" + dc_label + "_" + date_str + "_Orb_Cal_Deb_ML_Spk_TC"
                ##create default output folder for each selected polarization
                if out_dir is None:
                    out_dir_fp = "INT/" + dc_label
                    if os.path.isdir(out_dir_fp) == False:
                        os.makedirs(os.path.join(os.getcwd(), out_dir_fp))
                elif os.path.isdir(out_dir):
                    out_dir_fp = out_dir
                else:
                    raise RuntimeError("Please provide a valid filepath")
                
                out_path = os.path.join(out_dir_fp, out_name)
                
                write_tpm = parse_node("Write")
                write_tpm.parameters["file"] = out_path
                write_tpm.parameters["formatName"] = out_format
                workflow_tpm.insert_node(write_tpm, before=last_node)
                
                ##write graph and execute it
                workflow_tpm.write("HA_TPM_continued_proc_graph")
                
                execute('HA_TPM_continued_proc_graph.xml', gpt_args=gpt_paras)
        # exception for SNAP errors & creating error log
        except RuntimeError as e:
            print(str(e))
            with open("S1_INT_proc_ERROR_" + date_str + ".log", "w") as logf:
                logf.write(str(e))
            
            ##clean tmp folder to avoid overwriting errors even if exception is valid
            files = glob.glob(tmpdir + '/*')
            for f in files:
                if os.path.isfile(f) or os.path.islink(f):
                    os.unlink(f)
                elif os.path.isdir(f):
                    shutil.rmtree(f)
            
            continue
            
            ##clean tmp folder to avoid overwriting errors
        files = glob.glob(tmpdir + '/*')
        for f in files:
            if os.path.isfile(f) or os.path.islink(f):
                os.unlink(f)
            elif os.path.isdir(f):
                shutil.rmtree(f)
    ##delete tmp folder after processing
    if clean_tmpdir == True:
        shutil.rmtree(tmpdir)
