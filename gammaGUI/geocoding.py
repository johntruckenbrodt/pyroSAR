##############################################################
# Geocoding, topographic normalization and kml generation of
# intensity, coherence and polarimetric decomposition images
# module of software gammaGUI
# John Truckenbrodt 2015
##############################################################

"""
The following tasks are performed by executing this script:
-group all images belonging to one sar scene (i.e. all images in the scene folder plus interferometric products for which the individual scene was used as master)
 into an easily accessible object structure
-select a master mli from a priority list
-create lookup tables for matching the sar geometry to that of a reference sar image simulated from an earlier created dem file
-topographically normalize all mli images (including pauli decomp), fd3c and cloude decompositions
-perform backward geocoding for all mlis, coherence images and polarimetric decompositions
-transform all geocoded images to EQA projection and create bitmaps and kmls from them

please refer to the automatically created logfiles (subfolder LOG) for details on the process performance; an error is only thrown if the snr between master mli and the simulated
sar image is below the defined threshold (command init_offsetm)
# in the current implementation the dem is not being oversampled and, hence, needs to be in the required resolution before starting geocoding.
# the provided oversampling factor is for cross-correlation only
"""

import os
import re
import sys
from ancillary import dissolve, ReadPar
from envi import hdr
from auxiliary import Tuple, grouping, Environment
from gamma.util import ISPPar, correlate, init_offset, gamma

# create list of scene tuple objects
tuples = grouping()

# name of the dem to be used for geocoding
dem = sys.argv[1]

# define (and create) directory for logfile
path_log = os.path.join(os.getcwd(), "LOG/GEO/")
if not os.path.exists(path_log):
    os.makedirs(path_log)

# read processing parameter textfile
par = ReadPar(os.path.join(os.getcwd(), "PAR/geocoding.par"), type="exe")
par.interp_mode = Environment.dropoptions["interp_mode"].index(par.interp_mode)+1

# perform topographic normalization?
tnorm = True if par.topographic_normalization == "True" else False

for scene in tuples:

    # image ids in descending order of priority for master selection
    prioritytags = ["HV[_a-z]*_mli$", "VH[_a-z]*_mli$", "HH[_a-z]*_mli$", "VV[_a-z]*_mli$"]

    # select master
    master = ""
    for tag in prioritytags:
        for key in scene.__dict__.keys():
            if re.match(tag, key):
                master = getattr(scene, key)
                break

    if len(master) == 0:
        raise IOError("no appropriate master image found; consider multilooking")

    path_out = os.path.join(os.path.dirname(master), "GEO/")
    if not os.path.exists(path_out):
        os.makedirs(path_out)

    print "#############################################"
    print "processing of scene " + os.path.dirname(master).split("/")[-1] + " started"
    print "----------\nforward geocoding\n----------"
    print "selected master: " + os.path.basename(master)
    # concatenate names for files to be created
    basename = os.path.join(path_out+os.path.basename(master))
    coffs = basename+"_coffs"
    coffsets = basename+"_coffsets"
    dem_sub = basename+"_dem"
    dempar_eqa = basename+"_dem_eqa.par"
    diffpar = basename+"_diff.par"
    inc = basename+"_inc"
    lut_rough = basename+"_rough_map_to_rdc"
    lut_fine = basename+"_map_to_rdc"
    ls_map = basename+"_ls_map"
    offs = basename+"_offs"
    offsets = basename+"_offsets"
    # pix = basename+"_pix"
    # psi = basename+"_psi"
    sim_map = basename+"_sim_map"
    sim_rdc = basename+"_sim_rdc"
    snr = basename+"_snr"
    slope = basename+"_u"
    aspect = basename+"_v"

    # define parameter file of the master
    masterpar_name = master + ".par"
    masterpar = ISPPar(masterpar_name)

    # perform forward geocoding if the corresponding refined lookup table does not yet exist
    if not os.path.isfile(lut_fine):

        print"creating initial lookup table"
        gamma(["gc_map", masterpar_name, "-", dem + ".par", dem, dem_sub + ".par", dem_sub, lut_rough, "-", "-", sim_map, slope, aspect, inc, "-", "-", ls_map, par.frame, 2], path_out, path_log)

        for header in [inc+".hdr", ls_map+".hdr", slope+".hdr", aspect+".hdr"]:
            hdr(dem_sub+".par", header)

        if tnorm:
            print "initial pixel area estimation"
            topo_base = os.path.join(path_out, os.path.basename(master))
            pix_gamma = topo_base+"_pix_gamma"
            parfile = master+".par"
            gamma(["pixel_area", parfile, dem_sub + ".par", dem_sub, lut_rough, ls_map, inc, "-", pix_gamma], path_out, path_log)

        # read additional parameters
        dempar = ISPPar(dem_sub + ".par")
        samples_dem = dempar.width
        samples_mli = masterpar.range_samples
        lines_mli = masterpar.azimuth_lines

        print "refinement of lookup table"
        gamma(["geocode", lut_rough, sim_map, samples_dem, sim_rdc, samples_mli, lines_mli, par.interp_mode], path_out, path_log)

        print "initial computation of offsets"
        init_offset(master, sim_rdc, diffpar, path_log)

        print "cross-correlation offset and polynomial estimation"
        correlate(master, sim_rdc, diffpar, offs, offsets, snr, coffs, coffsets, path_log, maxwin=2048, overlap=.3, poly=par.polynomial, ovs=par.oversampling)

        print "refinement of lookup table with offset polynomials"
        gamma(["gc_map_fine", lut_rough, samples_dem, diffpar, lut_fine, 1], path_out, path_log)

        for item in [lut_rough, sim_map, sim_rdc]:
            os.remove(item)
    else:
        # read additional parameters
        samples_dem = ISPPar(dem_sub + ".par").width
        samples_mli = masterpar.range_samples
    # perform topographic normalization
    if tnorm:
        print "----------"
        # define regular expressions for images to be normalized (mli, fd3c decomp, cloude decomp)
        normtags = ["_mli$", "fdd_p[dsv]$", "ctd_[123]_mag$"]
        # create a list containing the actual image names including their absolute directory paths (only if the corresponding normalized processing result does not yet exist)
        normlist = dissolve([getattr(scene, x) for x in scene.index for tag in normtags if re.search(tag, x)])
        normlist = [x for x in normlist if not os.path.isfile(x+"_norm")]
        if len(normlist) > 0:
            print "topographic normalization\n----------"
            for image in normlist:
                print os.path.basename(image)
                topo_base = os.path.join(path_out, os.path.basename(image))
                pix_gamma = topo_base+"_pix_gamma"
                pix_gamma0ratio = topo_base+"_pix_gamma0ratio"
                ellip_pix_gamma0 = topo_base+"_ellip_pix_gamma0"
                parfile = image+".par" if re.search("_mli$", image) else scene.getTop("HH_(?:slc_|)(?:cal_|)mli$")+".par"
                # refined pixel area estimation
                gamma(["pixel_area", parfile, dem_sub + ".par", dem_sub, lut_fine, ls_map, inc, "-", pix_gamma], path_out, path_log)
                gamma(["radcal_MLI", image, parfile, "-", "temp", "-", "-", "-", 2, "-", "-", ellip_pix_gamma0], path_out, path_log)
                width_mli = ISPPar(parfile).range_samples
                gamma(["ratio", ellip_pix_gamma0, pix_gamma, pix_gamma0ratio, width_mli, 1, 1, 0], path_out, path_log)
                gamma(["product", image, pix_gamma0ratio, image + "_norm", width_mli, 1, 1, 0], path_out, path_log)
                # remove temporary files
                for item in [os.path.join(path_out, "temp"), pix_gamma, ellip_pix_gamma0, pix_gamma0ratio]:
                    os.remove(item)
                # update scene object
                scene = Tuple(scene.main)
        else:
            print "topographic normalization already performed"

    # define regular expressions for all images to be geocoded
    # [mli images (normalized mli if tnorm is True), coherence images, H-A-alpha decomp, cloude decomp, krogager decomp, fd3c decomp, kennaugh decomp]
    geotags = ["_norm$" if tnorm else "_mli$", "cc_(?:ad|wave)$", "cpd_(?:A|alpha|H)$", "ctd_[123]_mag$", "(?:diplane|helix|sphere)$", "fdd_p[dsv]$", "k[1234]{2}$"]

    # create list of processing candidates
    geolist = [x for x in dissolve([getattr(scene, x) for x in scene.index for tag in geotags if re.search(tag, x) and not re.search("reg", x)]) if not os.path.isfile(x+"_geo")]

    # continue with next scene if all candidates have already been processed
    if len(geolist) == 0:
        print "----------\ngeocoding already performed"
        continue
    else:
        print "----------\nbackward geocoding\n----------"

        # copy parameter file of dem in eqa projection
        # if not os.path.isfile(dempar_eqa):
        #     sp.check_call(["cp", os.path.join(os.getcwd(), "DEM/dem_final.par"), dempar_eqa])

        for image in geolist:
            print os.path.basename(image)
            interp_mode = par.interp_mode_coherence if "_cc" in image else par.interp_mode_intensity

            gamma(["geocode_back", image, samples_mli, lut_fine, image + "_geo", samples_dem, "-", interp_mode], path_out, path_log)

            # create envi header file
            hdr(dem_sub+".par", image+"_geo.hdr")

            # geoimage = image+"_geo"
            # eqaimage = geoimage+"_eqa"
            #
            # # transform image to eqa projection
            # gamma(["map_trans", dem_sub+".par", geoimage, dempar_eqa, eqaimage, "-", "-", "-", "-", 1], os.path.dirname(geoimage), path_log)
            #
            # width_eqa = ISPPar(dempar_eqa).width
            #
            # # create bitmap
            # if "cc" in geoimage:
            #     gamma(["rascc", eqaimage, "-", width_eqa, 1, 1, 0, 1, 1, .1, .9, 1.0, .35, 1, eqaimage + ".bmp"], os.path.dirname(geoimage), path_log)
            # else:
            #     gamma(["raspwr", eqaimage, width_eqa, 1, 0, 1, 1, 1, .35, 1.0, eqaimage + ".bmp", "0"], os.path.dirname(geoimage), path_log)
            #
            # # create kml
            # gamma(["mk_kml", dempar_eqa, eqaimage + ".bmp", eqaimage + ".kml"], os.path.dirname(geoimage), path_log)

print "#############################################\ngeocoding finished"
