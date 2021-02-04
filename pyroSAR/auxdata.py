###############################################################################
# tools for handling auxiliary data in software pyroSAR

# Copyright (c) 2019-2021, the pyroSAR Developers.

# This file is part of the pyroSAR Project. It is subject to the
# license terms in the LICENSE.txt file found in the top-level
# directory of this distribution and at
# https://github.com/johntruckenbrodt/pyroSAR/blob/master/LICENSE.txt.
# No part of the pyroSAR project, including this file, may be
# copied, modified, propagated, or distributed except according
# to the terms contained in the LICENSE.txt file.
###############################################################################
import os
import re
import sys
import ftplib
import zipfile as zf
from ftplib import FTP
from time import strftime, gmtime
from os.path import expanduser

if sys.version_info >= (3, 0):
    from io import StringIO
    from urllib.request import urlopen
    from urllib.error import HTTPError
else:
    from cStringIO import StringIO
    from urllib2 import urlopen, HTTPError

from osgeo import gdal

from pyroSAR import identify
from pyroSAR.examine import ExamineSnap
from spatialist import Raster
from spatialist.ancillary import dissolve, finder
from spatialist.auxil import gdalbuildvrt, crsConvert, gdalwarp


def dem_autoload(geometries, demType, vrt=None, buffer=None, username=None, password=None, product='dem'):
    """
    obtain all relevant DEM tiles for selected geometries

    Parameters
    ----------
    geometries: list
        a list of :class:`spatialist.vector.Vector` geometries to obtain DEM data for;
        CRS must be WGS84 LatLon (EPSG 4326)
    demType: str
        the type of DEM to be used; current options:

        - 'AW3D30' (ALOS Global Digital Surface Model "ALOS World 3D - 30m")

          * url: ftp://ftp.eorc.jaxa.jp/pub/ALOS/ext1/AW3D30/release_v1804

        - 'SRTM 1Sec HGT'

          * url: https://step.esa.int/auxdata/dem/SRTMGL1

        - 'SRTM 3Sec'

          * url: https://srtm.csi.cgiar.org/wp-content/uploads/files/srtm_5x5/TIFF
        
        - 'TDX90m'
        
          * registration:  https://geoservice.dlr.de/web/dataguide/tdm90
          * url: ftpes://tandemx-90m.dlr.de

    vrt: str or None
        an optional GDAL VRT file created from the obtained DEM tiles
    buffer: int or float
        a buffer in degrees to add around the individual geometries
    username: str or None
        (optional) the user name for services requiring registration
    password: str or None
        (optional) the password for the registration account
    product: str
        the sub-product to extract from the DEM product.
        The following options are available for the respective DEM types:
        
        - 'AW3D30'
         
          * 'dem': the actual Digital Elevation Model
          * 'msk': mask information for each pixel (Cloud/Snow Mask, Land water and
            low correlation mask, Sea mask, Information of elevation dataset used
            for the void-filling processing)
          * 'stk': number of DSM-scene files which were used to produce the 5m resolution DSM
          
        - 'SRTM 1Sec HGT'
         
          * 'dem': the actual Digital Elevation Model
          
        - 'SRTM 3Sec'
         
          * 'dem': the actual Digital Elevation Model
          
        - 'TDX90m'
         
          * 'dem': the actual Digital Elevation Model
          * 'am2': Amplitude Mosaic representing the minimum value
          * 'amp': Amplitude Mosaic representing the mean value
          * 'com': Consistency Mask
          * 'cov': Coverage Map
          * 'hem': Height Error Map
          * 'lsm': Layover and Shadow Mask, based on SRTM C-band and Globe DEM data
          * 'wam': Water Indication Mask
    
    Returns
    -------
    list or None
        the names of the obtained files or None if a VRT file was defined
    
    Examples
    --------
    download all SRTM 1 arcsec DEMs overlapping with a Sentinel-1 scene and mosaic them to a single GeoTiff file
    
    .. code-block:: python
        
        from pyroSAR import identify
        from pyroSAR.auxdata import dem_autoload
        from spatialist import gdalwarp
        
        # identify the SAR scene
        filename = 'S1A_IW_SLC__1SDV_20150330T170734_20150330T170801_005264_006A6C_DA69.zip'
        scene = identify(filename)
        
        # extract the bounding box as spatialist.Vector object
        bbox = scene.bbox()
        
        # download the tiles and virtually combine them in an in-memory
        # VRT file subsetted to the extent of the SAR scene plus a buffer of 0.01 degrees
        vrt = '/vsimem/srtm1.vrt'
        dem_autoload(geometries=[bbox], demType='SRTM 1Sec HGT',
                     vrt=vrt, buffer=0.01)
        
        # write the final GeoTiff file
        outname = scene.outname_base() + 'srtm1.tif'
        gdalwarp(src=vrt, dst=outname, options={'format': 'GTiff'})
        
        # alternatively use function dem_create and warp the DEM to UTM
        # including conversion from geoid to ellipsoid heights
        from pyroSAR.auxdata import dem_create
        outname = scene.outname_base() + 'srtm1_ellp.tif'
        dem_create(src=vrt, dst=outname, t_srs=32632, tr=(30, 30),
                   geoid_convert=True, geoid='EGM96')
    """
    with DEMHandler(geometries) as handler:
        return handler.load(demType=demType,
                            username=username,
                            password=password,
                            vrt=vrt,
                            buffer=buffer,
                            product=product)


def dem_create(src, dst, t_srs=None, tr=None, resampling_method='bilinear', geoid_convert=False, geoid='EGM96'):
    """
    create a new DEM GeoTiff file and optionally convert heights from geoid to ellipsoid
    
    Parameters
    ----------
    src: str
        the input dataset, e.g. a VRT from function :func:`dem_autoload`
    dst: str
        the output dataset
    t_srs: None, int, str or osr.SpatialReference
        A target geographic reference system in WKT, EPSG, PROJ4 or OPENGIS format.
        See function :func:`spatialist.auxil.crsConvert()` for details.
        Default (None): use the crs of ``src``.
    tr: None or tuple
        the target resolution as (xres, yres)
    resampling_method: str
        the gdalwarp resampling method; See `here <https://gdal.org/programs/gdalwarp.html#cmdoption-gdalwarp-r>`_
        for options.
    geoid_convert: bool
        convert geoid heights?
    geoid: str
        the geoid model to be corrected, only used if ``geoid_convert == True``; current options:
         * 'EGM96'

    Returns
    -------

    """
    
    with Raster(src) as ras:
        nodata = ras.nodata
        epsg_in = ras.epsg
    
    if t_srs is None:
        epsg_out = epsg_in
    else:
        epsg_out = crsConvert(t_srs, 'epsg')
    
    gdalwarp_args = {'format': 'GTiff', 'multithread': True,
                     'srcNodata': nodata, 'dstNodata': nodata,
                     'srcSRS': 'EPSG:{}'.format(epsg_in),
                     'dstSRS': 'EPSG:{}'.format(epsg_out),
                     'resampleAlg': resampling_method}
    
    if tr is not None:
        gdalwarp_args.update({'xRes': tr[0],
                              'yRes': tr[1]})
    
    if geoid_convert:
        if gdal.__version__ < '2.2':
            raise RuntimeError('geoid conversion requires GDAL >= 2.2;'
                               'see documentation of gdalwarp')
        if geoid == 'EGM96':
            gdalwarp_args['srcSRS'] += '+5773'
        else:
            raise RuntimeError('geoid model not yet supported')
    
    try:
        message = 'creating mosaic'
        crs = gdalwarp_args['dstSRS']
        if crs != 'EPSG:4326':
            message += ' and reprojecting to {}'.format(crs)
        print(message)
        gdalwarp(src, dst, gdalwarp_args)
    except RuntimeError as e:
        if os.path.isfile(dst):
            os.remove(dst)
        errstr = str(e)
        if 'Cannot open egm96_15.gtx' in errstr:
            addition = '\nplease refer to the following site for instructions ' \
                       'on how to use the file egm96_15.gtx (requires proj.4 >= 5.0.0):\n' \
                       'https://gis.stackexchange.com/questions/258532/' \
                       'noaa-vdatum-gdal-variable-paths-for-linux-ubuntu'
            raise RuntimeError(errstr + addition)
        else:
            raise e


class DEMHandler:
    """
    | An interface to obtain DEM data for selected geometries
    | The files are downloaded into the ESA SNAP auxdata directory structure
    
    Parameters
    ----------
    geometries: list of spatialist.vector.Vector
        a list of geometries
    """
    
    def __init__(self, geometries):
        if not isinstance(geometries, list):
            raise RuntimeError('geometries must be of type list')
        
        for geometry in geometries:
            if geometry.getProjection('epsg') != 4326:
                raise RuntimeError('input geometry CRS must be WGS84 LatLon (EPSG 4326)')
        
        self.geometries = geometries
        try:
            self.auxdatapath = ExamineSnap().auxdatapath
        except AttributeError:
            self.auxdatapath = os.path.join(os.path.expanduser('~'), '.snap', 'auxdata')
    
    def __enter__(self):
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        return
    
    @staticmethod
    def __applybuffer(extent, buffer):
        ext = dict(extent)
        if buffer is not None:
            ext['xmin'] -= buffer
            ext['xmax'] += buffer
            ext['ymin'] -= buffer
            ext['ymax'] += buffer
        return ext
    
    @staticmethod
    def __buildvrt(archives, vrtfile, pattern, vsi, extent, nodata=None, srs=None):
        locals = [vsi + x for x in dissolve([finder(x, [pattern]) for x in archives])]
        if nodata is None:
            with Raster(locals[0]) as ras:
                nodata = ras.nodata
        opts = {'outputBounds': (extent['xmin'], extent['ymin'],
                                 extent['xmax'], extent['ymax']),
                'srcNodata': nodata}
        if srs is not None:
            opts['outputSRS'] = crsConvert(srs, 'wkt')
        gdalbuildvrt(src=locals, dst=vrtfile,
                     options=opts)
    
    def __commonextent(self, buffer=None):
        ext_new = {}
        for geo in self.geometries:
            if len(ext_new.keys()) == 0:
                ext_new = geo.extent
            else:
                for key in ['xmin', 'ymin']:
                    if geo.extent[key] > ext_new[key]:
                        ext_new[key] = geo.extent[key]
                for key in ['xmax', 'ymax']:
                    if geo.extent[key] < ext_new[key]:
                        ext_new[key] = geo.extent[key]
        ext_new = self.__applybuffer(ext_new, buffer)
        return ext_new
    
    @staticmethod
    def __retrieve(url, filenames, outdir):
        files = list(set(filenames))
        if not os.path.isdir(outdir):
            os.makedirs(outdir)
        locals = []
        for file in files:
            infile = '{}/{}'.format(url, file)
            outfile = os.path.join(outdir, os.path.basename(file))
            if not os.path.isfile(outfile):
                try:
                    input = urlopen(infile)
                    print('{} <<-- {}'.format(outfile, infile))
                except HTTPError:
                    continue
                with open(outfile, 'wb') as output:
                    output.write(input.read())
                input.close()
            if os.path.isfile(outfile):
                locals.append(outfile)
        return sorted(locals)
    
    @staticmethod
    def __retrieve_ftp(url, filenames, outdir, username, password):
        files = list(set(filenames))
        if not os.path.isdir(outdir):
            os.makedirs(outdir)
        pattern = r'(ftp(?:es|))://([a-z0-9.\-]*)[/]*((?:[a-zA-Z0-9/_]*|))'
        protocol, url, path = re.search(pattern, url).groups()
        if protocol == 'ftpes':
            ftp = ftplib.FTP_TLS(url)
            try:
                ftp.login(username, password)  # login anonymously before securing control channel
            except ftplib.error_perm as e:
                raise RuntimeError(str(e))
            ftp.prot_p()  # switch to secure data connection.. IMPORTANT! Otherwise, only the user and password is encrypted and not all the file data.
        else:
            ftp = ftplib.FTP(url, timeout=100)
            ftp.login()
        if path != '':
            ftp.cwd(path)
        locals = []
        for product_remote in files:
            product_local = os.path.join(outdir, os.path.basename(product_remote))
            if not os.path.isfile(product_local):
                try:
                    targetlist = ftp.nlst(product_remote)
                except ftplib.error_temp:
                    continue
                address = '{}://{}/{}{}'.format(protocol, url, path + '/' if path != '' else '', product_remote)
                print('{} <<-- {}'.format(product_local, address))
                with open(product_local, 'wb') as myfile:
                    ftp.retrbinary('RETR {}'.format(product_remote), myfile.write)
            if os.path.isfile(product_local):
                locals.append(product_local)
        ftp.close()
        return sorted(locals)
    
    @property
    def config(self):
        return {
            'AW3D30': {'url': 'ftp://ftp.eorc.jaxa.jp/pub/ALOS/ext1/AW3D30/release_v1804',
                       'nodata': -9999,
                       'vsi': '/vsitar/',
                       'pattern': {'dem': '*DSM.tif',
                                   'msk': '*MSK.tif',
                                   'stk': '*STK.tif'}
                       },
            'SRTM 1Sec HGT': {'url': 'https://step.esa.int/auxdata/dem/SRTMGL1',
                              'nodata': -32768.0,
                              'vsi': '/vsizip/',
                              'pattern': {'dem': '*.hgt'}
                              },
            'SRTM 3Sec': {'url': 'https://srtm.csi.cgiar.org/wp-content/uploads/files/srtm_5x5/TIFF',
                          'nodata': -32768.0,
                          'vsi': '/vsizip/',
                          'pattern': {'dem': 'srtm*.tif'}
                          },
            'TDX90m': {'url': 'ftpes://tandemx-90m.dlr.de',
                       'nodata': -32767.0,
                       'vsi': '/vsizip/',
                       'pattern': {'dem': '*_DEM.tif',
                                   'am2': '*_AM2.tif',
                                   'amp': '*_AMP.tif',
                                   'com': '*_COM.tif',
                                   'cov': '*_COV.tif',
                                   'hem': '*_HEM.tif',
                                   'lsm': '*_LSM.tif',
                                   'wam': '*_WAM.tif'}
                       }
        }
    
    def load(self, demType, vrt=None, buffer=None, username=None, password=None, product='dem'):
        """
        obtain DEM tiles for the given geometries
        
        Parameters
        ----------
        demType: str
            the type fo DEM to be used
        vrt: str or None
            an optional GDAL VRT file created from the obtained DEM tiles
        buffer: int or float
            a buffer in degrees to add around the individual geometries
        username: str
            the download account user name
        password: str
            the download account password
        product: str
            the sub-product to extract from the DEM product
             * 'AW3D30'
             
              - 'dem': the actual Digital Elevation Model
              - 'msk': mask information for each pixel (Cloud/Snow Mask, Land water and
                low correlation mask, Sea mask, Information of elevation dataset used
                for the void-filling processing)
              - 'stk': number of DSM-scene files which were used to produce the 5m resolution DSM
              
             * 'SRTM 1Sec HGT'
             
              - 'dem': the actual Digital Elevation Model
              
             * 'SRTM 3Sec'
             
              - 'dem': the actual Digital Elevation Model
              
             * 'TDX90m'
             
              - 'dem': the actual Digital Elevation Model
              - 'am2': Amplitude Mosaic representing the minimum value
              - 'amp': Amplitude Mosaic representing the mean value
              - 'com': Consistency Mask
              - 'cov': Coverage Map
              - 'hem': Height Error Map
              - 'lsm': Layover and Shadow Mask, based on SRTM C-band and Globe DEM data
              - 'wam': Water Indication Mask

        Returns
        -------
        list or None
            the names of the obtained files or None if a VRT file was defined
        """
        keys = self.config.keys()
        if demType not in keys:
            raise RuntimeError("demType '{}' is not supported\n  "
                               "possible options: '{}'"
                               .format(demType, "', '".join(keys)))
        
        products = self.config[demType]['pattern'].keys()
        if product not in products:
            raise RuntimeError("product '{0}' not available for demType '{1}'\n"
                               "  options: '{2}'".format(product, demType, "', '".join(products)))
        
        outdir = os.path.join(self.auxdatapath, 'dem', demType)
        
        remotes = []
        for geo in self.geometries:
            corners = self.__applybuffer(geo.extent, buffer)
            remotes.extend(self.remote_ids(corners, demType=demType))
        
        if demType in ['AW3D30', 'TDX90m']:
            locals = self.__retrieve_ftp(self.config[demType]['url'], remotes, outdir,
                                         username=username, password=password)
        else:
            locals = self.__retrieve(self.config[demType]['url'], remotes, outdir)
        
        if product == 'dem':
            nodata = self.config[demType]['nodata']
        else:
            nodata = 0
        
        if vrt is not None:
            self.__buildvrt(archives=locals, vrtfile=vrt,
                            pattern=self.config[demType]['pattern'][product],
                            vsi=self.config[demType]['vsi'],
                            extent=self.__commonextent(buffer),
                            nodata=nodata)
            return None
        return locals
    
    @staticmethod
    def remote_ids(extent, demType):
        """
        parse the names of the remote files overlapping with an area of interest
        
        Parameters
        ----------
        extent: dict
            the extent of the area of interest with keys xmin, xmax, ymin, ymax
        demType: str
            the type fo DEM to be used
        
        Returns
        -------
        str
            the sorted names of the remote files
        """
        
        # generate sequence of integer coordinates marking the tie points of the individual tiles
        def intrange(extent, step):
            lat = range(int(float(extent['ymin']) // step) * step,
                        (int(float(extent['ymax']) // step) + 1) * step,
                        step)
            lon = range(int(float(extent['xmin']) // step) * step,
                        (int(float(extent['xmax']) // step) + 1) * step,
                        step)
            return lat, lon
        
        def index(x=None, y=None, nx=3, ny=3):
            if x is not None:
                xf = '{ew}{x:0{nx}d}'.format(ew='W' if x < 0 else 'E', x=abs(x), nx=nx)
            else:
                xf = ''
            if y is not None:
                yf = '{ns}{y:0{ny}d}'.format(ns='S' if y < 0 else 'N', y=abs(y), ny=ny)
            else:
                yf = ''
            out = yf + xf
            return out
        
        if demType in ['SRTM 1Sec HGT', 'TDX90m']:
            lat, lon = intrange(extent, step=1)
            
            if demType == 'SRTM 1Sec HGT':
                remotes = ['{}.SRTMGL1.hgt.zip'.format(index(x, y, nx=3, ny=2))
                           for x in lon for y in lat]
            else:
                remotes = []
                for x in lon:
                    xr = abs(x) // 10 * 10
                    for y in lat:
                        xf = index(x=x, nx=3)
                        yf = index(y=y, ny=2)
                        remotes.append('90mdem/DEM/{y}/{hem}{xr:03d}/TDM1_DEM__30_{y}{x}.zip'
                                       .format(x=xf, xr=xr, y=yf, hem=xf[0]))
        elif demType == 'AW3D30':
            remotes = []
            lat, lon = intrange(extent, step=1)
            for x in lon:
                for y in lat:
                    remotes.append(
                        '{}/{}.tar.gz'.format(index(x // 5 * 5, y // 5 * 5),
                                              index(x, y)))
        
        elif demType == 'SRTM 3Sec':
            lat = range(int((60 - float(extent['ymin'])) // 5) + 1,
                        int((60 - float(extent['ymax'])) // 5) + 2)
            lon = range(int((float(extent['xmin']) + 180) // 5) + 1,
                        int((float(extent['xmax']) + 180) // 5) + 2)
            
            remotes = ['srtm_{:02d}_{:02d}.zip'.format(x, y) for x in lon for y in lat]
        else:
            raise ValueError('unknown demType: {}'.format(demType))
        
        return sorted(remotes)


def getAuxdata(datasets, scenes):
    def getOrbitContentVersions(contentVersion):
        content = contentVersion.read().split('\n')
        items = [re.split(r'\s*=\s*', x.strip('\r')) for x in content if re.search('^[0-9]{4}', x)]
        return dict(items)
    
    auxDataPath = os.path.join(expanduser("~"), '.snap/auxdata')
    
    scenes = [identify(scene) if isinstance(scene, str) else scene for scene in scenes]
    sensors = list(set([scene.sensor for scene in scenes]))
    for dataset in datasets:
        if dataset == 'SRTM 1Sec HGT':
            files = [x.replace('hgt', 'SRTMGL1.hgt.zip') for x in
                     list(set(dissolve([scene.getHGT() for scene in scenes])))]
            for file in files:
                infile = os.path.join('https://step.esa.int/auxdata/dem/SRTMGL1', file)
                outfile = os.path.join(auxDataPath, 'dem/SRTM 1Sec HGT', file)
                if not os.path.isfile(outfile):
                    print(infile)
                    try:
                        input = urlopen(infile)
                    except HTTPError:
                        print('-> not available')
                        continue
                    with open(outfile, 'wb') as output:
                        output.write(input.read())
                    input.close()
        elif dataset == 'POEORB':
            for sensor in sensors:
                if re.search('S1[AB]', sensor):
                    
                    dates = [(scene.start[:4], scene.start[4:6]) for scene in scenes]
                    years = list(set([x[0] for x in dates]))
                    
                    remote_contentVersion = urlopen(
                        'https://step.esa.int/auxdata/orbits/Sentinel-1/POEORB/remote_contentVersion.txt')
                    versions_remote = getOrbitContentVersions(remote_contentVersion)
                    
                    for year in years:
                        dir_orb = os.path.join(auxDataPath, 'Orbits/Sentinel-1/POEORB', year)
                        
                        if not os.path.isdir(dir_orb):
                            os.makedirs(dir_orb)
                        contentVersionFile = os.path.join(dir_orb, 'contentVersion.txt')
                        
                        if os.path.isfile(contentVersionFile):
                            contentVersion = open(contentVersionFile, 'r+')
                            versions_local = getOrbitContentVersions(contentVersion)
                        else:
                            contentVersion = open(contentVersionFile, 'w')
                            versions_local = {}
                        
                        combine = dict(set(versions_local.items()) & set(versions_remote.items()))
                        
                        dates_select = [x for x in dates if x[0] == year]
                        months = list(set([x[1] for x in dates_select]))
                        
                        orb_ids = sorted(
                            [x for x in ['{}-{}.zip'.format(year, month) for month in months] if not x in combine])
                        
                        if len(orb_ids) > 0:
                            contentVersion.write('#\n#{}\n'.format(strftime('%a %b %d %H:%M:%S %Z %Y', gmtime())))
                            
                            for orb_id in orb_ids:
                                orb_remote = urlopen(
                                    'https://step.esa.int/auxdata/orbits/Sentinel-1/POEORB/{}'.format(orb_id))
                                orb_remote_stream = zf.ZipFile(StringIO(orb_remote.read()), 'r')
                                orb_remote.close()
                                
                                targets = [x for x in orb_remote_stream.namelist() if
                                           not os.path.isfile(os.path.join(dir_orb, x))]
                                orb_remote_stream.extractall(dir_orb, targets)
                                orb_remote_stream.close()
                                
                                versions_local[orb_id] = versions_remote[orb_id]
                                
                                for key, val in versions_local.iteritems():
                                    contentVersion.write('{}={}\n'.format(key, val))
                        
                        contentVersion.close()
                    remote_contentVersion.close()
                else:
                    print('not implemented yet')
        elif dataset == 'Delft Precise Orbits':
            path_server = 'dutlru2.lr.tudelft.nl'
            subdirs = {'ASAR:': 'ODR.ENVISAT1/eigen-cg03c', 'ERS1': 'ODR.ERS-1/dgm-e04', 'ERS2': 'ODR.ERS-2/dgm-e04'}
            ftp = FTP(path_server)
            ftp.login()
            for sensor in sensors:
                if sensor in subdirs.keys():
                    path_target = os.path.join('pub/orbits', subdirs[sensor])
                    path_local = os.path.join(auxDataPath, 'Orbits/Delft Precise Orbits', subdirs[sensor])
                    ftp.cwd(path_target)
                    for item in ftp.nlst():
                        ftp.retrbinary('RETR ' + item, open(os.path.join(path_local, item), 'wb').write)
            ftp.quit()
        else:
            print('not implemented yet')
