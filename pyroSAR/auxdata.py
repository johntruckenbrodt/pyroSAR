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
from .snap import ExamineSnap
from spatialist import Raster
from spatialist.ancillary import dissolve, finder
from spatialist.auxil import gdalbuildvrt, crsConvert, gdalwarp


def dem_autoload(geometries, demType, vrt=None, buffer=None, username=None, password=None):
    """
    obtain all relevant DEM tiles for selected geometries

    Parameters
    ----------
    geometries: list
        a list of :class:`spatialist.vector.Vector` geometries to obtain DEM data for; CRS must be WGS84 LatLon (EPSG 4326)
    demType: str
        the type fo DEM to be used; current options:

        - 'AW3D30' (ALOS Global Digital Surface Model "ALOS World 3D - 30m (AW3D30)")

          * url: ftp://ftp.eorc.jaxa.jp/pub/ALOS/ext1/AW3D30/release_v1804

        - 'SRTM 1Sec HGT'

          * url: https://step.esa.int/auxdata/dem/SRTMGL1

        - 'SRTM 3Sec'

          * url: http://srtm.csi.cgiar.org/wp-content/uploads/files/srtm_5x5/TIFF
        
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

    Returns
    -------
    list or str
        the names of the obtained files or the name of the VRT file
    
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
        vrt = dem_autoload(geometries=[bbox], demType='SRTM 1Sec HGT',
                           vrt='/vsimem/srtm1.vrt', buffer=0.01)
        
        # write the final GeoTiff file
        outname = scene.outname_base() + 'srtm1.tif'
        gdalwarp(src=vrt, dst=outname, options={'format': 'GTiff'})
        
        # alternatively use function dem_create and warp the DEM to UTM
        # including conversion from geoid to ellipsoid heights
        from pyroSAR.auxdata import dem_create
        outname = scene.outname_base() + 'srtm1_ellp.tif'
        dem_create(src=vrt, dst=outname, t_srs=32632, tr=(30, 30), geoid_convert=True, geoid='EGM96')
    """
    with DEMHandler(geometries) as handler:
        return handler.load(demType=demType,
                            username=username,
                            password=password,
                            vrt=vrt,
                            buffer=buffer)


def dem_create(src, dst, t_srs=None, tr=None, geoid_convert=False, geoid='EGM96'):
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
                     'dstSRS': 'EPSG:{}'.format(epsg_out)}
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
        gdalwarp(src, dst, gdalwarp_args)
    except RuntimeError as e:
        if os.path.isfile(dst):
            os.remove(dst)
        errstr = str(e)
        if 'Cannot open egm96_15.gtx' in errstr:
            addition = '\nplease refer to the following site for instructions ' \
                       'on how to use the file egm96_15.gtx:\n' \
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
    def remote_ids(extent, demType):
        if demType in ['SRTM 1Sec HGT', 'TDX90m']:
            # generate sequence of integer coordinates marking the tie points of the overlapping tiles
            lat = range(int(float(extent['ymin']) // 1),
                        int(float(extent['ymax']) // 1) + 1)
            lon = range(int(float(extent['xmin']) // 1),
                        int(float(extent['xmax']) // 1) + 1)
            
            # convert coordinates to string with leading zeros and hemisphere identification letter
            lat = [str(x).zfill(2 + len(str(x)) - len(str(x).strip('-'))) for x in lat]
            lat = [x.replace('-', 'S') if '-' in x else 'N' + x for x in lat]
            
            lon = [str(x).zfill(3 + len(str(x)) - len(str(x).strip('-'))) for x in lon]
            lon = [x.replace('-', 'W') if '-' in x else 'E' + x for x in lon]
            if demType == 'SRTM 1Sec HGT':
                remotes = [x + y + '.SRTMGL1.hgt.zip' for x in lat for y in lon]
            else:
                remotes = []
                for x in lon:
                    for y in lat:
                        xr = int(x[1:]) // 10 * 10
                        remotes.append('90mdem/DEM/{y}/{hem}{xr:03d}/TDM1_DEM__30_{y}{x}.zip'
                                       .format(x=x, xr=xr, y=y, hem=x[0]))
        elif demType == 'AW3D30':
            def index(x, y):
                return '{}{:03d}{}{:03d}'.format('S' if y < 0 else 'N', abs(y),
                                                 'W' if x < 0 else 'E', abs(x))
            
            remotes = []
            lat = range(int(float(extent['ymin']) // 5),
                        int(float(extent['ymax']) // 5) + 1)
            lon = range(int(float(extent['xmin']) // 5),
                        int(float(extent['xmax']) // 5) + 1)
            for x in lon:
                for y in lat:
                    remotes.append('{}_{}.tar.gz'.format(index(x * 5, y * 5),
                                                         index(x * 5 + 5, y * 5 + 5)))
        elif demType == 'SRTM 3Sec':
            lat = range(int((60 - float(extent['ymin'])) // 5) + 1,
                        int((60 - float(extent['ymax'])) // 5) + 2)
            lon = range(int((float(extent['xmin']) + 180) // 5) + 1,
                        int((float(extent['xmax']) + 180) // 5) + 2)
            
            remotes = ['srtm_{:02d}_{:02d}.zip'.format(x, y) for x in lon for y in lat]
        else:
            raise ValueError('unknown demType: {}'.format(demType))
        
        return sorted(remotes)
    
    def __commonextent(self, buffer=None):
        ext_new = {}
        for geo in self.geometries:
            if len(ext_new.keys()) == 0:
                ext_new = geo.extent
            else:
                for key in ['xmin', 'ymin']:
                    if geo.extent[key] < ext_new[key]:
                        ext_new[key] = geo.extent[key]
                for key in ['xmax', 'ymax']:
                    if geo.extent[key] > ext_new[key]:
                        ext_new[key] = geo.extent[key]
        ext_new = self.__applybuffer(ext_new, buffer)
        return ext_new
    
    @staticmethod
    def __buildvrt(archives, vrtfile, pattern, vsi, extent, nodata, srs=None):
        locals = [vsi + x for x in dissolve([finder(x, [pattern]) for x in archives])]
        opts = {'outputBounds': (extent['xmin'], extent['ymin'],
                                 extent['xmax'], extent['ymax']),
                'srcNodata': nodata}
        if srs is not None:
            opts['outputSRS'] = crsConvert(srs, 'wkt')
        gdalbuildvrt(src=locals, dst=vrtfile,
                     options=opts)
    
    @staticmethod
    def __retrieve(url, filenames, outdir):
        files = list(set(filenames))
        if not os.path.isdir(outdir):
            os.makedirs(outdir)
        locals = []
        for file in files:
            infile = '{}/{}'.format(url, file)
            outfile = os.path.join(outdir, file)
            if not os.path.isfile(outfile):
                try:
                    input = urlopen(infile)
                    print('{} -->> {}'.format(infile, outfile))
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
        ftps = ftplib.FTP_TLS(url)
        try:
            ftps.login(username, password)  # login anonymously before securing control channel
        except ftplib.error_perm as e:
            raise RuntimeError(str(e))
        ftps.prot_p()  # switch to secure data connection.. IMPORTANT! Otherwise, only the user and password is encrypted and not all the file data.
        
        locals = []
        for product_remote in files:
            product_local = os.path.join(outdir, os.path.basename(product_remote))
            if not os.path.isfile(product_local):
                try:
                    targetlist = ftps.nlst(product_remote)
                except ftplib.error_temp:
                    continue
                print('ftpes://{}/{} -->> {}'.format(url, product_remote, product_local))
                with open(product_local, 'wb') as myfile:
                    ftps.retrbinary('RETR {}'.format(product_remote), myfile.write)
            if os.path.isfile(product_local):
                locals.append(product_local)
        ftps.close()
        return sorted(locals)
    
    @property
    def config(self):
        return {
            'AW3D30': {'url': 'ftp://ftp.eorc.jaxa.jp/pub/ALOS/ext1/AW3D30/release_v1804',
                       'nodata': -9999,
                       'vsi': '/vsitar/',
                       'pattern': '*DSM.tif'},
            'SRTM 1Sec HGT': {'url': 'https://step.esa.int/auxdata/dem/SRTMGL1',
                              'nodata': -32768.0,
                              'vsi': '/vsizip/',
                              'pattern': '*.hgt'},
            'SRTM 3Sec': {'url': 'http://srtm.csi.cgiar.org/wp-content/uploads/files/srtm_5x5/tiff',
                          'nodata': -32768.0,
                          'vsi': '/vsizip/',
                          'pattern': '*.tif'},
            'TDX90m': {'url': 'tandemx-90m.dlr.de',
                       'nodata': -32767.0,
                       'vsi': '/vsizip/',
                       'pattern': '*_DEM.tif'}
        }
    
    def load(self, demType, vrt=None, buffer=None, username=None, password=None):
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

        Returns
        -------
        list or str
            the names of the obtained files or the name of the VRT file
        """
        keys = self.config.keys()
        if demType not in keys:
            raise RuntimeError("demType '{}' is not supported\n  "
                               "possible options: '{}'"
                               .format(demType, "', '".join(keys)))
        
        outdir = os.path.join(self.auxdatapath, 'dem', demType)
        
        remotes = []
        for geo in self.geometries:
            corners = self.__applybuffer(geo.extent, buffer)
            remotes.extend(self.remote_ids(corners, demType=demType))
        
        if demType == 'TDX90m':
            locals = self.__retrieve_ftp(self.config[demType]['url'], remotes, outdir,
                                         username=username, password=password)
        else:
            locals = self.__retrieve(self.config[demType]['url'], remotes, outdir)
        
        if vrt is not None:
            self.__buildvrt(archives=locals, vrtfile=vrt,
                            pattern=self.config[demType]['pattern'],
                            vsi=self.config[demType]['vsi'],
                            extent=self.__commonextent(buffer),
                            nodata=self.config[demType]['nodata'])
            return vrt
        return locals


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
