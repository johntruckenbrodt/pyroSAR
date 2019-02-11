import os
import sys
import ftplib

if sys.version_info >= (3, 0):
    from urllib.request import urlopen
    from urllib.error import HTTPError
else:
    from urllib2 import urlopen, HTTPError

from .snap import ExamineSnap
from spatialist.ancillary import dissolve, finder
from spatialist.auxil import gdalbuildvrt


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

        - 'SRTM 3sec'

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
    
    >>> from pyroSAR import identify
    >>> from pyroSAR.auxdata import dem_autoload
    >>> from spatialist import gdalwarp
    # identify the SAR scene
    >>> scene = identify('S1A_IW_GRDH_1SDV_20141012T162337_20141012T162402_002799_00326F_8197.zip')
    # extract the bounding box as spatialist.Vector object
    >>> bbox = scene.bbox()
    # download the tiles and virtually combine them in an in-memory VRT file subsetted to the extent of the SAR scene
    >>> vrt = dem_autoload(geometries=[bbox], demType='SRTM 1Sec HGT', vrt='/vsimem/srtm1.vrt', buffer=0.01)
    # write the final GeoTiff file
    >>> outname = scene.outname_base() + 'srtm1.tif'
    >>> gdalwarp(src=vrt, dst='', options={'format': 'GTiff'})
    """
    with DEMHandler(geometries) as handler:
        if demType == 'AW3D30':
            return handler.aw3d30(vrt, buffer)
        elif demType == 'SRTM 1Sec HGT':
            return handler.srtm_1sec_hgt(vrt, buffer)
        elif demType == 'SRTM 3Sec':
            return handler.srtm_3sec(vrt, buffer)
        elif demType == 'TDX90m':
            return handler.tdx90m(username=username,
                                  password=password,
                                  vrt=vrt,
                                  buffer=buffer)
        else:
            raise RuntimeError('demType unknown: {}'.format(demType))


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
    def __buildvrt(archives, vrtfile, pattern, vsi, extent, nodata):
        locals = [vsi + x for x in dissolve([finder(x, [pattern]) for x in archives])]
        gdalbuildvrt(src=locals, dst=vrtfile,
                     options={'outputBounds': (extent['xmin'], extent['ymin'],
                                               extent['xmax'], extent['ymax']),
                              'srcNodata': nodata})
    
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
        return locals
    
    @staticmethod
    def __retrieve_ftp(url, filenames, outdir, username, password):
        files = list(set(filenames))
        if not os.path.isdir(outdir):
            os.makedirs(outdir)
        ftps = ftplib.FTP_TLS(url)
        ftps.login(username, password)  # login anonymously before securing control channel
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
        return locals
    
    def aw3d30(self, vrt=None, buffer=None):
        """
        obtain ALOS Global Digital Surface Model "ALOS World 3D - 30m (AW3D30)"
        
        Parameters
        ----------
        vrt: str or None
            an optional GDAL VRT file created from the obtained DEM tiles
        buffer: int or float
            a buffer in degrees to add around the individual geometries
        
        Returns
        -------
        list or str
            the names of the obtained files or the name of the VRT file
        """
        url = 'ftp://ftp.eorc.jaxa.jp/pub/ALOS/ext1/AW3D30/release_v1804'
        outdir = os.path.join(self.auxdatapath, 'dem', 'AW3D30')
        
        def index(x, y):
            return '{}{:03d}{}{:03d}'.format('S' if y < 0 else 'N', abs(y),
                                             'W' if x < 0 else 'E', abs(x))
        
        files = []
        for geo in self.geometries:
            corners = self.__applybuffer(geo.extent, buffer)
            xmin = int(corners['xmin'] // 5)
            xmax = int(corners['xmax'] // 5)
            ymin = int(corners['ymin'] // 5)
            ymax = int(corners['ymax'] // 5)
            
            for x in range(xmin, xmax + 1):
                for y in range(ymin, ymax + 1):
                    files.append('{}_{}.tar.gz'.format(index(x * 5, y * 5),
                                                       index(x * 5 + 5, y * 5 + 5)))
        
        locals = self.__retrieve(url, files, outdir)
        if vrt is not None:
            self.__buildvrt(archives=locals, vrtfile=vrt, pattern='*DSM.tif',
                            vsi='/vsitar/', extent=self.__commonextent(buffer),
                            nodata=-9999)
            return vrt
        return locals
    
    def srtm_1sec_hgt(self, vrt=None, buffer=None):
        """
        obtain SRTM 1arcsec DEM tiles in HGT format
        
        Parameters
        ----------
        vrt: str or None
            an optional GDAL VRT file created from the obtained DEM tiles
        buffer: int or float
            a buffer in degrees to add around the individual geometries
        
        Returns
        -------
        list or str
            the names of the obtained files or the name of the VRT file
        """
        url = 'https://step.esa.int/auxdata/dem/SRTMGL1'
        outdir = os.path.join(self.auxdatapath, 'dem', 'SRTM 1Sec HGT')
        
        files = []
        
        for geo in self.geometries:
            corners = self.__applybuffer(geo.extent, buffer)
            
            # generate sequence of integer coordinates marking the tie points of the overlapping hgt tiles
            lat = range(int(float(corners['ymin']) // 1),
                        int(float(corners['ymax']) // 1) + 1)
            lon = range(int(float(corners['xmin']) // 1),
                        int(float(corners['xmax']) // 1) + 1)
            
            # convert coordinates to string with leading zeros and hemisphere identification letter
            lat = [str(x).zfill(2 + len(str(x)) - len(str(x).strip('-'))) for x in lat]
            lat = [x.replace('-', 'S') if '-' in x else 'N' + x for x in lat]
            
            lon = [str(x).zfill(3 + len(str(x)) - len(str(x).strip('-'))) for x in lon]
            lon = [x.replace('-', 'W') if '-' in x else 'E' + x for x in lon]
            
            # concatenate all formatted latitudes and longitudes with each other as final product
            files.extend([x + y + '.SRTMGL1.hgt.zip' for x in lat for y in lon])
        
        locals = self.__retrieve(url, files, outdir)
        
        if vrt is not None:
            self.__buildvrt(archives=locals, vrtfile=vrt, pattern='*.hgt',
                            vsi='/vsizip/', extent=self.__commonextent(buffer),
                            nodata=-32768.0)
            return vrt
        return locals
    
    def srtm_3sec(self, vrt=None, buffer=None):
        """
        obtain SRTM 3arcsec DEM tiles in GeoTiff format
        
        Parameters
        ----------
        vrt: str or None
            an optional GDAL VRT file created from the obtained DEM tiles
        buffer: int or float
            a buffer in degrees to add around the individual geometries
        
        Returns
        -------
        list or str
            the names of the obtained files or the name of the VRT file
        """
        url = 'http://srtm.csi.cgiar.org/wp-content/uploads/files/srtm_5x5/TIFF'
        outdir = os.path.join(self.auxdatapath, 'dem', 'SRTM 3Sec')
        files = []
        for geo in self.geometries:
            corners = self.__applybuffer(geo.extent, buffer)
            x_id = [int((corners[x] + 180) // 5) + 1 for x in ['xmin', 'xmax']]
            y_id = [int((60 - corners[x]) // 5) + 1 for x in ['ymin', 'ymax']]
            files.extend(['srtm_{:02d}_{:02d}.zip'.format(x, y) for x in x_id for y in y_id])
        locals = self.__retrieve(url, files, outdir)
        if vrt is not None:
            self.__buildvrt(archives=locals, vrtfile=vrt, pattern='*.tif',
                            vsi='/vsizip/', extent=self.__commonextent(buffer),
                            nodata=-32768.0)
            return vrt
        return locals
    
    def tdx90m(self, username, password, vrt=None, buffer=None):
        """
        | obtain TanDEM-X 90 m DEM tiles in GeoTiff format
        | registration via DLR is necessary, see https://geoservice.dlr.de/web/dataguide/tdm90
        
        Parameters
        ----------
        username: str
            the DLR user name
        password: str
            the user password
        vrt: str or None
            an optional GDAL VRT file created from the obtained DEM tiles
        buffer: int or float
            a buffer in degrees to add around the individual geometries

        Returns
        -------
        list or str
            the names of the obtained files or the name of the VRT file
        """
        url = 'tandemx-90m.dlr.de'
        outdir = os.path.join(self.auxdatapath, 'dem', 'TDX90m')
        
        remotes = []
        for geo in self.geometries:
            corners = self.__applybuffer(geo.extent, buffer)
            lat = range(int(float(corners['ymin']) // 1),
                        int(float(corners['ymax']) // 1) + 1)
            lon = range(int(float(corners['xmin']) // 1),
                        int(float(corners['xmax']) // 1) + 1)
            
            # convert coordinates to string with leading zeros and hemisphere identification letter
            lat_id = [str(x).zfill(2 + len(str(x)) - len(str(x).strip('-'))) for x in lat]
            lat_id = [x.replace('-', 'S') if '-' in x else 'N' + x for x in lat_id]
            
            lon_id = [str(x).zfill(3 + len(str(x)) - len(str(x).strip('-'))) for x in lon]
            lon_id = [x.replace('-', 'W') if '-' in x else 'E' + x for x in lon_id]
            
            for x in lon_id:
                for y in lat_id:
                    xr = int(x[1:]) // 10 * 10
                    remotes.append('90mdem/DEM/{y}/{hem}{xr:03d}/TDM1_DEM__30_{y}{x}.zip'
                                   .format(x=x, xr=xr, y=y, hem=x[0]))
        
        locals = self.__retrieve_ftp(url, remotes, outdir, username=username, password=password)
        if vrt is not None:
            self.__buildvrt(archives=locals, vrtfile=vrt, pattern='*_DEM.tif',
                            vsi='/vsizip/', extent=self.__commonextent(buffer),
                            nodata=-32767.0)
            return vrt
        return locals
