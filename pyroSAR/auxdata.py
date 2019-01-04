import os
import sys

if sys.version_info >= (3, 0):
    from urllib.request import urlopen
    from urllib.error import HTTPError
else:
    from urllib2 import urlopen, HTTPError

from . import identify
from .snap import ExamineSnap

from spatialist.ancillary import dissolve


class Handler:
    def __init__(self, scenes):
        self.scenes = [identify(scene) if isinstance(scene, str) else scene for scene in scenes]
        try:
            self.auxdatapath = ExamineSnap().auxdatapath
        except AttributeError:
            self.auxdatapath = os.path.join(os.path.expanduser('~'), '.snap', 'auxdata')
    
    @staticmethod
    def __retrieve(url, filenames, outdir):
        files = list(set(filenames))
        locals = []
        for file in files:
            infile = os.path.join(url, file)
            outfile = os.path.join(outdir, file)
            if not os.path.isfile(outfile):
                try:
                    input = urlopen(infile)
                    print('downloading file {}'.format(infile))
                except HTTPError:
                    continue
                with open(outfile, 'wb') as output:
                    output.write(input.read())
                input.close()
            if os.path.isfile(outfile):
                locals.append(outfile)
        return locals
    
    def srtm_1sec_hgt(self):
        url = 'https://step.esa.int/auxdata/dem/SRTMGL1'
        outdir = os.path.join(self.auxdatapath, 'dem', 'SRTM 1Sec HGT')
        files = [x.replace('hgt', 'SRTMGL1.hgt.zip') for x in
                 list(set(dissolve([scene.getHGT() for scene in self.scenes])))]
        locals = self.__retrieve(url, files, outdir)
        return locals
    
    def srtm_3sec(self):
        url = 'http://srtm.csi.cgiar.org/wp-content/uploads/files/srtm_5x5/TIFF'
        outdir = os.path.join(self.auxdatapath, 'dem', 'SRTM 3Sec')
        files = []
        for scene in self.scenes:
            corners = scene.getCorners()
            x_id = [int((corners[x] + 180) // 5) + 1 for x in ['xmin', 'xmax']]
            y_id = [int((60 - corners[x]) // 5) + 1 for x in ['ymin', 'ymax']]
            files.extend(['srtm_{:02d}_{:02d}.zip'.format(x, y) for x in x_id for y in y_id])
        locals = self.__retrieve(url, files, outdir)
        return locals
