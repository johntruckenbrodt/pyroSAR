import os
import sys
import time
import pytest
from pyroSAR import identify
from pyroSAR.S1 import OSV


def test_scene_osv(tmpdir, testdata):
    id = identify(testdata['s1'])
    osvdir = os.path.join(str(tmpdir), 'osv')
    if sys.version_info >= (2, 7, 9):
        id.getOSV(osvdir)
        with OSV(osvdir) as osv:
            with pytest.raises(IOError):
                osv.catch(sensor='S1A', osvtype='XYZ')
            res = osv.catch(sensor='S1A', osvtype='RES', start=osv.mindate('POE'), stop=osv.maxdate('POE'))
            assert len(res) == 0
            
            assert len(osv.getLocals('POE')) == 3
            assert len(osv.getLocals('RES')) == 0
            assert osv.match(sensor=id.sensor, timestamp=id.start, osvtype='POE') is not None
            assert osv.match(sensor=id.sensor, timestamp=id.start, osvtype=['POE', 'RES']) is not None
            assert osv.match(sensor=id.sensor, timestamp=id.start, osvtype='RES') is None
            for item in osv.getLocals('POE')[1:3]:
                os.remove(item)
            assert len(osv.getLocals('POE')) == 1
            res = osv.catch(sensor='S1A', osvtype='RES', start='20180101T120000', stop='20180102T120000')
            assert len(res) == 33
            osv.retrieve(res[0:3])
            assert len(osv.getLocals('RES')) == 3
            res = osv.catch(sensor='S1A', osvtype='POE', start=time.strftime('%Y%m%dT%H%M%S'))
    else:
        with pytest.raises(RuntimeError):
            id.getOSV(osvdir, osvType='POE')
