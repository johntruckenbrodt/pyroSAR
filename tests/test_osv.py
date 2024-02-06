import os
import time
import pytest
from pyroSAR import identify
from pyroSAR.S1 import OSV
from datetime import datetime, timedelta


def test_osv_cleanres(tmpdir):
    with OSV(str(tmpdir)) as osv:
        assert osv.getLocals('POE') == []
        assert osv.getLocals('RES') == []
        now = (datetime.now() - timedelta(hours=10)).strftime('%Y%m%dT%H%M%S')
        res = osv.catch(sensor='S1A', osvtype='RES', start=now)
        nfiles = len(res)
        osv.retrieve(res)
        osv.clean_res()
        assert len(osv.getLocals('RES')) == nfiles


def test_scene_osv(tmpdir, testdata):
    id = identify(testdata['s1_orbit'])
    osvdir = os.path.join(str(tmpdir), 'osv')
    id.getOSV(osvdir)
    with OSV(osvdir) as osv:
        with pytest.raises(RuntimeError):
            osv.catch(sensor='S1A', osvtype='XYZ')
        res = osv.catch(sensor='S1A', osvtype='RES', start=osv.mindate('POE'), stop=osv.maxdate('POE'))
        assert len(res) == 0
        
        assert len(osv.getLocals('POE')) == 1
        assert len(osv.getLocals('RES')) == 0
        assert osv.match(sensor=id.sensor, timestamp=id.start, osvtype='POE') is not None
        assert osv.match(sensor=id.sensor, timestamp=id.start, osvtype=['POE', 'RES']) is not None
        assert osv.match(sensor=id.sensor, timestamp=id.start, osvtype='RES') is None
        for item in osv.getLocals('POE')[1:3]:
            os.remove(item)
        assert len(osv.getLocals('POE')) == 1
        res = osv.catch(sensor='S1A', osvtype='RES', start='20210201T00000', stop='20210201T150000', url_option=1)
        assert len(res) == 11
        osv.retrieve(res[0:3])
        assert len(osv.getLocals('RES')) == 3
        # check retrieving files for the current day (e.g. to ensure that search is not extended to the future)
        poe = osv.catch(sensor='S1A', osvtype='POE', start=time.strftime('%Y%m%dT%H%M%S'))
        assert len(poe) == 0
        # check retrieving files whose start is in the previous month of the search start
        poe = osv.catch(sensor='S1A', osvtype='POE', start='20220201T163644', stop='20220201T163709')
        assert len(poe) == 1
