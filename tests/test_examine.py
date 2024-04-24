import os
import pytest
from pyroSAR.examine import ExamineSnap, SnapProperties


def test_snap_config(tmpdir):
    conf_snap = ExamineSnap()
    conf = SnapProperties(path=os.path.dirname(conf_snap.etc))
    tmp = conf.userpath
    conf.userpath = tmpdir
    assert conf.userpath == tmpdir
    with pytest.raises(KeyError):
        conf['foobar'] = tmpdir
    
    conf['snap.jai.tileCacheSize'] = 2048
    assert conf['snap.jai.tileCacheSize'] == 2048
    assert isinstance(conf['snap.jai.tileCacheSize'], int)
    
    conf = SnapProperties(path=os.path.dirname(conf_snap.etc))
    assert conf['snap.jai.tileCacheSize'] == 2048
    assert isinstance(conf['snap.jai.tileCacheSize'], int)
    
    conf['snap.jai.tileCacheSize'] = 2048.
    assert isinstance(conf['snap.jai.tileCacheSize'], float)
    
    conf = SnapProperties(path=os.path.dirname(conf_snap.etc))
    assert conf['snap.jai.tileCacheSize'] == 2048.
    assert isinstance(conf['snap.jai.tileCacheSize'], float)
    
    conf['snap.jai.tileCacheSize'] = None
    assert conf['snap.jai.tileCacheSize'] is None
    
    conf = SnapProperties(path=os.path.dirname(conf_snap.etc))
    assert conf['snap.jai.tileCacheSize'] is None
    
    conf['snap.jai.tileCacheSize'] = True
    assert conf['snap.jai.tileCacheSize'] is True

    conf = SnapProperties(path=os.path.dirname(conf_snap.etc))
    assert conf['snap.jai.tileCacheSize'] is True
    
    conf.userpath = tmp
