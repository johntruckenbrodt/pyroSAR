import os
import pytest
import subprocess as sp
import pyroSAR.ancillary as anc


def test_dissolve_with_lists():
    assert anc.dissolve([[1, 2], [3, 4]]) == [1, 2, 3, 4]
    assert anc.dissolve([[[1]]]) == [1]
    assert anc.dissolve(((1, 2,), (3, 4))) == [1, 2, 3, 4]
    assert anc.dissolve(((1, 2), (1, 2))) == [1, 2, 1, 2]


def test_union():
    assert anc.union([1], [1]) == [1]


def test_dictmerge():
     assert anc.dictmerge({'a': 1, 'b': 2}, {'c': 3, 'd': 4}) == {'a': 1, 'b': 2, 'c': 3, 'd': 4}


def test_parse_literal():
    assert anc.parse_literal(['1', '2.2', 'a'])== [1, 2.2, 'a']
    with pytest.raises(IOError):
        anc.parse_literal(1)


def test_seconds():
    assert anc.seconds('test_20151212T234411') == 3658952651.0


def test_run(tmpdir):
    log = os.path.join(str(tmpdir), 'test_run.log')
    out, err = anc.run(cmd=['gdalinfo',
                            'pyroSAR/tests/data/S1A__IW___A_20150309T173017_VV_grd_mli_geo_norm_db.tif'],
                       logfile=log, void=False)
    with pytest.raises(OSError):
        anc.run(['foobar'])
    with pytest.raises(sp.CalledProcessError):
        anc.run(['gdalinfo', 'foobar'])


def test_which():
    assert os.path.isfile(anc.which('gdalinfo'))
