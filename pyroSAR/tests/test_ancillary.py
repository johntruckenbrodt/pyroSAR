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
    assert anc.parse_literal(['1', '2.2', 'a']) == [1, 2.2, 'a']
    with pytest.raises(IOError):
        anc.parse_literal(1)


def test_seconds():
    assert anc.seconds('test_20151212T234411') == 3658952651.0


def test_run(tmpdir, testdata):
    log = os.path.join(str(tmpdir), 'test_run.log')
    out, err = anc.run(cmd=['gdalinfo', testdata['tif']],
                       logfile=log, void=False)
    with pytest.raises(OSError):
        anc.run(['foobar'])
    with pytest.raises(sp.CalledProcessError):
        anc.run(['gdalinfo', 'foobar'])


def test_which():
    program = anc.which('gdalinfo')
    assert os.path.isfile(program)
    assert anc.which(program) == program
    assert anc.which('foobar') is None


def test_multicore():
    add = lambda x, y, z: x + y + z
    assert anc.multicore(add, cores=2, multiargs={'x': [1, 2]}, y=5, z=9) == [15, 16]
    assert anc.multicore(add, cores=2, multiargs={'x': [1, 2], 'y': [5, 6]}, z=9) == [15, 17]
    with pytest.raises(AttributeError):
        anc.multicore(add, cores=2, multiargs={'foobar': [1, 2]}, y=5, z=9)
    with pytest.raises(AttributeError):
        anc.multicore(add, cores=2, multiargs={'x': [1, 2]}, y=5, foobar=9)
    with pytest.raises(AttributeError):
        anc.multicore(add, cores=2, multiargs={'x': [1, 2], 'y': [5, 6, 7]}, foobar=9)


def test_finder(tmpdir):
    dir = str(tmpdir)
    dir_sub1 = os.path.join(dir, 'testdir1')
    dir_sub2 = os.path.join(dir, 'testdir2')
    os.makedirs(dir_sub1)
    os.makedirs(dir_sub2)
    with open(os.path.join(dir_sub1, 'testfile1.txt'), 'w') as t1:
        t1.write('test')
    with open(os.path.join(dir_sub2, 'testfile2.txt'), 'w') as t2:
        t2.write('test')
    assert len(anc.finder(dir, ['test*'], foldermode=0)) == 2
    assert len(anc.finder(dir, ['test*'], foldermode=0, recursive=False)) == 0
    assert len(anc.finder(dir, ['test*'], foldermode=1)) == 4
    assert len(anc.finder(dir, ['test*'], foldermode=2)) == 2
    assert len(anc.finder([dir_sub1, dir_sub2], ['test*'])) == 2
    with pytest.raises(TypeError):
        anc.finder(1, [])


def test_rescale():
    assert anc.rescale([1000, 2000, 3000], [1, 3]) == [1, 2, 3]
    with pytest.raises(RuntimeError):
        anc.rescale([1000, 1000])


def test_Queue():
    st = anc.Queue()
    st.push('a')
    assert st.pop() == 'a'
    assert st.length() == 0


def test_Stack():
    st = anc.Stack()
    assert st.empty() is True
    st = anc.Stack(['a', 'b'])
    assert st.length() == 2
    st = anc.Stack('a')
    st.push('b')
    st.push(['c', 'd'])
    assert st.pop() == 'd'
    st.flush()
    assert st.empty() is True
