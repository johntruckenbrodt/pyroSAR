from pyroSAR._dev_config import Storage, LOOKUP, URL, STORAGE, ConfigHandler
import os
import pytest


class TestStorage:
    def test_insert(self):
        storage = Storage(a=1, b=2)
        assert storage.a == 1
        assert storage.b == 2
    
    def test_keys(self):
        storage = Storage(a=1, b=2)
        key = storage.keys()
        key = list(key)
        assert key[0] == 'a' or 'b'
        assert key[1] == 'b' or 'a'


class TestLookup:
    def test_attributes(self):
        assert LOOKUP.attributes['sensor'] == 'TEXT'
        assert LOOKUP.attributes['vh'] == 'INTEGER'


class TestSTORAGE:
    def test_STORAGE_URL(self):
        assert STORAGE.URL.dem.ace == URL.dem.ace
        assert STORAGE.URL.orbit.doris == URL.orbit.doris
        assert STORAGE.URL.auxcal.ers == URL.auxcal.ers
    
    def test_STORAGE_LOOKUP(self):
        assert LOOKUP.attributes['sensor'] == STORAGE.LOOKUP.attributes['sensor']
        assert LOOKUP.attributes['vh'] == STORAGE.LOOKUP.attributes['vh']


class TestConfigHandler:
    
    def test_make_dir_and_config(self, tmpdir):
        conf = ConfigHandler()
        
        path_pyrosar = os.path.exists(conf._ConfigHandler__GLOBAL['path'])
        path_config = os.path.isfile(conf._ConfigHandler__GLOBAL['config'])
        
        assert path_pyrosar is True
        assert path_config is True
    
    def test_add_section(self):
        conf = ConfigHandler()
        conf.add_section('FOO')
        
        assert 'FOO' in conf.sections
    
    def test_options(self, tmpdir):
        conf = ConfigHandler()
        conf.set('FOO', 'bar', 'foobar')
        
        # cannot set attribute for section that does not exist
        with pytest.raises(AttributeError):
            conf.set('SNAPp', 'etc', 'temp/dir')
        
        assert conf['FOO']['bar'] == 'foobar'
        assert conf['FOO'] == {'bar': 'foobar'}
    
    def test_overwrite(self, tmpdir):
        conf = ConfigHandler()
        
        with pytest.raises(RuntimeError):
            conf.set('FOO', 'bar', 'loremipsum')
        
        conf.set('FOO', 'bar', 'loremipsum', overwrite=True)
        assert conf['FOO']['bar'] == 'loremipsum'
    
    def test_remove(self, tmpdir):
        conf = ConfigHandler()
        
        with pytest.raises(AttributeError):
            conf.remove_option('SNAP', 'kex')
        
        with pytest.raises(AttributeError):
            conf.remove_option('SNApP', 'etc')
        
        conf.remove_option('FOO', 'bar')
        assert list(conf['FOO'].keys()) == []
        
        conf.remove_section('FOO')
