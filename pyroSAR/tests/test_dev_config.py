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
    def test_suffix(self):
        assert LOOKUP.snap.suffix['Apply-Orbit-File'] == 'Orb'
        assert LOOKUP.snap.suffix['Terrain-Correction'] == 'TC'
    
    def test_attributes(self):
        assert LOOKUP.attributes['sensor'] == 'TEXT'
        assert LOOKUP.attributes['vh'] == 'INTEGER'


class TestSTORAGE:
    def test_STORAGE_URL(self):
        assert STORAGE.URL.dem.ace == URL.dem.ace
        assert STORAGE.URL.orbit.doris == URL.orbit.doris
        assert STORAGE.URL.auxcal.ers == URL.auxcal.ers
    
    def test_STORAGE_LOOKUP(self):
        assert LOOKUP.snap.suffix['Apply-Orbit-File'] == STORAGE.LOOKUP.snap.suffix['Apply-Orbit-File']
        assert LOOKUP.snap.suffix['Terrain-Correction'] == STORAGE.LOOKUP.snap.suffix['Terrain-Correction'] == 'TC'
        assert LOOKUP.attributes['sensor'] == STORAGE.LOOKUP.attributes['sensor']
        assert LOOKUP.attributes['vh'] == STORAGE.LOOKUP.attributes['vh']


class TestConfigHandler:
    
    def test_make_dir_and_config(self, tmpdir):
        config = os.path.join(str(tmpdir), 'unit_test_config.ini')
        conf = ConfigHandler(config_fname=config)
        
        path_pyrosar = os.path.exists(conf._ConfigHandler__GLOBAL['path'])
        path_config = os.path.isfile(conf._ConfigHandler__GLOBAL['config'])
        
        assert path_pyrosar is True
        assert path_config is True
    
    def test_add_section(self, tmpdir):
        config = os.path.join(str(tmpdir), 'unit_test_config.ini')
        conf = ConfigHandler(config_fname=config)
        conf.add_section('SNAP')
        
        assert conf.sections[0] == 'SNAP'
    
    def test_options(self, tmpdir):
        config = os.path.join(str(tmpdir), 'unit_test_config.ini')
        conf = ConfigHandler(config_fname=config)
        conf.add_section('SNAP')
        conf.set('SNAP', 'etc', 'temp/dir')
        
        with pytest.raises(AttributeError):
            conf.set('SNAPp', 'etc', 'temp/dir')
        
        assert conf['SNAP']['etc'] == 'temp/dir'
        assert conf['SNAP'] == {'etc': 'temp/dir'}
        assert conf.sections == ['SNAP']
    
    def test_overwrite(self, tmpdir):
        config = os.path.join(str(tmpdir), 'unit_test_config.ini')
        conf = ConfigHandler(config_fname=config)
        conf.add_section('SNAP')
        conf.set('SNAP', 'etc', 'temp/dir')
        
        with pytest.raises(RuntimeError):
            conf.set('SNAP', 'etc', 'temp/dir2')
        
        conf.set('SNAP', 'etc', 'temp/dir2', True)
        assert conf['SNAP']['etc'] == 'temp/dir2'
    
    def test_remove(self, tmpdir):
        config = os.path.join(str(tmpdir), 'unit_test_config.ini')
        conf = ConfigHandler(config_fname=config)
        conf.add_section('SNAP')
        conf.set('SNAP', 'etc', 'temp/dir')
        
        with pytest.raises(AttributeError):
            conf.remove_option('SNAP', 'kex')
        
        with pytest.raises(AttributeError):
            conf.remove_option('SNApP', 'etc')
        
        conf.remove_option('SNAP', 'etc')
        assert list(conf['SNAP'].keys()) == []
    
    def test_delete_unit_data(self, tmpdir):
        config = os.path.join(str(tmpdir), 'unit_test_config.ini')
        conf = ConfigHandler(config_fname=config)
        
        os.remove(conf._ConfigHandler__GLOBAL['config'])
