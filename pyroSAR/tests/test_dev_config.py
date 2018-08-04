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

    def test_make_dir_and_config(self):
        conf = ConfigHandler(config_fname='unit_test_config.ini')

        path_pyrosar = os.path.exists(conf._ConfigHandler__GLOBAL['path'])
        path_config = os.path.isfile(conf._ConfigHandler__GLOBAL['config'])

        assert path_pyrosar == True
        assert path_config == True

    def test_add_section(self):
        conf = ConfigHandler(config_fname='unit_test_config.ini')
        conf.add_section('SNAP')

        assert conf.sections[0] == 'SNAP'

        with pytest.raises(AttributeError) as excinfo:
            conf.add_section('SNAPP')
        assert "['SNAP']" in str(excinfo)

    def test_options(self):
        conf = ConfigHandler(config_fname='unit_test_config.ini')
        conf.add_section('SNAP')
        conf.set('SNAP', 'etc', 'temp/dir')

        with pytest.raises(AttributeError) as excinfo:
            conf.set('SNAPp', 'etc', 'temp/dir')

        assert conf.SNAP['etc'] == 'temp/dir'
        assert conf.get('SNAP')['etc'] == 'temp/dir'
        assert conf.get('SNAP') == {'etc': 'temp/dir'}
        assert conf.keys('SNAP') == ['etc']
        assert conf.sections == ['SNAP']

    def test_overwrite(self):
        conf = ConfigHandler(config_fname='unit_test_config.ini')
        conf.add_section('SNAP')
        conf.set('SNAP', 'etc', 'temp/dir')
        assert conf.SNAP['etc'] == 'temp/dir'

        conf.set('SNAP', 'etc', 'temp/dir2', True)
        assert conf.SNAP['etc'] == 'temp/dir2'

    def test_remove(self):
        conf = ConfigHandler(config_fname='unit_test_config.ini')
        conf.add_section('SNAP')
        conf.set('SNAP', 'etc', 'temp/dir')

        with pytest.raises(AttributeError) as excinfo:
            conf.remove_option('SNAP', 'kex')

        with pytest.raises(AttributeError) as excinfo:
            conf.remove_option('SNApP', 'etc')

        conf.remove_option('SNAP', 'etc')
        conf.keys('SNAP')
        assert conf.keys('SNAP') == []

    def test_delete_unit_data(self):
        conf = ConfigHandler(config_fname='unit_test_config.ini')

        os.remove(conf._ConfigHandler__GLOBAL['config'])
