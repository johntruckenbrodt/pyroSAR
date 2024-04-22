from pyroSAR.config import ConfigHandler
import os
import pytest


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
