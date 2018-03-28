import pytest

from pyroSAR._dev_config import Storage, LOOKUP, URL, STORAGE


class TestStorage:
    def test_insert(self):
        storage = Storage(a=1, b=2)
        assert storage.a == 1
        assert storage.b == 2

class TestLookup:
    def test_suffix(self):
        assert LOOKUP.snap.suffix[0]['Apply-Orbit-File'] == 'Orb'
        assert LOOKUP.snap.suffix[0]['Terrain-Correction'] == 'TC'

    def test_attributes(self):
        assert LOOKUP.attributes['sensor'] == 'TEXT'
        assert LOOKUP.attributes['vh'] == 'INTEGER'

class TestSTORAGE:
    def test_STORAGE_URL(self):
        assert STORAGE.URL.dem.ace == URL.dem.ace
        assert STORAGE.URL.orbit.doris == URL.orbit.doris
        assert STORAGE.URL.auxcal.ers == URL.auxcal.ers

    def test_STORAGE_LOOKUP(self):
        assert LOOKUP.snap.suffix[0]['Apply-Orbit-File'] == STORAGE.LOOKUP.snap.suffix[0]['Apply-Orbit-File']
        assert LOOKUP.snap.suffix[0]['Terrain-Correction'] == STORAGE.LOOKUP.snap.suffix[0]['Terrain-Correction'] == 'TC'
        assert LOOKUP.attributes['sensor'] == STORAGE.LOOKUP.attributes['sensor']
        assert LOOKUP.attributes['vh'] == STORAGE.LOOKUP.attributes['vh']