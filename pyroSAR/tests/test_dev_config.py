import pytest

from pyroSAR._dev_config import Storage, LOOKUP, URL, STORAGE


class TestStorage:
    def test_insert(self):
        storage = Storage(a=1, b=2)
        assert storage.a == 1
        assert storage.b == 2

    def test_keys(self):
        storage = Storage(a=1, b=2)
        key = storage.keys()
        assert key[0] == 'a'
        assert key[1] == 'b'

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
