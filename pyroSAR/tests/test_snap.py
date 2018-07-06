from pyroSAR.snap import geocode


def test_geocode(tmpdir, testdata):
    scene = testdata['s1']
    geocode(scene, str(tmpdir), test=True)
