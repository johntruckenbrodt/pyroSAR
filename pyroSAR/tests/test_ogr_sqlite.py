from osgeo import ogr

driver = ogr.GetDriverByName('SQLite')

if driver is None:
    raise RuntimeError('OGR was built without SQLite driver')
