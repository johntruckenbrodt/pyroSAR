from pyroSAR.ancillary import seconds, groupbyTime, groupby


def test_seconds():
    assert seconds('test_20151212T234411') == 3658952651.0


def test_groupby():
    """
    Test correct grouping of filenames by their attributes
    Methodology is to provide a list of partially overlapping filenames
    and ensure the resultant list of lists contains the correct entry numbers
    """
    filenames = ['S1A__IW___A_20150309T173017_VV_grd_mli_geo_norm_db.tif',
                 'S1A__IW___A_20150309T173017_HH_grd_mli_geo_norm_db.tif',
                 'S2A__IW___A_20180309T173017_HH_grd_mli_geo_norm_db.tif']
    sensor_groups = groupby(filenames, 'sensor')
    print(sensor_groups)
    assert len(sensor_groups) == 2
    assert isinstance(sensor_groups[0], list)
    assert len(sensor_groups[0]) == 2

    filenames += ['S2A__IW___A_20180309T173017_VV_grd_mli_geo_norm_db.tif']

    polarization_groups = groupby(filenames, 'polarization')
    print(polarization_groups)
    assert len(polarization_groups) == 2
    assert isinstance(polarization_groups[0], list)
    assert isinstance(polarization_groups[1], list)
    assert len(polarization_groups[0]) == 2
    assert len(polarization_groups[1]) == 2

    filenames += ['S2A__IW___A_20180309T173017_HV_grd_mli_geo_norm_db.tif']

    polarization_groups = groupby(filenames, 'polarization')
    print(polarization_groups)
    assert len(polarization_groups) == 3
    assert isinstance(polarization_groups[0], list)
    assert isinstance(polarization_groups[1], list)
    assert isinstance(polarization_groups[2], list)
    assert len(polarization_groups[0]) == 2
    assert len(polarization_groups[1]) == 1
    assert len(polarization_groups[2]) == 2


def test_groupbyTime():
    filenames = ['S1__IW___A_20151212T120000',
                 'S1__IW___A_20151212T120100',
                 'S1__IW___A_20151212T120300']
    groups = groupbyTime(filenames, seconds, 60)
    print(groups)
    assert len(groups) == 2
    assert isinstance(groups[0], list)
    assert len(groups[0]) == 2

    filenames = ['S1__IW___A_20151212T120000',
                 'S1__IW___A_20151212T120100',
                 'S1__IW___A_20151212T120200']
    groups = groupbyTime(filenames, seconds, 60)
    print(groups)
    assert len(groups[0]) == 3
