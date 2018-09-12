import re
import yaml
from spatialist.raster import Raster, Dtype
from .ancillary import parse_datasetname


def meta_collect(dataset):
    sensor_lookup = {'ASAR': ('ENVISAT', 'ASAR'),
                     'ERS1': ('ERS-1', 'SAR'),
                     'ERS2': ('ERS-2', 'SAR'),
                     'PSR1': ('ALOS-1', 'PALSAR'),
                     'PSR2': ('ALOS-2', 'PALSAR-2'),
                     'S1A': ('SENTINEL-1', 'C-SAR'),
                     'S1B': ('SENTINEL-1', 'C-SAR'),
                     'TSX1': ('TERRASAR-X_1', 'SAR'),
                     'TDX1': ('TANDEM-X_1', 'SAR')}
    
    meta = parse_datasetname(dataset)
    
    if meta is None:
        raise ValueError('could not identify dataset: {}'.format(dataset))
    
    meta['measurement'] = meta['polarization']
    
    if meta['sensor'] not in sensor_lookup.keys():
        raise ValueError('unknown sensor: {}'.format(meta['sensor']))
    
    meta['platform'], meta['instrument'] = sensor_lookup[meta['sensor']]
    with Raster(dataset) as ras:
        meta['dtype'] = Dtype(ras.dtype).numpystr
        meta['nodata'] = ras.nodata
        meta['format'] = ras.format
        meta['xres'], meta['yres'] = ras.res
        meta['crs'] = ras.epsg
        meta['is_projected'] = ras.projcs is not None
    
    # check whether the data type is supported
    pattern = '(?:(?:u|)int(?:8|16|32|64)|float(?:32|64))'
    if not re.search(pattern, meta['dtype']):
        raise ValueError('unsupported data type {}'.format(meta['dtype']))
    
    return meta


def __check_consistency(meta, reference, attributes):
    for attribute in attributes:
        if meta[attribute] != reference[attribute]:
            raise ValueError('inconsistency for attribute {0}: '
                             '{1} vs. {2}'.format(attribute,
                                                  reference[attribute],
                                                  meta[attribute]))


def parse_product_yaml(name, type, description, units, datasets, ymlfile):
    """
    create a product definition yaml file from a list of pyroSAR processed datasets.
    This yml file is intended for creating a new datacube product.
    
    Parameters
    ----------
    name: str
        the name of the product to be created
    type:
        the type of the product. E.g. `surface reflectance` or `gamma0`
    description: str
        a description of the product characteristics
    units: str
        the measurement units, e.g. `DN` or `percent`
    datasets: list of str
        a list of files which belong to the product.
        This list is easiest to be created with function :func:`~pyroSAR.ancillary.find_datasets`.
    ymlfile: str
        the name of the yml file
    
    Returns
    -------
    
    """
    # define parameters, which are to be consistent across all datasets and measurements
    fixtures_product = ['format', 'xres', 'yres', 'crs', 'platform', 'instrument']
    
    # define parameters, which are to be consistent across an individual measurement
    fixtures_measurement = ['acquisition_mode', 'proc_steps', 'dtype', 'nodata']
    
    # extract metadata for a reference dataset against which all other datasets are checked for consistency
    ref_product = meta_collect(datasets[0])
    
    # extract metadata from all datasets and check consistency with reference
    attributes = []
    for dataset in datasets:
        meta = meta_collect(dataset)
        __check_consistency(meta, ref_product, fixtures_product)
        attributes.append(meta)
    
    # group the attribute entries by measurements
    measurements = {}
    for item in attributes:
        measurement = item['measurement']
        if measurement not in measurements.keys():
            measurements[measurement] = []
        measurements[measurement].append(item)
    
    # check consistency of measurement-specific metadata across a selection of attributes
    for entries in measurements.values():
        # define reference dataset
        ref_measurement = entries[0]
        
        for dataset in entries:
            __check_consistency(dataset, ref_measurement, fixtures_measurement)
    
    # create dictionary for resolution metadata depending on CRS characteristics
    resolution_keys = ('x', 'y') if ref_product['is_projected'] else ('longitude', 'latitude')
    resolution_dict = dict(zip(resolution_keys, (ref_product['xres'], ref_product['yres'])))
    
    # create dictionary containing input for YAML file
    metadata = {'name': name,
                'description': description,
                'metadata_type': 'eo',
                'metadata': {'format': {'name': ref_product['format']},
                             'instrument': {'name': ref_product['instrument']},
                             'platform': {'code': ref_product['platform']},
                             'product_type': type},
                'measurements': [{'name': key,
                                  'units': units,
                                  'dtype': measurements[key][0]['dtype'],
                                  'nodata': measurements[key][0]['nodata']}
                                 for key in measurements.keys()],
                'storage': {'crs': 'EPSG:{}'.format(ref_product['crs']),
                            'resolution': resolution_dict}
                }
    
    # write YAML file
    with open(ymlfile, 'w') as yml:
        yaml.dump(metadata, yml, default_flow_style=False)
