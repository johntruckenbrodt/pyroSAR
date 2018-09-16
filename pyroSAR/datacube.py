import os
import re
import yaml
import uuid
from spatialist.raster import Raster, Dtype
from spatialist.ancillary import union
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


class Dataset(object):
    def __init__(self, filename):
        
        self.filename = os.path.realpath(filename)
        
        sensor_lookup = {'ASAR': ('ENVISAT', 'ASAR'),
                         'ERS1': ('ERS-1', 'SAR'),
                         'ERS2': ('ERS-2', 'SAR'),
                         'PSR1': ('ALOS-1', 'PALSAR'),
                         'PSR2': ('ALOS-2', 'PALSAR-2'),
                         'S1A': ('SENTINEL-1', 'C-SAR'),
                         'S1B': ('SENTINEL-1', 'C-SAR'),
                         'TSX1': ('TERRASAR-X_1', 'SAR'),
                         'TDX1': ('TANDEM-X_1', 'SAR')}
        
        meta = parse_datasetname(filename)
        
        if meta is None:
            raise ValueError('could not identify dataset: {}'.format(filename))
        
        for key, val in meta.items():
            setattr(self, key, val)
        
        self.measurement = self.polarization
        
        if self.sensor not in sensor_lookup.keys():
            raise ValueError('unknown sensor: {}'.format(self.sensor))
        
        self.platform, self.instrument = sensor_lookup[self.sensor]
        
        with Raster(filename) as ras:
            self.dtype = Dtype(ras.dtype).numpystr
            self.nodata = ras.nodata
            self.format = ras.format
            self.xres, self.yres = ras.res
            self.crs = 'EPSG:{}'.format(ras.epsg)
            self.is_projected = ras.projcs is not None
            self.extent = self.__extent_convert(ras.geo, 'x', 'y')
            with ras.bbox() as bbox:
                bbox.reproject(4236)
                self.extent_4326 = self.__extent_convert(bbox.extent, 'lon', 'lat')
        
        # create dictionary for resolution metadata depending on CRS characteristics
        resolution_keys = ('x', 'y') if self.is_projected else ('longitude', 'latitude')
        self.resolution = dict(zip(resolution_keys, (self.xres, self.yres)))
        
        # check whether the data type is supported
        pattern = '(?:(?:u|)int(?:8|16|32|64)|float(?:32|64))'
        if not re.search(pattern, self.dtype):
            raise ValueError('unsupported data type {}'.format(self.dtype))
    
    @staticmethod
    def __extent_convert(extent, xkey, ykey):
        return {'ll': {xkey: extent['xmin'],
                       ykey: extent['ymin']},
                'lr': {xkey: extent['xmax'],
                       ykey: extent['ymin']},
                'ul': {xkey: extent['xmin'],
                       ykey: extent['ymax']},
                'ur': {xkey: extent['xmax'],
                       ykey: extent['ymax']}}
    
    def __enter__(self):
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
    
    def close(self):
        return


class Product(object):
    def __init__(self, definition=None, name=None, product_type=None,
                 description=None, measurement_units=None):
        self.measurement_units = measurement_units if measurement_units is not None else {}
        
        missing_message = "when initializing from {}, parameters " \
                          "'name', 'product_type' and 'description' must be defined"
        
        if isinstance(definition, str) and os.path.isfile(definition):
            with open(definition, 'r') as yml:
                try:
                    self.meta = yaml.load(yml)
                except yaml.YAMLError:
                    raise RuntimeError('the provided file does not seem to be a YAML file')
        
        elif isinstance(definition, list):
            if None in [name, product_type, description]:
                raise ValueError(missing_message.format('list'))
            self.__initialize(name, product_type, description)
            for dataset in definition:
                with Dataset(dataset) as DS:
                    self.add(DS)
        
        elif definition is None:
            self.__initialize(name, product_type, description)
        else:
            raise TypeError('tpye of parameter definition must be either str, list or None')
    
    def __enter__(self):
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
    
    def __str__(self):
        return yaml.dump(self.meta, default_flow_style=False)
    
    def close(self):
        return
    
    def __add_measurement(self, name, dtype, nodata):
        if name not in self.measurement_units.keys():
            self.measurement_units[name] = 'DN'
        if name in self.measurements.keys():
            raise IndexError('measurement {} already exists'.format(name))
        self.meta['measurements'].append({'name': name,
                                          'dtype': dtype,
                                          'units': self.measurement_units[name],
                                          'nodata': nodata})
    
    def __initialize(self, name, product_type, description):
        self.meta = {'description': description,
                     'measurements': [],
                     'metadata': {'product_type': product_type},
                     'metadata_type': 'eo',
                     'name': name,
                     'storage': {}}
    
    @staticmethod
    def __check_dict_keys(keys, reference):
        return len(union(keys, reference)) == len(keys)
    
    @property
    def __fixture_fields(self):
        return ['description', 'measurements', 'metadata', 'metadata_type', 'name', 'storage']
    
    @property
    def __fixture_measurement(self):
        return ['dtype', 'nodata']
    
    @property
    def __fixture_metadata(self):
        return ['format', 'instrument', 'platform']
    
    @property
    def __fixture_storage(self):
        return ['crs', 'resolution']
    
    def __validate(self):
        try:
            assert isinstance(self.meta, dict)
            assert self.__check_dict_keys(self.__fixture_fields, self.meta.keys())
            assert 'product_type' in self.meta['metadata'].keys()
            for measurement in self.meta['measurements']:
                assert self.__check_dict_keys(self.__fixture_measurement, measurement.keys())
        except AssertionError as e:
            print(e)
            raise RuntimeError('product invalid')
    
    def add(self, dataset):
        if not isinstance(dataset, Dataset):
            raise TypeError('input must be of type pyroSAR.datacube.Dataset')
        self.check_integrity(dataset, allow_new_measurements=True)
        for attr in self.__fixture_metadata:
            self.meta['metadata'][attr] = {'name': getattr(dataset, attr)}
        
        for attr in self.__fixture_storage:
            self.meta['storage'][attr] = getattr(dataset, attr)
        
        if dataset.measurement not in self.measurements.keys():
            self.__add_measurement(dataset.measurement, dataset.dtype, dataset.nodata)
    
    def check_integrity(self, dataset, allow_new_measurements=False):
        
        # define parameters, which are to be consistent across an individual measurement
        # fixtures_measurement = ['acquisition_mode', 'proc_steps', 'dtype', 'nodata']
        
        # check metadata field
        for attr in self.__fixture_metadata:
            if attr in self.meta['metadata'].keys():
                if getattr(dataset, attr) != self.meta['metadata'][attr]['name']:
                    raise ValueError('mismatch of attribute {}'.format(attr))
        
        # check storage field
        for attr in self.__fixture_storage:
            if attr in self.meta['storage'].keys():
                if getattr(dataset, attr) != self.meta['storage'][attr]:
                    raise ValueError('mismatch of attribute {}'.format(attr))
        
        if dataset.measurement not in self.measurements.keys():
            if not allow_new_measurements:
                raise ValueError('measurement mismatch')
        else:
            match = self.measurements[dataset.measurement]
            for attr in self.__fixture_measurement:
                if match[attr] != getattr(dataset, attr):
                    raise ValueError('mismatch of attribute {}'.format(attr))
    
    def export_ingestion_yml(self, dataset):
        self.check_integrity(dataset)
        out = {'id': str(uuid.uuid4()),
               'image': {'bands': {}},
               'grid_spatial': {'projection': {}},
               'extent': {'coord': {}},
               'lineage': {'source_datasets': {}}}
        
        out['image']['bands'][dataset.measurement] = {'path': dataset.filename}
        
        for attr in self.__fixture_metadata:
            out[attr] = {'name': getattr(dataset, attr)}
        
        out['grid_spatial']['projection']['geo_ref_points'] = dataset.extent
        out['grid_spatial']['projection']['spatial_reference'] = dataset.crs
        
        out['extent']['coord'] = dataset.extent_4326
        
        out['product_type'] = self.meta['metadata']['product_type']
        
        print(yaml.dump(out, default_flow_style=False))
    
    @property
    def measurements(self):
        return dict([(x['name'], x) for x in self.meta['measurements']])
    
    def write(self, ymlfile):
        with open(ymlfile, 'w') as yml:
            yaml.dump(self.meta, yml, default_flow_style=False)
