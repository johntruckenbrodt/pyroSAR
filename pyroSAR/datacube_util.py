import os
import re
import yaml
import uuid
from spatialist.raster import Raster, Dtype
from spatialist.ancillary import union
from .ancillary import parse_datasetname


class Dataset(object):
    def __init__(self, filename, units='DN'):
        
        if isinstance(filename, list):
            combined = sum([Dataset(x, units) for x in filename])
            self.__init__(combined)
        
        elif isinstance(filename, Dataset):
            for attr, value in vars(filename).items():
                setattr(self, attr, value)
        
        elif isinstance(filename, str):
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
                    bbox.reproject(4326)
                    self.extent_4326 = self.__extent_convert(bbox.extent, 'lon', 'lat')
            
            # create dictionary for resolution metadata depending on CRS characteristics
            resolution_keys = ('x', 'y') if self.is_projected else ('longitude', 'latitude')
            self.resolution = dict(zip(resolution_keys, (self.xres, self.yres)))
            
            # check whether the data type is supported
            pattern = '(?:(?:u|)int(?:8|16|32|64)|float(?:32|64))'
            if not re.search(pattern, self.dtype):
                raise ValueError('unsupported data type {}'.format(self.dtype))
            
            # determine the dataset units
            if isinstance(units, str):
                units = units
            elif isinstance(units, dict):
                try:
                    units = units[self.polarization]
                except KeyError:
                    raise KeyError("parameter 'units' does not contain key '{}'".format(self.polarization))
            else:
                raise TypeError("parameter 'units' must be of type str or dict")
            
            # create the measurement entry from collected metadata; this is intended for easy access by class Product
            self.measurements = {self.polarization: {'dtype': self.dtype,
                                                     'name': self.polarization,
                                                     'nodata': self.nodata,
                                                     'filename': filename,
                                                     'units': units}}
        else:
            raise TypeError('filename must be of type str, list or Dataset')
    
    def __add__(self, dataset):
        for attr in ['extent', 'crs', 'sensor', 'acquisition_mode', 'proc_steps', 'outname_base']:
            if getattr(self, attr) != getattr(dataset, attr):
                raise ValueError('value mismatch: {}'.format(attr))
        # self.filename.append(dataset.filename)
        for key in dataset.measurements.keys():
            if key in self.measurements.keys():
                raise RuntimeError('only different measurements can be combined to one dataset')
        self.measurements.update(dataset.measurements)
        return self
    
    def __radd__(self, dataset):
        if dataset == 0:
            return self
        else:
            return self.__add__(dataset)
    
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
    
    def __get_measurement_attr(self, attr):
        return dict([(key, self.measurements[key][attr]) for key in self.measurements.keys()])
    
    @property
    def filenames(self):
        return self.__get_measurement_attr('filename')
    
    @property
    def identifier(self):
        return '{}_{}'.format(self.outname_base, '_'.join(self.proc_steps))
    
    @property
    def units(self):
        return self.__get_measurement_attr('units')
    
    @units.setter
    def units(self, value):
        keys = list(self.measurements.keys())
        if isinstance(value, str):
            if len(keys) == 1:
                self.measurements[keys[0]]['units'] = value
            else:
                raise TypeError('the dataset contains multiple measurements; '
                                'in this case a dictionary is needed for setting the measurement units')
        elif isinstance(value, dict):
            for name, unit in value.items():
                if name in keys:
                    self.measurements[name]['units'] = unit
                else:
                    raise KeyError("the dataset does not contain a measurement '{}'".format(name))
    
    def close(self):
        return


class Product(object):
    def __init__(self, definition=None, name=None, product_type=None,
                 description=None):
        
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
    
    def __add_measurement(self, name, dtype, nodata, units):
        if name in self.measurements.keys():
            raise IndexError('measurement {} already exists'.format(name))
        self.meta['measurements'].append({'name': name,
                                          'dtype': dtype,
                                          'units': units,
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
        return ['dtype', 'nodata', 'units']
    
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
            key = 'code' if attr == 'platform' else 'name'
            self.meta['metadata'][attr] = {key: getattr(dataset, attr)}
        
        for attr in self.__fixture_storage:
            self.meta['storage'][attr] = getattr(dataset, attr)
        
        for measurement, content in dataset.measurements.items():
            if measurement not in self.measurements.keys():
                self.__add_measurement(content['name'], content['dtype'], content['nodata'], content['units'])
    
    def check_integrity(self, dataset, allow_new_measurements=False):
        
        # define parameters, which are to be consistent across an individual measurement
        # fixtures_measurement = ['acquisition_mode', 'proc_steps', 'dtype', 'nodata']
        
        # check metadata field
        for attr in self.__fixture_metadata:
            if attr in self.meta['metadata'].keys():
                subkey = 'code' if attr == 'platform' else 'name'
                if getattr(dataset, attr) != self.meta['metadata'][attr][subkey]:
                    raise ValueError('mismatch of attribute {}'.format(attr))
        
        # check storage field
        for attr in self.__fixture_storage:
            if attr in self.meta['storage'].keys():
                if getattr(dataset, attr) != self.meta['storage'][attr]:
                    raise ValueError('mismatch of attribute {}'.format(attr))
        
        # check measurement fields
        for measurement, content in dataset.measurements.items():
            if measurement not in self.measurements.keys():
                if not allow_new_measurements:
                    raise ValueError('measurement mismatch')
            else:
                match = self.measurements[measurement]
                for attr in self.__fixture_measurement:
                    if match[attr] != content[attr]:
                        raise ValueError("mismatch of attribute '{0}': "
                                         "{1}, {2}".format(attr, match[attr], content[attr]))
    
    def export_ingestion_yml(self, dataset, outdir):
        if not os.path.isdir(outdir):
            os.makedirs(outdir)
        
        self.check_integrity(dataset)
        out = {'id': str(uuid.uuid4()),
               'image': {'bands': {}},
               'grid_spatial': {'projection': {}},
               'extent': {'coord': {}},
               'lineage': {'source_datasets': {}}}
        
        for measurement, content in dataset.measurements.items():
            out['image']['bands'][measurement] = {'path': content['filename']}
        
        for attr in self.__fixture_metadata:
            out[attr] = {'name': getattr(dataset, attr)}
        
        out['grid_spatial']['projection']['geo_ref_points'] = dataset.extent
        out['grid_spatial']['projection']['spatial_reference'] = dataset.crs
        
        out['extent']['coord'] = dataset.extent_4326
        
        out['product_type'] = self.meta['metadata']['product_type']
        
        outname = os.path.join(outdir, dataset.identifier+'.yml')
        
        with open(outname, 'w') as yml:
            yaml.dump(out, yml, default_flow_style=False)
    
    @property
    def measurements(self):
        return dict([(x['name'], x) for x in self.meta['measurements']])
    
    def write(self, ymlfile):
        with open(ymlfile, 'w') as yml:
            yaml.dump(self.meta, yml, default_flow_style=False)
