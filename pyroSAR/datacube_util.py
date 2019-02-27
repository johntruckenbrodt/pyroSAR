"""
This (still experimental) module is intended to easily prepare SAR scenes processed
by pyroSAR for ingestion into an Open Data Cube.

.. code-block:: python

    from pyroSAR.datacube_util import Product, Dataset
    from pyroSAR.ancillary import find_datasets
    
    # find pyroSAR files by metadata attributes
    archive_s1 = '/.../sentinel1/GRD/processed'
    scenes_s1 = find_datasets(archive_s1, sensor=('S1A', 'S1B'), acquisition_mode='IW')
    
    # define the polarization units describing the data sets
    units = {'VV': 'backscatter VV', 'VH': 'backscatter VH'}
    
    # create a new product
    with Product(name='S1_GRD_index',
                 product_type='gamma0',
                 description='Gamma Naught RTC backscatter') as prod:
        
        for dataset in scenes_s1:
            with Dataset(dataset, units=units) as ds:
                
                # add the dataset to the product
                prod.add(ds)
                
                # parse datacube indexing YMLs from product and data set metadata
                prod.export_indexing_yml(ds, 'yml_index_outdir')
        
        # write the product YML
        prod.write('yml_product')
        
        # print the product metadata which is written to the product YML
        print(prod)
"""

import os
import re
import yaml
import uuid
from time import strftime, strptime
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
            # map pyroSAR sensor identifiers to platform and instrument codes
            sensor_lookup = {'ASAR': ('ENVISAT', 'ASAR'),
                             'ERS1': ('ERS-1', 'SAR'),
                             'ERS2': ('ERS-2', 'SAR'),
                             'PSR1': ('ALOS-1', 'PALSAR'),
                             'PSR2': ('ALOS-2', 'PALSAR-2'),
                             'S1A': ('SENTINEL-1', 'C-SAR'),
                             'S1B': ('SENTINEL-1', 'C-SAR'),
                             'TSX1': ('TERRASAR-X_1', 'SAR'),
                             'TDX1': ('TANDEM-X_1', 'SAR')}
            
            # extract basic metadata attributes from the filename and register them to the object
            meta = parse_datasetname(filename)
            
            if meta is None:
                raise ValueError('could not identify dataset: {}'.format(filename))
            
            for key, val in meta.items():
                setattr(self, key, val)
            
            # define acquisition start and end time; Currently both are set to the acquisition start time,
            # which is contained in the filename
            # Time will only be correct if the full scene was processed, start and end time of s subset will
            # differ. Thus, accurately setting both is not seen as too relevant.
            self.from_dt = strftime('%Y-%m-%dT%H:%M:%S', strptime(self.start, '%Y%m%dT%H%M%S'))
            self.to_dt = strftime('%Y-%m-%dT%H:%M:%S', strptime(self.start, '%Y%m%dT%H%M%S'))
            
            # match the sensor ID from the filename to a platform and instrument
            if self.sensor not in sensor_lookup.keys():
                raise ValueError('unknown sensor: {}'.format(self.sensor))
            
            self.platform, self.instrument = sensor_lookup[self.sensor]
            
            # extract general geo metadata from the GTiff information
            with Raster(filename) as ras:
                self.dtype = Dtype(ras.dtype).numpystr
                self.nodata = ras.nodata
                self.format = ras.format
                self.xres, self.yres = ras.res
                self.crs = 'EPSG:{}'.format(ras.epsg)
                self.is_projected = ras.projcs is not None
                self.extent = self.__extent_convert(ras.geo, 'x', 'y')
                # reproject the raster bounding box to EPSG 4326 and store its extent
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
            
            # create the measurement entry from collected metadata;
            # this is intended for easy access by class Product
            self.measurements = {self.polarization: {'dtype': self.dtype,
                                                     'name': self.polarization,
                                                     'nodata': self.nodata,
                                                     'filename': filename,
                                                     'units': units}}
        else:
            raise TypeError('filename must be of type str, list or Dataset')
    
    def __add__(self, dataset):
        """
        override the + operator. This is intended to easily combine two Dataset objects, which were
        created from different files belonging to the same measurement, e.g. two GeoTiffs with one polarization
        each.
        
        Parameters
        ----------
        dataset: Dataset
            the dataset to add to the current one

        Returns
        -------
        Dataset
            the combination of the two
        """
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
        """
        similar to Dataset.__add__ but for function sum, e.g. sum([Dataset1, Dataset2])
        
        Parameters
        ----------
        dataset: Dataset
            the dataset to add to the current one

        Returns
        -------
        Dataset
            the combination of the two
        """
        if dataset == 0:
            return self
        else:
            return self.__add__(dataset)
    
    @staticmethod
    def __extent_convert(extent, xkey, ykey):
        """
        convert the extent of a :class:`~spatialist.raster.Raster` object to a
        datacube-compliant dictionary.
        
        Parameters
        ----------
        extent: dict
            the extent as returned by a :class:`~spatialist.raster.Raster` object
        xkey: {'longitude', 'x'}
            the key of the x dimension
        ykey: {'latitude', 'y'}
            the key of the y dimension

        Returns
        -------
        dict
            a dictionary with keys `ll`, `lr`, `ul` and ``ur
        """
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
        """
        get a certain measurement attribute from all measurements
        
        Parameters
        ----------
        attr: str
            the attribute to get

        Returns
        -------
        dict
            a dictionary with the measurement names as keys and the respective attribute as value
        """
        return dict([(key, self.measurements[key][attr]) for key in self.measurements.keys()])
    
    @property
    def filenames(self):
        """
        
        Returns
        -------
        dict
            all file names registered in the dataset
        """
        return self.__get_measurement_attr('filename')
    
    @property
    def identifier(self):
        """
        
        Returns
        -------
        str
            a unique dataset identifier
        """
        return '{}_{}'.format(self.outname_base, '_'.join(self.proc_steps))
    
    @property
    def units(self):
        """
        
        Returns
        -------
        dict
            all measurement unit names registered in the dataset
        """
        return self.__get_measurement_attr('units')
    
    @units.setter
    def units(self, value):
        """
        (re)set the units of all measurements
        
        Parameters
        ----------
        value: str or dict
            the unit(s) to be set; if multiple measurements are present,
            a dictionary with measurement names as keys needs to be defined

        Returns
        -------

        """
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
        
        missing_message = "when initializing {}, parameters " \
                          "'name', 'product_type' and 'description' must be defined"
        
        if isinstance(definition, str):
            if os.path.isfile(definition):
                with open(definition, 'r') as yml:
                    try:
                        self.meta = yaml.load(yml)
                    except yaml.YAMLError:
                        raise RuntimeError('the provided file does not seem to be a YAML file')
            else:
                raise RuntimeError('definition file does not exist')
        
        elif isinstance(definition, list):
            if None in [name, product_type, description]:
                raise ValueError(missing_message.format(' a product from list'))
            self.__initialize(name, product_type, description)
            for dataset in definition:
                with Dataset(dataset) as DS:
                    self.add(DS)
        
        elif definition is None:
            if None in [name, product_type, description]:
                raise ValueError(missing_message.format('a blank product'))
            self.__initialize(name, product_type, description)
        else:
            raise TypeError('type of parameter definition must be either str, list or None')
    
    def __enter__(self):
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
    
    def __str__(self):
        return yaml.dump(self.meta, default_flow_style=False)
    
    def __getattr__(self, item):
        if item in self.__fixture_storage:
            return self.meta['storage'][item]
        elif item in self.__fixture_metadata:
            subkey = 'code' if item == 'platform' else 'name'
            return self.meta['metadata'][item][subkey]
        elif item == 'product_type':
            return self.meta['metadata']['product_type']
        else:
            return object.__getattribute__(self, item)
    
    def __setattr__(self, key, value):
        if key in self.__fixture_storage:
            self.meta['storage'][key] = value
        elif key in self.__fixture_metadata:
            subkey = 'code' if key == 'platform' else 'name'
            self.meta['metadata'][key][subkey] = value
        elif key == 'product_type':
            self.meta['metadata']['product_type'] = value
        else:
            super(Product, self).__setattr__(key, value)
    
    def close(self):
        return
    
    def __add_measurement(self, name, dtype, nodata, units):
        """
        create a new measurement entry
        
        Parameters
        ----------
        name: str
            the measurement name
        dtype: str
            the data type, e.g. float32
        nodata: int or float
            the nodata value of the data
        units: str
            the measurement units

        Returns
        -------

        """
        if name in self.measurements.keys():
            raise IndexError('measurement {} already exists'.format(name))
        self.meta['measurements'].append({'name': name,
                                          'dtype': dtype,
                                          'units': units,
                                          'nodata': nodata})
    
    def __initialize(self, name, product_type, description):
        """
        create a new blank product
        
        Parameters
        ----------
        name: str
            the name of the product
        product_type: str
            the product type, e.g. `gamma0`
        description: str
            a description of the product content/purpose

        Returns
        -------

        """
        self.meta = {'description': description,
                     'measurements': [],
                     'metadata': {'platform': {'code': None},
                                  'instrument': {'name': None},
                                  'format': {'name': None},
                                  'product_type': product_type},
                     'metadata_type': 'eo',
                     'name': name,
                     'storage': {'crs': None,
                                 'resolution': None}}
    
    @staticmethod
    def __check_dict_keys(keys, reference):
        return len(union(keys, reference)) == len(keys)
    
    @property
    def __fixture_fields(self):
        """
        
        Returns
        -------
        list
            the names of the top-level metadata fields, which must be defined
        """
        return ['description', 'measurements', 'metadata', 'metadata_type', 'name', 'storage']
    
    @property
    def __fixture_measurement(self):
        """
        
        Returns
        -------
        list
            the names of the metadata fields, which must be defined for all measurements
        """
        return ['dtype', 'nodata', 'units']
    
    @property
    def __fixture_metadata(self):
        """
        
        Returns
        -------
        list
            the names of the metadata fields, which must be defined in the general metadata section
        """
        return ['format', 'instrument', 'platform']
    
    @property
    def __fixture_storage(self):
        """
        
        Returns
        -------
        list
            the names of the metadata fields, which must be defined for the storage section
        """
        return ['crs', 'resolution']
    
    def __validate(self):
        """
        assert whether the Product is valid
        
        Returns
        -------
        
        Raises
        ------
        RuntimeError
        """
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
        """
        Add a dataset to the abstracted product description. This first performs a check
        whether the dataset is compatible with the product and its already existing measurements.
        If a measurement in the dataset does not yet exist in the product description it is added.
        
        Parameters
        ----------
        dataset: Dataset
            the dataset whose description is to be added

        Returns
        -------

        """
        if not isinstance(dataset, Dataset):
            raise TypeError('input must be of type pyroSAR.datacube.Dataset')
        self.check_integrity(dataset, allow_new_measurements=True)
        
        # set the general product definition attributes if they are None
        for attr in self.__fixture_metadata + self.__fixture_storage:
            if getattr(self, attr) is None:
                setattr(self, attr, getattr(dataset, attr))
        
        # if it is not yet present, add the dataset measurement definition to that of the product
        for measurement, content in dataset.measurements.items():
            if measurement not in self.measurements.keys():
                self.__add_measurement(dtype=content['dtype'],
                                       name=content['name'],
                                       nodata=content['nodata'],
                                       units=content['units'])
    
    def check_integrity(self, dataset, allow_new_measurements=False):
        """
        check if a dataset is compatible with the product definition.
        
        Parameters
        ----------
        dataset: Dataset
            the dataset to be checked
        allow_new_measurements: bool
            allow new measurements to be added to the product definition?
            If not and the dataset contains measurements,
            which are not defined in the product, an error is raised.

        Returns
        -------
        
        Raises
        ------
        RuntimeError
        """
        # check general metadata and storage fields
        for attr in self.__fixture_metadata + self.__fixture_storage:
            val_ds = getattr(dataset, attr)
            val_prod = getattr(self, attr)
            if val_prod is not None and val_ds != val_prod:
                raise RuntimeError("mismatch of attribute '{0}': {1}, {2}".format(attr, val_ds, val_prod))
        
        # check measurement fields
        for measurement, content in dataset.measurements.items():
            if measurement not in self.measurements.keys():
                if not allow_new_measurements:
                    raise RuntimeError("measurement '{}' is not present in the product definition "
                                       "and allow_new_measurements is set to False".format(measurement))
            else:
                match = self.measurements[measurement]
                for attr in self.__fixture_measurement:
                    if match[attr] != content[attr]:
                        raise RuntimeError("mismatch of measurement '{0}', "
                                           "attribute '{1}': {2}, {3}".
                                           format(measurement, attr, match[attr], content[attr]))
    
    def export_indexing_yml(self, dataset, outdir):
        
        self.__validate()
        
        outname = os.path.join(outdir, dataset.identifier + '_dcindex.yml')
        
        if os.path.isfile(outname):
            raise RuntimeError('indexing YML already exists: \n   {}'.format(outname))
        
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
            subkey = 'code' if attr == 'platform' else 'name'
            out[attr] = {subkey: getattr(dataset, attr)}
        
        out['grid_spatial']['projection']['geo_ref_points'] = dataset.extent
        out['grid_spatial']['projection']['spatial_reference'] = dataset.crs
        
        out['extent']['coord'] = dataset.extent_4326
        out['extent']['from_dt'] = dataset.from_dt
        out['extent']['to_dt'] = dataset.to_dt
        
        out['product_type'] = self.meta['metadata']['product_type']
        
        with open(outname, 'w') as yml:
            yaml.dump(out, yml, default_flow_style=False)
    
    def export_ingestion_yml(self, outname, product_name, ingest_location):
        
        if os.path.isfile(outname):
            raise RuntimeError('product definition YML already exists: \n   {}'.format(outname))
        
        self.__validate()
        
        if product_name == self.meta['name']:
            raise ValueError('source and target product names must be different')

        outdir = os.path.dirname(outname)
        if not os.path.isdir(outdir):
            os.makedirs(outdir)
        
        file_path_template = '{0}/{1}_{2}_{3}_{4}_' \
                             '{{tile_index[0]}}_' \
                             '{{tile_index[1]}}_' \
                             '{{start_time}}.nc'.format(product_name,
                                                        self.platform,
                                                        self.instrument,
                                                        self.product_type,
                                                        self.crs.replace('EPSG:', ''))
        
        global_attributes = {'instrument': self.instrument,
                             'platform': self.platform,
                             'institution': 'ESA',
                             'achknowledgment': 'Sentinel-1 data is provided by the European Space Agency '
                                                'on behalf of the European Commission via download.'}
        
        storage = self.meta['storage']
        storage['driver'] = 'NetCDF CF'
        storage['tile_size'] = {}
        for key in storage['resolution']:
            storage['tile_size'][key] = storage['resolution'][key] * 500
            storage['tile_size'][key] = storage['resolution'][key] * 500
        storage['chunking'] = {'time': 1}
        for key in storage['resolution']:
            storage['chunking'][key] = 500
            storage['chunking'][key] = 500
        storage['dimension_order'] = ['time', 'y', 'x']
        
        measurements = self.meta['measurements']
        for measurement in measurements:
            measurement['resampling_method'] = 'nearest'
            measurement['src_varname'] = measurement['name']
        
        out = {'source_type': self.meta['name'],
               'output_type': product_name,
               'description': self.meta['description'],
               'location': ingest_location,
               'file_path_template': file_path_template,
               'storage': self.meta['storage'],
               'measurements': self.meta['measurements'],
               'global_attributes': global_attributes}
        
        with open(outname, 'w') as yml:
            yaml.dump(out, yml, default_flow_style=False)
    
    @property
    def measurements(self):
        """
        
        Returns
        -------
        dict of dict
            a dictionary with measurement names as keys
        """
        return dict([(x['name'], x) for x in self.meta['measurements']])
    
    def write(self, ymlfile):
        """
        write the product definition to a YML file
        
        Parameters
        ----------
        ymlfile: str
            the file to write to
        
        Returns
        -------
        
        """
        if os.path.isfile(ymlfile):
            raise RuntimeError('ingestion YML already exists: \n   {}'.format(ymlfile))
        
        self.__validate()
        with open(ymlfile, 'w') as yml:
            yaml.dump(self.meta, yml, default_flow_style=False)
