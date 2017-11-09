# -*- coding: utf-8 -*-
"""
Created on Thu Nov 09 11:42:39 2017

@author: ibari
"""
from __future__ import division
import numpy as np

class Memorize(dict):
    """ Memorize results.
    """
    def __getattr__(self, name):
        try:
            return self[name]
        except KeyError:
            raise AttributeError(name)

    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__

    def __repr__(self):
        if self.keys():
            m = max(map(len, list(self.keys()))) + 1
            return '\n'.join([k.rjust(m) + ': ' + repr(v)
                              for k, v in sorted(self.items())])
        else:
            return self.__class__.__name__ + "()"

    def __dir__(self):
        return list(self.keys())
    
class SpatialResults(dict):
    """ Represents the statistical result.
    
    
    Returns
    ----------
    Note:
            The returns are attributes.
            
    slope :         float or ndarray
                    Slope of the regression line.

    intersect1 :    float or ndarray
                    If the input are complex, int1 represents the first
                    intersection point with the complex plain.

    intersect2 :    float or ndarray
                    If the input are complex, int2 represents the second
                    intersection  point with the complex plain.
    
    intercept :     float or ndarray
                    Intercept of the regression line.
                    
    rvalue :        float or ndarray
                    correlation coefficient
                    
    pvalue :        float or ndarray
                    two-sided p-value for a hypothesis test whose null 
                    hypothesis is that the slope is zero.
                    
    stderr :        float or ndarray
                    Standard error of the estimated gradient.

    Notes
    -----
    There may be additional attributes not listed above depending of the
    specific solver. Since this class is essentially a subclass of dict
    with attribute accessors, one can see which attributes are available
    using the `keys()` method.
    """
    def __getattr__(self, name):
        try:
            return self[name]
        except KeyError:
            raise AttributeError(name)

    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__

    def __repr__(self):
        if self.keys():
            m = max(map(len, list(self.keys()))) + 1
            return '\n'.join([k.rjust(m) + ': ' + repr(v)
                              for k, v in sorted(self.items())])
        else:
            return self.__class__.__name__ + "()"

    def __dir__(self):
        return list(self.keys())

# ------- GIS operations ------- #
def nan_to_num(arrays, no_data=-99999):
    """
    Convert nan values to a value
    """
    arrays[np.isnan(arrays)] = no_data

def subset (data, x, y):
    """
    Note
    ----------        
    Subsetting an Image with coordinates.
    
    Parameters
    ----------
    data:            array
                     Data to subset.
    
    area:            tuple with lists
                     Subset coordinates like ([450,477], [0,10]).
    Returns
    -------
    array_like
    
    """
# Subsetting an Image with a
    if isinstance(data, tuple):
        arrays = []
        for i in range(len(data)):
            arrays_subset = data[i][x[0]:x[1], y[0]:y[1]]
            nan_to_num(arrays_subset)
            arrays.append(arrays_subset)

        data = tuple(data)

    else:
        data = data[x[0]:x[1], y[0]:y[1]]
        nan_to_num(data)

    return data

# ------- Data Tests ------- #
def test_data(data, verbose=1):
    """Data Test
    Test if the imported data fulfilled the requirements. Requires a ImportResult class instance

    """
    col = data.shape[1]
    _col = col[1:] == col[:-1]

    row = data.shape[0]
    _row = row[1:] == row[:-1]

    size = data.size
    _size = size[1:] == size[:-1]
    
    if verbose is 2:
        if (_col and _row and _size):
            print "Status: All requirements are fulfilled"
    
    if (_col != True or _row != True):
        raise AssertionError("Status: Input dimensions must agree", 
                             'shapes: cols = {1}, rows = {2}'.format(col, row))        
    
    if (_col is not True or _row is not True or _size is not True):
        raise AssertionError("Status: Input dimensions and size must agree", 
                             'shapes: cols = {1}, rows = {2}, size: {3}'.format(col, row, size))
        
def test_import_type(file_name):
    string_check = test_string(file_name)
    tuple_check = test_tuple(file_name)
    
    if tuple_check:
        tuple_str_check = tuple([type(item) == str for item in file_name])
        tuple_str_check = np.all(np.all(tuple_str_check) == True)
    else:
        tuple_str_check = None
        
    return string_check, tuple_check, tuple_str_check

def test_string(data):
    return isinstance(data, str)

def test_tuple(data):
    return isinstance(data, tuple)