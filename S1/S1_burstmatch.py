#!/usr/bin/env python

"""
******************************************************************************
*
* Program for matching overlapping bursts of two Sentinel-1 SLC sub-swaths using Gamma.
* Author:   John Truckenbrodt, john.truckenbrodt@uni-jena.de
*
******************************************************************************
* Copyright (c) 2016, John Truckenbrodt <john dot truckenbrodt at uni-jena dot de>
*
* Permission is hereby granted, free of charge, to any person obtaining a
* copy of this software and associated documentation files (the "Software"),
* to deal in the Software without restriction, including without limitation
* the rights to use, copy, modify, merge, publish, distribute, sublicense,
* and/or sell copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following conditions:
*
* The above copyright notice and this permission notice shall be included
* in all copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
* OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
* THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
* DEALINGS IN THE SOFTWARE.
 ****************************************************************************
"""

import os
import re
import argparse
from osgeo import ogr
import subprocess as sp
from collections import OrderedDict

ogr.UseExceptions()


def slc_burst_corners(parfile, tops_parfile):
    proc = sp.Popen(['SLC_burst_corners', parfile, tops_parfile], stdin=sp.PIPE, stdout=sp.PIPE, stderr=sp.PIPE)
    out, err = proc.communicate()
    if proc.returncode != 0:
        raise RuntimeError(err + '\nreading of burst corners failed')
    matches = [re.findall('[0-9.]+', x) for x in re.findall('Burst:[0-9 .]*', out)]
    bursts = {}
    for line in matches:
        lat = map(float, line[1::2])
        lon = map(float, line[2::2])
        bursts[int(line[0])] = zip(lat, lon)
    return OrderedDict(sorted(bursts.items()))


def coord2geom(coordinatelist):
    coordstrings = [' '.join(map(str, x)) for x in coordinatelist]
    coordstrings.append(coordstrings[0])
    return ogr.CreateGeometryFromWkt('POLYGON(({}))'.format(', '.join(coordstrings)))


def burstoverlap(slc1_parfile, slc1_tops_parfile, slc2_parfile, slc2_tops_parfile):
    b1 = slc_burst_corners(slc1_parfile, slc1_tops_parfile)
    b2 = slc_burst_corners(slc2_parfile, slc2_tops_parfile)
    master_select = []
    slave_select = []
    for mburst in b1.keys():
        for sburst in b2.keys():
            poly1 = coord2geom(b1[mburst])
            poly2 = coord2geom(b2[sburst])
            intersection = poly1.Intersection(poly2)
            area_ratio = intersection.GetArea() / poly1.GetArea()
            if area_ratio > 0.9:
                master_select.append(mburst)
                slave_select.append(sburst)
    return master_select, slave_select


def main(slc1_par, slc2_par, slc1_topspar, slc2_topspar, outfile):
    print ''
    print 'SLC1_PAR:      ', slc1_par
    print 'SLC1_TOPS_PAR: ', slc1_topspar
    print 'SLC2_PAR:      ', slc2_par
    print 'SLC2_TOPS_PAR: ', slc2_topspar
    slc1_select, slc2_select = burstoverlap(slc1_par, slc2_par, slc1_topspar, slc2_topspar)

    if len(slc1_select + slc2_select) == 0:
        outstring = 'no overlap'
    else:
        outstring = '\nIDs of overlapping bursts:\n' \
                   'SLC1:{0}\nSLC2:{1}\n'.format(', '.join(map(str, slc1_select)),
                                                 ', '.join(map(str, slc2_select)))

    if outfile:
        with open(outfile, 'w') as out:
            out.write(outstring)
    else:
        print outstring

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='This program returns the IDs of overlapping bursts of two S1 TOPS burst SLC swaths. '
                    'The returned IDs are to be used by the Gamma program SLC_copy_S1_TOPS to crop the two swaths to their common overlap.')
    parser.add_argument("-of", "--outfile", default=None, help="a textfile which the result is written to", type=str)
    parser.add_argument('SLC1_PAR', help='(input) SLC-1 parameter file for the TOPS burst SLC')
    parser.add_argument('SLC1_TOPS_PAR', help='(input) SLC-1 TOPS parameter file for the TOPS burst SLC')
    parser.add_argument('SLC2_PAR', help='(input) SLC-2 parameter file for the TOPS burst SLC')
    parser.add_argument('SLC2_TOPS_PAR', help='(input) SLC-2 TOPS parameter file for the TOPS burst SLC')
    args = parser.parse_args()
    main(args.SLC1_PAR, args.SLC1_TOPS_PAR, args.SLC2_PAR, args.SLC2_TOPS_PAR, args.outfile)
