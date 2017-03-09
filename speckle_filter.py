#!/usr/bin/env python2.7
# ______________________________________________________________________________________________________________
# |
# |	NAME:
# |	     quegan_multitemp_filter
# |
# |	PURPOSE:
# |	     To apply a filter to a multitemporal set of coregistered SAR images to produce output images that are:
# |              - individually unbiased
# |              - have minimum variance (so minimum speckle)
# |              - all have the same equivalent number of looks (ENL)
# |              - have essentially the same resolution as the original data
# |  DATE:
# |	     23 October 2015 (Joao Carreiras, NCEO/U.Sheffield, UK# j.carreiras@sheffield.ac.uk)
# |
# |  Modified:
# |		 28 October 2015 (Felix Cremer & Marcel Urban, FSU Jena)
# |
# |_____________________________________________________________________________________________________________
#

from time import asctime
from osgeo import gdal
from spatial import raster
import numpy as np
from gdalconst import *
from astropy.convolution import convolve, Box2DKernel
import argparse
import pyhht.emd as hht


# @profile
def quegan(infile, outfile, kernel, is_list=False,
           nodata=0.0, max_memory=10000):
    print "##########################################"
    print "Start:", asctime()

    if not is_list:
        data = raster.Raster(infile)
        band_count = data.bands
    else:
        band_count = len(infile)
        data = raster.Raster(infile[0])

    # Import Layerstack

    print "RASTER BAND COUNT: ", band_count
    print "Cols, Rows:"
    print data.cols, data.rows
    mean_image = [[] for i in range(band_count)]
    # print mean_image
    ratio = [[] for i in range(band_count)]
    output = [[] for i in range(band_count)]
    maxlines = max_memory // (data.cols * band_count * 5 * 4 / 1048576.)
    print int(round(data.rows / (maxlines - kernel + 1) + 0.5))
    box = Box2DKernel(kernel)
    driver = gdal.GetDriverByName("ENVI")
    out = driver.Create(outfile, data.cols, data.rows, data.bands, GDT_Float32)
    for i in range(int(round(data.rows / (maxlines - kernel + 1) + 0.5))):
        # Loop through Bands and rows
        print 'Streifennummer: '
        print i
        print asctime()

        for band in range(band_count):
            band += 1

            if is_list:
                data = raster.Raster(infile[band])
                if data.nodata:
                    nodata = data.nodata

            print 'BandNummer:'
            print band
            print asctime()
            srcband = np.split(data.matrix(band), [(maxlines - kernel + 1) * i, (maxlines - kernel + 1) * i + maxlines], axis=0)[1]
            srcband[srcband == nodata] = np.nan

            # print srcband.shape
            # print 'Test'
            # mean_image[band-1] = np.zeros(srcband)
            # print 'Test'
            # 1st STEP: apply a mean filter (in this case 7x7) to each separate image to estimate local mean
            # (mean_image(*, *, k) = mean_filter(input(*, *, k), 7, 7, /arithmetic, invalid = 0.0))

            mean_image[band - 1] = convolve(srcband, box)
            mean_image[band - 1][np.isnan(srcband)] = np.nan

            # print mean_image[band-1]

            # 2nd STEP: divide the original band by the filtered band
            # (ratio(*, *, k) = input(*, *, k)/mean_image(*, *, k))
            # print mean_image[band-1]
            # print 'Ratio'
            ratio[band - 1] = srcband / mean_image[band - 1]
            # print ratio[band-1]
            # print ratio[band-1]

            # 3rd step: calculate the average image from the ratio images
            # (mean_ratio(*, *) = mean(ratio(*, *, *), dimension = 3))
            # print 'Done'
        mean_ratio = np.nanmean(ratio, axis=0)
        # print mean_ratio.shape
        # 4th/final STEP works also on a band-by-band basis
        # for l = 0, (num_layers-1) do begin
        # output(*, *, l) = mean_image(*, *, l) * mean_ratio(*, *)
        # endfor
        print 'Output'
        for band in range(band_count):

            band += 1
            print 'BandNummer:'
            print band
            print asctime()
            output[band - 1] = mean_image[band - 1] * mean_ratio

            if i != 0:
                output[band - 1] = np.delete(output[band - 1], range((kernel - 1) / 2), axis=0)
                # print 'Output:'
                # print output[1].shape

        for band in range(band_count):
            band += 1
            # print band
            # print i
            yoffset = int((maxlines - (kernel - 1)) * i + (kernel - 1) / 2)
            if i == 0:
                yoffset = 0
            maskout = out.GetRasterBand(band)
            maskout.WriteArray(output[band - 1], 0, yoffset)
            maskout.FlushCache()
        output = [[] for i in range(band_count)]
    out.SetGeoTransform(data.raster.GetGeoTransform())
    out.SetProjection(data.raster.GetProjection())

    out = None

    print "End:", asctime()
    print "##########################################"


def quegan_cube(infile, outfile, kernel, time_kernel, is_list=False, nodata=0.0, max_memory=10000):
    print "##########################################"
    print "Start:", asctime()

    if not is_list:
        data = raster.Raster(infile)
        band_count = data.bands
    else:
        band_count = len(infile)
        data = raster.Raster(infile[0])

    # Import Layerstack

    print "RASTER BAND COUNT: ", band_count
    print "Cols, Rows:"
    print data.cols, data.rows
    mean_image = [[] for i in range(band_count)]
    # print mean_image
    ratio = [[] for i in range(band_count)]
    output = [[] for i in range(band_count)]
    maxlines = max_memory // (data.cols * band_count * 3 * 4 / 1048576.)
    print int(round(data.rows / (maxlines - kernel + 1) + 0.5))
    box = Box2DKernel(kernel)
    driver = gdal.GetDriverByName("ENVI")
    out = driver.Create(outfile, data.cols, data.rows, data.bands, GDT_Float32)
    for i in range(int(round(data.rows / (maxlines - kernel + 1) + 0.5))):
        # Loop through Bands and rows
        print 'Streifennummer: '
        print i
        print asctime()

        for band in range(band_count):
            band += 1

            if is_list:
                data = raster.Raster(infile[band])
                if data.nodata:
                    nodata = data.nodata

            print 'BandNummer:'
            print band
            print asctime()
            srcband = np.split(data.matrix(band), [(maxlines - kernel + 1) * i, (maxlines - kernel + 1) * i + maxlines], axis=0)[1]
            srcband[srcband == nodata] = np.nan

            # print srcband.shape
            # print 'Test'
            # mean_image[band-1] = np.zeros(srcband)
            # print 'Test'
            # 1st STEP: apply a mean filter (in this case 7x7) to each separate image to estimate local mean
            # (mean_image(*, *, k) = mean_filter(input(*, *, k), 7, 7, /arithmetic, invalid = 0.0))

            mean_image[band - 1] = convolve(srcband, box)
            mean_image[band - 1][np.isnan(srcband)] = np.nan

            # print mean_image[band-1]

            # 2nd STEP: divide the original band by the filtered band
            # (ratio(*, *, k) = input(*, *, k)/mean_image(*, *, k))
            # print mean_image[band-1]
            # print 'Ratio'
            ratio[band - 1] = srcband / mean_image[band - 1]
            # print ratio[band-1]
            # print ratio[band-1]


            # 3rd step: calculate the average image from the ratio images
            # (mean_ratio(*, *) = mean(ratio(*, *, *), dimension = 3))
            # print 'Done'
            # mean_ratio = np.nanmean(ratio, axis=0)
            # print mean_ratio.shape
            # 4th/final STEP works also on a band-by-band basis
            # for l = 0, (num_layers-1) do begin
            # output(*, *, l) = mean_image(*, *, l) * mean_ratio(*, *)
        # endfor
        print 'Output'
        for band in range(band_count):
            mean_ratio = ratio[band - 1]
            num = 1
            for time in range(1, (int(time_kernel) - 1) / 2):
                if band + time in range(band_count):
                    mean_ratio += ratio[band + time]
                    num += 1
                if band - time in range(band_count):
                    mean_ratio += ratio[band - time]
                    num += 1
            mean_ratio /= num
            band += 1
            print 'BandNummer:'
            print band
            print asctime()
            output[band - 1] = mean_image[band - 1] * mean_ratio

            if i != 0:
                output[band - 1] = np.delete(output[band - 1], range((kernel - 1) / 2), axis=0)
                # print 'Output:'
                # print output[1].shape

        for band in range(band_count):
            band += 1
            # print band
            # print i
            yoffset = int((maxlines - (kernel - 1)) * i + (kernel - 1) / 2)
            if i == 0:
                yoffset = 0
            maskout = out.GetRasterBand(band)
            maskout.WriteArray(output[band - 1], 0, yoffset)
            maskout.FlushCache()
        output = [[] for i in range(band_count)]
    out.SetGeoTransform(data.raster.GetGeoTransform())
    out.SetProjection(data.raster.GetProjection())

    out = None

    print "End:", asctime()
    print "##########################################"

def emd(infile, outfile,nodata=0.0, max_memory=10000):
    print "##########################################"
    print "Start:", asctime()

    if infile is not list:
        data = raster.Raster(infile)
        band_count = data.bands
    else:
       band_count = len(infile)
       data=raster.Raster(infile[0])
    decomposition=[]
    print data.dim
    things =  data.raster.ReadAsArray()
    things[things==nodata] = np.nan
    output = np.zeros_like(things)
    print things.shape
    decomposition = np.empty_like(things)
    print data.rows, data.cols, data.bands
    non_filtered = 0
    for row in range(data.rows):
        for col in range(data.cols):
            print "Row and column:" row, col
            timeseries = things[:, row, col]

            timeseries_masked = timeseries[np.logical_not(np.isnan(timeseries))]

            deco = hht.EMD(timeseries_masked)
            try:
                imfs = deco.decompose()
            except ValueError:
                imfs=[timeseries_masked]
                non_filtered+=1
            except TypeError:
                imfs=[timeseries_masked]
                non_filtered+=1
            num_nan =0

            for i in range(len(timeseries)):
                if np.isnan(timeseries[i]):
                    num_nan +=1
                    output[i,row,col] = np.nan
                else:
                    output[i,row,col] = timeseries[i]-imfs[0][i-num_nan]

            #hhtvis.plot_imfs(timeseries_masked, np.asarray(range(len(timeseries_masked))), imfs)
    #return hht.EMD(timeseries).decompose()
    print out
    print non_filtered
    driver = gdal.GetDriverByName("ENVI")
    out = driver.Create(outfile, data.cols, data.rows, data.bands, GDT_Float32)


    for band in range(band_count):
        band+=1
           # print band
            #print i
        maskout = out.GetRasterBand(band)
        maskout.WriteArray(output[band-1,:,:])
        maskout.FlushCache()
    out.SetGeoTransform(data.raster.GetGeoTransform())
    out.SetProjection(data.raster.GetProjection())

    out = None

    print "End:", asctime()
    print "##########################################"
    return None


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('input', help='The input file.')
    parser.add_argument('--kernel', '-k', default=3, type=float, help='The size of the moving window')
    parser.add_argument('--memory', '-m', default=10240, type=float, help='The size of the maximal ram')
    parser.add_argument('--nodata', '-n', default=0, type=float, help='The no data value')
    parser.add_argument('--islist', '-l', default=False, help='Is it a list?')
    parser.add_argument('--time', '-t', default=None, help='The number of scenes that should be averaged in the time axis. If it is not specified, all time steps are averaged.')
    args = vars(parser.parse_args())

    infile = args['input']
    kernel = args['kernel']
    outfile = infile + 'mtf' + str(kernel)

    if args['time']:
        quegan_cube(infile, outfile, kernel, args['time'], False, args['nodata'], args['memory'])
    else:
        quegan(infile, outfile, kernel, False, args['nodata'], args['memory'])
