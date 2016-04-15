#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import math
import os
import os.path
from os.path import join
import shutil
import subprocess

from gamma.util import ISPPar

GEOCODE_SUBDIR = 'geocode'
COREGISTRATION_SUBDIR = 'coregistration'
ENHANCEMENT_SUBDIR = 'enhancement'
PRODUCT_SUBDIR = 'products'


class SubprocessWrapper(object):
    def __init__(self):
        self.stdin = None
        self.stdout = None
        self.stderr = None

    def __getattr__(self, name):
        return lambda *args: self.gamma(name, *args)

    def gamma(self, command, *args):
        args = [str(command)] + [str(a) for a in args]  # Convert all arguments into strings
        return subprocess.check_call(args, stdin=self.stdin, stdout=self.stdout, stderr=self.stderr)


def recreate_directory(filename):
    try:
        os.mkdir(filename)
    except OSError:
        shutil.rmtree(filename)
        os.mkdir(filename)


def main(args):
    gamma = SubprocessWrapper()
    gamma.stdout = open('logfile', 'w')
    gamma.stderr = subprocess.STDOUT  # Redirect stderr to stdout
    try:
        #
        # Generate multi-looked images
        #
        active_par = ISPPar(args.SLC1_par)
        az_factor = int(round(args.target_res / active_par.azimuth_pixel_spacing, ndigits=0))
        rg_factor = int(round(args.target_res /
                              (active_par.range_pixel_spacing / math.sin(math.radians(active_par.incidence_angle))),
                              ndigits=0))
        acquisition_date = '{:04d}{:02d}{:02d}'.format(active_par.date[0], active_par.date[1], active_par.date[2])
        gamma.multi_look(args.SLC1, args.SLC1_par, 'active.mli', 'active.mli.par', rg_factor, az_factor)
        gamma.multi_look(args.SLC2, args.SLC2_par, 'passive.mli', 'passive.mli.par', rg_factor, az_factor)
        active_par = ISPPar('active.mli.par')
        if args.quicklooks:
            gamma.raspwr('active.mli', active_par.range_samples, 1, active_par.azimuth_lines, 1, 1, 1., .35, 1,
                         'active.mli.bmp')
        #
        # Transform the DEM into radar-doppler coordinates (with refinement step)
        #
        recreate_directory(GEOCODE_SUBDIR)
        dem_par = ISPPar(args.DEM_par)
        oversampling_factor = dem_par.post_east / float(args.target_res)
        gamma.gc_map('active.mli.par', '-', args.DEM_par, args.DEM,
                     join(GEOCODE_SUBDIR, 'ref.map.dem_par'), join(GEOCODE_SUBDIR, 'ref.map.dem'),
                     join(GEOCODE_SUBDIR, 'ref.rough.map_to_rdc'), oversampling_factor, oversampling_factor,
                     join(GEOCODE_SUBDIR, 'ref.map.sim_sar'), '-', '-', join(GEOCODE_SUBDIR, 'ref.map.inc'),
                     '-', join(GEOCODE_SUBDIR, 'ref.map.pix'), join(GEOCODE_SUBDIR, 'ref.map.ls_map'))
        os.chdir(GEOCODE_SUBDIR)
        dem_par = ISPPar('ref.map.dem_par')
        gamma.geocode('ref.rough.map_to_rdc', 'ref.map.sim_sar', dem_par.width, 'ref.rdc.sim_sar',
                      active_par.range_samples, active_par.azimuth_lines, 2)
        gamma.create_diff_par('../active.mli.par', '-', 'diff_par', 1, 0)
        if args.offs_init:
            try:
                gamma.init_offsetm('../active.mli', 'ref.rdc.sim_sar', 'diff_par')
            except subprocess.CalledProcessError:
                pass
        for size in args.window_sizes:
            # Calculate number of windows required to meet the overlap factor
            rg_count = active_par.range_samples * args.window_overlap // size
            az_count = active_par.azimuth_lines * args.window_overlap // size
            # Estimate the offsets
            gamma.offset_pwrm('../active.mli', 'ref.rdc.sim_sar', 'diff_par', 'offs', 'snr', size, size, 'offsets', 1,
                              rg_count, az_count, 7., 0)
            gamma.offset_fitm('offs', 'snr', 'diff_par', 'coffs', 'coffsets', 7., 3)
        gamma.gc_map_fine('ref.rough.map_to_rdc', dem_par.width, 'diff_par', 'ref.fine.map_to_rdc', 1)
        gamma.geocode('ref.fine.map_to_rdc', 'ref.map.dem', dem_par.width, 'ref.rdc.dem', active_par.range_samples,
                      active_par.azimuth_lines, 2)
        os.chdir('..')
        #
        # Create interferometric phase and derive a differential interferogram
        #
        gamma.create_offset(args.SLC2_par, args.SLC1_par, 'offset_par', 1, rg_factor, az_factor, 0)
        gamma.phase_sim_orb(args.SLC1_par, args.SLC2_par, 'offset_par',
                            join(GEOCODE_SUBDIR, 'ref.rdc.dem'), 'tdm.ph_sim_orb', args.SLC1_par, '-', '-', 0)
        gamma.SLC_diff_intf(args.SLC1, args.SLC2, args.SLC1_par, args.SLC2_par, 'offset_par', 'tdm.ph_sim_orb',
                            'tdm.diff', rg_factor, az_factor, 1, 0, .25)
        if args.quicklooks:
            gamma.rasmph_pwr24('tdm.diff', 'active.mli', active_par.range_samples, 1, 1, active_par.azimuth_lines,
                               1, 1, 1., .35, 1, 'tdm.diff.bmp')
        #
        # Unwrap the phase with MCF
        #
        gamma.cc_ad('tdm.diff', 'active.mli', 'passive.mli', '-', '-', 'tdm.cc_ad', active_par.range_samples)
        if args.quicklooks:
            gamma.rascc('tdm.cc_ad', 'active.mli', active_par.range_samples, 1, 1, active_par.azimuth_lines,
                        1, 1, .1, .9, 1., .35, 1, 'tdm.cc_ad.bmp')
        gamma.rascc_mask('tdm.cc_ad', 'active.mli', active_par.range_samples, 1, 1, active_par.azimuth_lines, 1, 1, .3,
                         0., 0., 1., 1., .35, 1, 'tdm.cc_ad.ras')
        gamma.mcf('tdm.diff', 'tdm.cc_ad', '-', 'tdm.diff.unw', active_par.range_samples, 0, 0, 0, '-',
                  '-', 1, 1, '-', active_par.range_samples // 2, active_par.azimuth_lines // 2, 0)
        if args.quicklooks:
            gamma.rasrmg('tdm.diff.unw', 'active.mli', active_par.range_samples, 1, 1, active_par.azimuth_lines, 1,
                         1, .5, 1., .35, 0., 1, 'tdm.diff.unw.bmp')
        #
        # Remove large scale tilts
        #
        gamma.create_diff_par('active.mli.par', '-', 'diff_par', 1, 0)
        gamma.quad_fit('tdm.diff.unw', 'diff_par', 16, 16, 'tdm.cc_ad.ras', '-', 3)
        gamma.quad_sub('tdm.diff.unw', 'diff_par', 'tdm.diff.unw.shifted', 0, 0)
        #
        # Convert unwrapped phases to relative heights
        #
        gamma.dh_map_orb(args.SLC1_par, args.SLC2_par, 'offset_par', join(GEOCODE_SUBDIR, 'ref.rdc.dem'),
                         'tdm.diff.unw.shifted', 'tdm.dpdh', 'tdm.dh', args.SLC1_par, 0)
        #
        #
        #
        recreate_directory(ENHANCEMENT_SUBDIR)
        if not args.master is None:
            master_dir = args.master
            recreate_directory(COREGISTRATION_SUBDIR)
            gamma.rdc_trans(join(master_dir, 'active.mli.par'),
                            join(master_dir, GEOCODE_SUBDIR, 'ref.rdc.dem'), 'active.mli.par',
                            join(COREGISTRATION_SUBDIR, 'lut.rough.mas_to_slv'))
            master_par = ISPPar(join(master_dir, 'active.mli.par'))
            dem_par = ISPPar(join(master_dir, GEOCODE_SUBDIR, 'ref.map.dem_par'))
            os.chdir(COREGISTRATION_SUBDIR)
            # Coregistrate the active intensity image to the master scene
            gamma.geocode('lut.rough.mas_to_slv', join(master_dir, 'active.mli'), master_par.range_samples,
                          'master.rmli', active_par.range_samples, active_par.azimuth_lines, 2, 0)
            gamma.create_diff_par('../active.mli.par', '-', 'diff_par', 1, 0)
            try:
                gamma.init_offsetm('master.rmli', '../active.mli', 'diff_par', 1, 1)
            except subprocess.CalledProcessError:
                pass
            gamma.offset_pwrm('master.rmli', '../active.mli', 'diff_par', 'offs', 'snr', 256, 256, 'offsets', 2, 32, 32, 7., 0)
            gamma.offset_fitm('offs', 'snr', 'diff_par', 'coffs', 'coffsets', 7., 4)
            gamma.gc_map_fine('lut.rough.mas_to_slv', master_par.range_samples, 'diff_par', 'lut.fine.mas_to_slv')
            # Resample intensity, coherence and heights into the master RDC geometry
            gamma.geocode_back('../active.mli', active_par.range_samples, 'lut.fine.mas_to_slv', 'active.mli', master_par.range_samples, master_par.azimuth_lines, 1, 0)
            gamma.geocode_back('../tdm.cc_ad', active_par.range_samples, 'lut.fine.mas_to_slv', 'tdm.cc_ad', master_par.range_samples, master_par.azimuth_lines, 1, 0)
            gamma.geocode_back('../tdm.dh', active_par.range_samples, 'lut.fine.mas_to_slv', 'tdm.dh', master_par.range_samples, master_par.azimuth_lines, 1, 0)
            os.chdir('..')
            # Resample all products from the RDC in the MAP geometry
            gamma.geocode_back(join(COREGISTRATION_SUBDIR, 'active.mli'), active_par.range_samples,
                               join(master_dir, GEOCODE_SUBDIR, 'ref.fine.map_to_rdc'),
                               join(ENHANCEMENT_SUBDIR, 'active.mli0'), dem_par.width, dem_par.nlines, 1, 0)
            gamma.geocode_back(join(COREGISTRATION_SUBDIR, 'tdm.cc_ad'), active_par.range_samples,
                               join(master_dir, GEOCODE_SUBDIR, 'ref.fine.map_to_rdc'),
                               join(ENHANCEMENT_SUBDIR, 'tdm.cc_ad0'), dem_par.width, dem_par.nlines, 1, 0)
            gamma.geocode_back(join(COREGISTRATION_SUBDIR, 'tdm.dh'), active_par.range_samples,
                               join(master_dir, GEOCODE_SUBDIR, 'ref.fine.map_to_rdc'),
                               join(ENHANCEMENT_SUBDIR, 'tdm.dh'), dem_par.width, dem_par.nlines, 1, 0)
        else:
            master_dir = os.getcwd()
            master_par = active_par
            gamma.geocode_back('active.mli', active_par.range_samples,
                               join(GEOCODE_SUBDIR, 'ref.fine.map_to_rdc'),
                               join(ENHANCEMENT_SUBDIR, 'active.mli0'), dem_par.width, dem_par.nlines, 1, 0)
            gamma.geocode_back('tdm.cc_ad', active_par.range_samples,
                               join(GEOCODE_SUBDIR, 'ref.fine.map_to_rdc'),
                               join(ENHANCEMENT_SUBDIR, 'tdm.cc_ad0'), dem_par.width, dem_par.nlines, 1, 0)
            gamma.geocode_back('tdm.dh', active_par.range_samples,
                               join(GEOCODE_SUBDIR, 'ref.fine.map_to_rdc'),
                               join(ENHANCEMENT_SUBDIR, 'tdm.dh'), dem_par.width, dem_par.nlines, 1, 0)
        #
        # Apply quality enhancements to the surface model
        #
        dem_par = ISPPar(join(master_dir, GEOCODE_SUBDIR, 'ref.map.dem_par'))
        gamma.float_math(join(master_dir, GEOCODE_SUBDIR, 'ref.map.dem'),
                         join(ENHANCEMENT_SUBDIR, 'tdm.dh'), join(ENHANCEMENT_SUBDIR, 'tdm.hgt0'),
                         dem_par.width, 0)
        os.chdir(ENHANCEMENT_SUBDIR)
        # Mask areas of low coherence
        gamma.rascc_mask('tdm.cc_ad0', 'active.mli0', dem_par.width, 1, 1, dem_par.nlines, 1, 1,
                         args.coh_thres, 0., 0., 1., 1., .35, 1, 'tdm.cc_ad0.ras')
        gamma.mask_class('tdm.cc_ad0.ras', 'tdm.hgt0', 'tdm.hgt1', 0, 1, 1, 1, 0, 0.)
        # Mask outliers
        gamma.interp_ad('tdm.hgt0', 'tdm.hgt0.filt', dem_par.width, 4, 9, 25, 2, 2, 0)
        gamma.lin_comb(2, 'tdm.hgt1', 'tdm.hgt0.filt', 0., 1., -1., 'tdm.dh1', dem_par.width, 1, dem_par.nlines,
                       1, 1, 0)
        gamma.single_class_mapping(1, 'tdm.dh1', -10., 10., 'tdm.dh1a.ras', dem_par.width, 1, dem_par.nlines, 1, 1)
        gamma.mask_class('tdm.dh1a.ras', 'tdm.hgt1', 'tdm.hgt2', 0, 1, 1, 1, 0, 0.)
        # Interpolate small holes
        gamma.interp_ad('tdm.hgt2', 'tdm.hgt3', dem_par.width, 4, 8, 2, 2, 1)
        # Apply spatial filtering to reduce noise
        gamma.interp_ad('tdm.hgt3', 'tdm.hgt', dem_par.width, 1, 4, 15, 2, 2, 0)
        # Mask areas with low pixel normalization areas to avoid interpolation artifacts
        gamma.single_class_mapping(1, join(master_dir, GEOCODE_SUBDIR, 'ref.map.pix'), 0.,
                                   args.pix_thres, 'pixnorm_mask.ras', dem_par.width, 1, dem_par.nlines)
        gamma.mask_class('pixnorm_mask.ras', 'active.mli0', 'active.mli', 0, 1, -1, 1, 0, 0.)
        gamma.mask_class('pixnorm_mask.ras', 'tdm.cc_ad0', 'tdm.cc_ad', 0, 1, -1, 1, 0, 0.)
        gamma.linear_to_dB('active.mli0', 'active.mli0.db', dem_par.width, 0, -99.)
        gamma.linear_to_dB('active.mli', 'active.mli.db', dem_par.width, 0, -99.)
        os.chdir('..')
        #
        # Fill the final product directory
        #
        recreate_directory(PRODUCT_SUBDIR)
        gamma.data2geotiff(join(master_dir, GEOCODE_SUBDIR, 'ref.map.dem_par'),
                           join(ENHANCEMENT_SUBDIR, 'active.mli0.db'), 2,
                           join(PRODUCT_SUBDIR, acquisition_date + '_int_db.tif'), -99.)
        gamma.data2geotiff(join(master_dir, GEOCODE_SUBDIR, 'ref.map.dem_par'),
                           join(ENHANCEMENT_SUBDIR, 'active.mli.db'), 2,
                           join(PRODUCT_SUBDIR, acquisition_date + '_int_mask_db.tif'), -99.)
        gamma.data2geotiff(join(master_dir, GEOCODE_SUBDIR, 'ref.map.dem_par'), join(ENHANCEMENT_SUBDIR, 'tdm.cc_ad0'),
                           2, join(PRODUCT_SUBDIR, acquisition_date + '_coh.tif'))
        gamma.data2geotiff(join(master_dir, GEOCODE_SUBDIR, 'ref.map.dem_par'), join(ENHANCEMENT_SUBDIR, 'tdm.cc_ad'),
                           2, join(PRODUCT_SUBDIR, acquisition_date + '_coh_mask.tif'))
        gamma.data2geotiff(join(master_dir, GEOCODE_SUBDIR, 'ref.map.dem_par'), join(ENHANCEMENT_SUBDIR, 'tdm.hgt0'),
                           2, join(PRODUCT_SUBDIR, acquisition_date + '_dsm.tif'))
        gamma.data2geotiff(join(master_dir, GEOCODE_SUBDIR, 'ref.map.dem_par'), join(ENHANCEMENT_SUBDIR, 'tdm.hgt'),
                           2, join(PRODUCT_SUBDIR, acquisition_date + '_dsm_mask.tif'))
        #
        # Generate additional products
        #
        if args.hillshade:
            gamma.gdaldem('hillshade', join(PRODUCT_SUBDIR, acquisition_date + '_dsm_mask.tif'),
                          join(PRODUCT_SUBDIR, acquisition_date + '_dsm_mask_hillshade.tif'))
        if args.slope:
            gamma.gdaldem('slope', join(PRODUCT_SUBDIR, acquisition_date + '_dsm_mask.tif'),
                          join(PRODUCT_SUBDIR, acquisition_date + '_dsm_mask_slope.tif'))
        if args.normalization:
            os.chdir(ENHANCEMENT_SUBDIR)
            gamma.sigma2gamma('active.mli', join(master_dir, GEOCODE_SUBDIR, 'ref.map.inc'), 'gamma0', dem_par.width)
            gamma.product('gamma0', join(master_dir, GEOCODE_SUBDIR, 'ref.map.pix'), 'gamma0_pix', dem_par.width, 1)
            os.remove('gamma0')
            gamma.lin_comb(1, 'gamma0_pix', 0., math.cos(math.radians(master_par.incidence_angle)), 'active.mli.norm',
                           dem_par.width)
            os.remove('gamma0_pix')
            gamma.linear_to_dB('active.mli.norm', 'active.mli.norm.db', dem_par.width, 0, -99.)
            os.chdir('..')
            gamma.data2geotiff(join(master_dir, GEOCODE_SUBDIR, 'ref.map.dem_par'),
                               join(ENHANCEMENT_SUBDIR, 'active.mli.norm.db'), 2,
                               join(PRODUCT_SUBDIR, acquisition_date + '_int_mask_norm_db.tif'), -99.)
    finally:
        gamma.stdout.close()


if __name__ == '__main__':
    # Parse command-line arguments
    parser = argparse.ArgumentParser(
        description='Generates a surface model from bistatic Tandem-X CoSSC experimental data. The approach used is '
                    'based on a reference elevation model of which the differential phase to the Tandem-X data is '
                    'estimated. For further information see the documentation distributed with this program. Note that '
                    'all intermediate products are written the current  directory and it\'s subdirectories.')
    parser.add_argument('-m', '--master', type=str, metavar='DIR',
                        help='Working directory of a scene used as master for the current one. All generated products '
                             'will be registrated to this scene\'s geometry. Be aware that the master scene should '
                             'have the same orbit direction and approximately the same extent.')
    parser.add_argument('-tr', '--target-resolution', type=float, default=12., metavar='RES', dest='target_res',
                        help='Targeted ground resolution in meters (default: 12 m).')
    opt = parser.add_argument_group(title='additional products')
    opt.add_argument('-ph', '--hillshade', default=False, const=True, action='store_const',
                     help='Creates a shaded relief map (hillshade) of the derived surface model.')
    opt.add_argument('-pn', '--normalization', default=False, const=True, action='store_const',
                     help='Apply a topographic normalization based on the reference DEM to the backscatter intensity '
                          'images.')
    opt.add_argument('-pq', '--quicklooks', default=False, const=True, action='store_const', dest='quicklooks',
                     help='Creates bitmap quicklooks for several intermediate products in SAR geometry.')
    opt.add_argument('-ps', '--slope', default=False, const=True, action='store_const',
                     help='Generate a slope map of the derived surface model.')
    mask = parser.add_argument_group(title='masking options')
    mask.add_argument('-mc', '--minimum-coherence', type=float, default=.3, metavar='COH', dest='coh_thres',
                      help='Coherence value which is considered the minimum value for valid height estimations, '
                           'pixels with a lower coherence will be masked (default: 0.3).')
    mask.add_argument('-mp', '--minimum-pixelarea', type=float, default=.3, metavar='PIX', dest='pix_thres',
                      help='Pixel area normalization factor threshold used for masking the backscatter intensity and '
                           'the coherence products (default: 0.3).')
    geo = parser.add_argument_group(title='geocoding options')
    geo.add_argument('-si', '--skip-initial-offset', default=True, const=False, action='store_const', dest='offs_init',
                     help='Skips the initial offset estimation during the geocoding procedure.')
    geo.add_argument('-wo', '--window-overlap', type=float, default=2.5, metavar='FAC', dest='window_overlap',
                     help='Overlapping factor applied to the scene\'s dimensions before calculating the number of '
                          'windows for offset estimation (default: 2.5).')
    geo.add_argument('-ws', '--window-size', type=int, nargs='+', default=[256], metavar='SIZE', dest='window_sizes',
                     help='Window size to be used during the geocoding procedure. To apply more iterations multiple '
                          'sizes must be supplied. The default is on iteration with 256 px wide windows.')
    parser.add_argument('SLC1', help='Single-look complex image (active).', type=str)
    parser.add_argument('SLC1_par', help='SLC parameter file (active).', type=str)
    parser.add_argument('SLC2', help='Single-look complex image (passive).', type=str)
    parser.add_argument('SLC2_par', help='SLC parameter file (passive).', type=str)
    parser.add_argument('DEM', help='DEM data file (reference).', type=str)
    parser.add_argument('DEM_par', help='DEM parameter file (reference).', type=str)
    main(parser.parse_args())

