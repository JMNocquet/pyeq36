#!/usr/local/geodesy/anaconda38/bin/python

###################################################################
# SCRIPT    : pyeq_plot_kinematic_results.py
# AUTHOR    : JM NOCQUET
# DATE      : June 2020
###################################################################

###################################################################
# MODULES IMPORT
###################################################################

import os
import argparse
import pyeq.plot
import pyeq.lib.objects.pyeq_model
import pyacs.lib.shapefile
import glob
import pyeq.plot
from str2bool import str2bool
from os import mkdir, path

import logging
import pyeq.message.message as MESSAGE
import pyeq.message.verbose_message as VERBOSE
import pyeq.message.error as ERROR
import pyeq.message.warning as WARNING
import pyeq.message.debug_message as DEBUG

#TODO STress

###################################################################
# PARSE ARGUMENT LINE
###################################################################

SETTINGS = pyeq.plot.plot_settings()

prog_info="Makes various plots from from pyeq_kinematic_inversion.py directory results"
prog_epilog="J.-M. Nocquet (Geoazur-UCA-IRD-CNRS-OCA) - April 2020"

parser = argparse.ArgumentParser(description=prog_info,epilog=prog_epilog,prefix_chars='-')
parser.add_argument('-odir', action='store', type=str, dest='odir',default=None,help='output directory')
parser.add_argument('-conf', action='store', type=str, dest='conf',default=None,help='output directory')
parser.add_argument('-defaults', action='count',default=0,help='print defaults. Useful to make a template')
parser.add_argument('-verbose', '-v', action='count',default=0,help='verbose mode')
parser.add_argument('-debug', action='count',default=0,help='debug mode')

# parser.add_argument('-shp_line', action='append', type=str, dest='shp_line',default=[],help='line shapefile. Can be repeated')
# parser.add_argument('-shp_poly', action='append', type=str, dest='shp_poly',default=[],help='poly shapefile. Can be repeated')
# parser.add_argument('-all', '-a', action='count',default=0,help='make all plots')
# parser.add_argument('-time_series', '-ts', action='count',default=0,help='make time series plots')
# parser.add_argument('-cum', action='count',default=0,help='make cumulative slip plots')
# parser.add_argument('-ccum', action='count',default=0,help='make contour cumulative slip plots')
# parser.add_argument('-rate', action='count',default=0,help='make slip rate plots')
# parser.add_argument('-crate', action='count',default=0,help='make contour rate slip plots')
# parser.add_argument('-obs', action='count',default=0,help='make cumulative slip plots with obs/model GPS arrows')
# parser.add_argument('-res', action='count',default=0,help='make cumulative slip plots with residuals GPS arrows')
# parser.add_argument('-stf', action='count',default=0,help='make stf/cstf plots')

args = parser.parse_args()

# verbose & debug
logging.getLogger("my_logger").setLevel(logging.WARNING)

if args.verbose>0:
    verbose=True
    logging.getLogger("my_logger").setLevel(logging.INFO)
else:
    verbose=False
if args.debug>0:
    debug=True
    logging.getLogger("my_logger").setLevel(logging.DEBUG)
else:debug=False

if args.defaults:
    SETTINGS.print()
    exit()

if args.odir is None:
    parser.print_help()
    exit()


VERBOSE("Current settings for plots")
if verbose: SETTINGS.print()

###################################################################
# LOAD MODEL PCK
###################################################################

pck = args.odir+'/npy/model.mpck'

model = pyeq.lib.objects.pyeq_model.load( pck )
model.verbose = args.verbose
model.debug = args.debug

if debug:model.print_info()

###################################################################
# READ PLOT CONF FILE IF PROVIDED
###################################################################

if args.conf is not None:
    SETTINGS.load(args.conf)

model.plot_settings = SETTINGS

###################################################################
# MAKE PLOT DIRECTORIES
###################################################################

model.odir = args.odir

###################################################################
if not path.exists(model.odir + '/plots'): mkdir(model.odir + '/plots')
if not path.exists(model.odir + '/plots/map'): mkdir(model.odir + '/plots/map')
if not path.exists(model.odir + '/shapefile'): mkdir(model.odir + '/shapefile')

###########################################################################
# MAKE PLOT FOR TIME SERIES
###########################################################################
if str2bool(SETTINGS.time_series):
    MESSAGE("plotting time series")
    pyeq.plot.plot_time_series(model)


###########################################################################
# GENERATE SHAPEFILES OF CUMULATIVE & RATE SLIP FOR MODEL DISPLAY IN QGIS
###########################################################################

# change 05/04/2021 - always generate shapefiles
#if str2bool(model.plot_settings.cum_slip) or str2bool(model.plot_settings.ccum_slip) \
#        or str2bool(model.plot_settings.rate) or str2bool(model.plot_settings.crate):

MESSAGE("writing cumulative slip and slip rate GMT and shapefiles for visualization in QGIS")
one_degree = 111.1
TRIANGLE = True

lslip_dat = glob.glob(model.odir + "/slip/cumulative/*_cumulative_slip.dat")
pyeq.plot.model2shp_gmt(model.geometry, 'tde', lslip_dat, out_dir_shp=args.odir + '/shapefile/slip_cumulative',
                        out_dir_gmt=model.odir + '/gmt/slip_cumulative', verbose=model.debug)

lslip_dat = glob.glob( model.odir + "/slip/rate/*_slip_rate.dat")
pyeq.plot.model2shp_gmt(model.geometry, 'tde', lslip_dat, out_dir_shp=args.odir + '/shapefile/slip_rate',
                        out_dir_gmt=model.odir + '/gmt/slip_rate', verbose=model.debug)

MESSAGE("writing spatial resolution map as shapefile")
pyeq.plot.model2shp_gmt(model.geometry, 'tde', [model.odir + "/info/spatial_resolution.dat"], out_dir_shp=args.odir + '/shapefile',
                        out_dir_gmt=model.odir + '/gmt', verbose=model.debug)

###########################################################################
# GENERATE SHAPEFILES OF CUMULATIVE & RATE GPS DISPLACEMENTS FOR DISPLAY IN QGIS
###########################################################################

MESSAGE("writing GPS displacements GMT and shapefiles for visualization in QGIS")

ldisp_dat = glob.glob(model.odir + "/displacement/cumulative/model/*disp.dat")
for disp in ldisp_dat:
    shp = model.odir + "/shapefile/disp_cumulative/"+disp.split('/')[-1].split('.')[0]
    pyacs.lib.shapefile.psvelo_to_shapefile(disp, shp, verbose=False)

###########################################################################
# GENERATE SHAPEFILE FOR SLIP DIRECTION
###########################################################################

# unit slip vector to shapefile
shp = model.odir + "/shapefile/geometry/slip_dir_en"
disp = model.odir + '/info/slip_dir_en.dat'

pyacs.lib.shapefile.psvelo_to_shapefile(disp, shp, verbose=False)

###################################################################
# UPDATE MODEL WITH EXTERNAL SHAPEFILE FOR PLOT
###################################################################
model.external_shapefile_line = model.plot_settings.shp_line
model.external_shapefile_poly = model.plot_settings.shp_poly

###########################################################################
# MAKE PLOT FOR MODELS
###########################################################################


MESSAGE("plotting model")
# if model.interseismic !=0:
#    pyeq.lib.plot.plot_model_interseismic( model )
# else:
pyeq.plot.plot_model(model, interpolation=False)

###########################################################################
# MAKE PLOT FOR STF
###########################################################################

if str2bool(model.plot_settings.stf):
    MESSAGE("making plots for stf/cstf results in %s" % (model.odir + '/plots/stf'))
    pyeq.plot.plot_stf(model)
