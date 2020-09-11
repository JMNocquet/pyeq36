#!/usr/bin/env python

###################################################################
# SCRIPT    : pyeq_plot_kinematic_results.py
# AUTHOR    : JM NOCQUET
# DATE      : June 2020
###################################################################

###################################################################
# MODULES IMPORT
###################################################################

import sys
import os
import argparse
import pyeq.lib.plot
import pyacs.lib.shapefile
import pickle
from colors import red
import glob


###################################################################
# PARSE ARGUMENT LINE
###################################################################

prog_info="Makes model, stf, time series plots from from pyeq_kinematic_inversion.py directory results"
prog_epilog="J.-M. Nocquet (Geoazur-UCA-IRD-CNRS-OCA) - April 2020"


parser = argparse.ArgumentParser(description=prog_info,epilog=prog_epilog,prefix_chars='-')
parser.add_argument('-odir', action='store', type=str, dest='odir',required=True,help='output directory')
parser.add_argument('-shp_line', action='append', type=str, dest='shp_line',default=[],help='line shapefile. Can be repeated')
parser.add_argument('-shp_poly', action='append', type=str, dest='shp_poly',default=[],help='poly shapefile. Can be repeated')
parser.add_argument('-verbose', '-v', action='count',default=0,help='verbose mode')

args = parser.parse_args()

###################################################################
# LOAD MODEL PCK
###################################################################

pck = args.odir+'/npy/model.pck'
try:
    print("-- Loading %s (%.2f Gb) " % ( pck , os.path.getsize( pck ) /1024 / 1024 / 1024 ) )
    with open( pck, "rb") as f:
        model = pickle.load( f )
    f.close()
    print("-- model object loaded.")
except:
    print( red("[PYEQ ERROR] Could not load: %s " % ( pck ) ))
    sys.exit()

###################################################################
# UPDATE VERBOSE & ODIR
###################################################################

model.verbose = args.verbose
model.odir = args.odir

###########################################################################
# MAKE PLOT FOR TIME SERIES
###########################################################################

print("-- plotting time series")
pyeq.lib.plot.plot_time_series(model)

###########################################################################
# GENERATE SHAPEFILES OF CUMULATIVE & RATE SLIP FOR MODEL DISPLAY IN QGIS
###########################################################################

print('-- writing GMT and shapefiles for visualization in QGIS')
one_degree = 111.1
TRIANGLE = True

lslip_dat = glob.glob(model.odir + "/slip/cumulative/*_cumulative_slip.dat")
pyeq.lib.plot.model2shp_gmt(model.geometry, 'tde', lslip_dat, out_dir_shp=model.odir + '/shapefile/slip_cumulative',
                            out_dir_gmt=model.odir + '/gmt/slip_cumulative', verbose=model.verbose)

lslip_dat = glob.glob(model.odir + "/slip/rate/*_slip_rate.dat")
pyeq.lib.plot.model2shp_gmt(model.geometry, 'tde', lslip_dat, out_dir_shp=model.odir + '/shapefile/slip_rate',
                            out_dir_gmt=model.odir + '/gmt/slip_rate', verbose=model.verbose)


###########################################################################
# GENERATE SHAPEFILES OF CUMULATIVE & RATE GPS DISPLACEMENTS FOR DISPLAY IN QGIS
###########################################################################

ldisp_dat = glob.glob(model.odir + "/displacement/cumulative/model/*disp.dat")
for disp in ldisp_dat:
    shp = model.odir + "/shapefile/disp_cumulative/"+disp.split('/')[-1].split('.')[0]
    pyacs.lib.shapefile.psvelo_to_shapefile(disp, shp, verbose=False)

###################################################################
# UPDATE MODEL WITH EXTERNAL SHAPEFILE FOR PLOT
###################################################################
model.external_shapefile_line = args.shp_line
model.external_shapefile_poly = args.shp_poly

###################################################################
# MAKE PLOT
###################################################################
#
pyeq.lib.plot.make_plot( model )
