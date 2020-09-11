#!/usr/bin/env python

###################################################################
# SCRIPT    : pyeq_plot_kinematic_results.py
# AUTHOR    : JM NOCQUET
# DATE      : April 2020
###################################################################

###################################################################
# MODULES IMPORT
###################################################################

import sys
import os
import argparse
import pyeq.lib.plot
import pickle
from colors import red

###################################################################
# PARSE ARGUMENT LINE
###################################################################

prog_info="Makes model, stf, time series plots from from pyeq_kinematic_inversion.py directory results"
prog_epilog="J.-M. Nocquet (Geoazur-UCA-IRD-CNRS-OCA) - April 2020"


parser = argparse.ArgumentParser(description=prog_info,epilog=prog_epilog,prefix_chars='-')
parser.add_argument('-odir', action='store', type=str, dest='odir',required=True,help='output directory')
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

###################################################################
# MAKE PLOT
###################################################################

pyeq.lib.plot.make_plot( model )
