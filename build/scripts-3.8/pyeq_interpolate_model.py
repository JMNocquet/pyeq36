#!/usr/local/geodesy/anaconda38/bin/python

###################################################################
# SCRIPT    : pyeq_interpolate_model.py
# AUTHOR    : JM NOCQUET
# DATE      : April 2020
###################################################################

###################################################################
# MODULES IMPORT
###################################################################

import sys
import os
import argparse
import numpy as np
import pickle
from colors import red

import pyeq.plot

###################################################################
# PARSE ARGUMENT LINE
###################################################################

prog_info="Interpolate models from pyeq_kinematic_inversion.py directory results"
prog_epilog="J.-M. Nocquet (Geoazur-UCA-IRD-CNRS-OCA) - April 2020"


parser = argparse.ArgumentParser(description=prog_info,epilog=prog_epilog,prefix_chars='-')
parser.add_argument('-odir', action='store', type=str, dest='odir',required=True,help='output directory')
parser.add_argument('-geometry', action='store', type=str, dest='geometry',required=True,help='geometry file as npy')
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
# LOAD GEOMETRY
###################################################################

try:
    geometry = np.load( args.geometry )
except:
    print( red("[PYEQ ERROR] Could not load: %s " % ( args.geometry ) ))
    sys.exit()

print("-- interpolate %d faults on %d faults" % ( model.geometry.shape[0] , geometry.shape[0] ) )
    
###################################################################
# MAKE INTERPOLATION
###################################################################

pyeq.plot.interpolate_model(model, geometry)
