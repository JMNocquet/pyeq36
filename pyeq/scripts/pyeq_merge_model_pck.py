#!/usr/bin/env python
'''
###################################################################
# SCRIPT    : pyeq_merge_model_model_pck.py
# AUTHOR    : Jean-Mathieu Nocquet
# DATE      : June 2020
# INPUT     :
# OUTPUT    :
# NOTE      : Major refactoring in October 2019 & March 2020
###################################################################
'''

###################################################################
# MODULES IMPORT
###################################################################

# GENERAL

import sys, os
import argparse
from argparse import RawTextHelpFormatter
import numpy as np
import pkg_resources
from time import time
import pyacs.lib.astrotime as at
import copy
import pickle
from colors import red

# PYACS

import pyacs.lib.glinalg
import pyeq.lib.date
import pyeq.lib.regularization
import pyacs.lib.utils

# PYEQ
import pyeq.lib.forward_model
import pyeq.lib.objects
import pyeq.lib.log
import pyeq.lib.conf
import pyeq.lib.elastic_tensor
import pyeq.lib.green_tensor
import pyeq.lib.gps_time_series
import pyeq.lib.obs_tensor.set_zero_at_first_obs
import pyeq.lib.make_inversion

###################################################################
# PARSE COMMAND LINE
###################################################################

prog_info="Merges pyeq_kinematic_inversion.py results"
prog_epilog="J.-M. Nocquet (Geoazur-UCA-IRD-CNRS-OCA) - June 2020"

parser = argparse.ArgumentParser(description=prog_info,epilog=prog_epilog,prefix_chars='-')
parser.add_argument('-odir', action='store', type=str, dest='odir',required=True,help='output directory')
parser.add_argument('-ipck', action='store', type=str, dest='ipck',required=True,help='list of model pck files')
parser.add_argument('-verbose', '-v', action='count',default=0,help='verbose mode')

args = parser.parse_args()

###################################################################
# INIT MODEL
# MODEL WILL INCLUDE ALL INFORMATION
###################################################################

merged_model = pyeq.lib.objects.pyeq_model()

print("###############################################################################")
print("STARTING MERGING MODELS")
print("###############################################################################")

# SET STARTING TIME
merged_model.start_time = time()

###################################################################
# READS THE MODEL PCK
###################################################################

print("-- Reading model pck files")

try:
    lfile = open( args.ipck , 'r')
    lpck_filename = lfile.read().strip()
except:
    print( red("[PYEQ ERROR] could not read: %s" % (args.ipck )))
    sys.exit()

lpck = []

for pck_filename in lpck_filename:
    ###################################################################
    # LOAD MODEL PCK
    ###################################################################

    try:
        print("-- Loading %s (%.2f Gb) " % (pck_filename, os.path.getsize(pck_filename) / 1024 / 1024 / 1024))
        with open(pck_filename, "rb") as f:
            lpck.append( pickle.load(f) )
        f.close()
    except:
        print(red("[PYEQ ERROR] Could not load: %s " % (pck_filename)))
        sys.exit()

###################################################################
# CHECK PCK ARE CONSISTENT
###################################################################

print("-- Checking pck consistency")

ref_model = lpck[0]

for model in lpck[1:]:
    # test model geometry
    if not np.allclose( ref_model , model ):
        print(red("[PYEQ ERROR] pck.geometry are different: %s vs %s " % ( ref_model.name , model.name )))
    # test input_npz
    if ref_model.input_npz != ref_model.model:
        print(red("[PYEQ ERROR] pck.input_npz are different: %s vs %s " % ( ref_model.name , model.name )))

print("-- OK.")

###################################################################
# CHECK PCK DATES ARE CONSISTENT AND DEFINE THE COMBINED MODEL DATES
###################################################################




###################################################################
# COMBINE SLIP
###################################################################

###################################################################
# PRINT RESULTS
###################################################################


