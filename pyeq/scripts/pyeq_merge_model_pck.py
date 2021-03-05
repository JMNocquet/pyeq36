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
import numpy as np
from time import time
import pickle
from colors import red

# PYACS

# PYEQ
import pyeq.lib.objects
import pyeq.optimization.wrapper.make_inversion

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
# HANDLING DATES
###################################################################
# Here, we define the list of np_mid_model_date_s (dates in seconds)
# for the merged model. Taking np_mid_model_date_s should enable
# direct interpolation of slip rate value


# get the list of model dates to be updated
np_mid_mmodel_date_s = np.array( set( ref_model.np_mid_model_date_s.tolist() + i_model.np_mid_model_date_s.tolist()) )

###################################################################
# INTERPOLATES MODELS AT NEW DATES
###################################################################

slip_ref = np.zeros( ( model.nfaults, np_mid_mmodel_date_s ) )


###################################################################
# COMBINE SLIP
###################################################################

###################################################################
# INSERT SLIP
###################################################################

import copy

mmodel = copy.deepcopy( ref_model )


###################################################################
# PRINT RESULTS
###################################################################


