#!/usr/local/geodesy/anaconda38/bin/python
'''
###################################################################
# SCRIPT    : pyeq_insert_model_into_model.py
# AUTHOR    : Jean-Mathieu Nocquet
# DATE      : November 2020
# INPUT     :
# OUTPUT    :
# NOTE      :
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
import copy
import pickle
from colors import red

# PYACS

# PYEQ
import pyeq.lib.objects
import pyeq.optimization.wrapper.make_inversion

###################################################################
# PARSE COMMAND LINE
###################################################################

prog_info="Inserts a model into a reference model"
prog_epilog="J.-M. Nocquet (Geoazur-UCA-IRD-CNRS-OCA) - June 2020"

parser = argparse.ArgumentParser(description=prog_info,epilog=prog_epilog,prefix_chars='-')
parser.add_argument('-odir', action='store', type=str, dest='odir',required=True,help='output directory')
parser.add_argument('-ref_pck', action='store', type=str, dest='ref_pck',required=True,help='reference model as pck')
parser.add_argument('-i_pck', action='store', type=str, dest='i_pck',required=True,help='model to be inserted')
parser.add_argument('-verbose', '-v', action='count',default=0,help='verbose mode')

args = parser.parse_args()

# verbose
if args.verbose>0:
    print("-- Verbose mode")
    verbose=True
else:
    verbose=False


###################################################################
# INIT MODEL
# MODEL WILL INCLUDE ALL INFORMATION
###################################################################

mmodel = pyeq.lib.objects.pyeq_model()

print("###############################################################################")
print("READING IMPUT MODELS")
print("###############################################################################")

# SET STARTING TIME
mmodel.start_time = time()


###################################################################
# READS THE MODEL PCK
###################################################################

print("-- Reading reference model pck files")

try:
    print("-- Loading %s (%.2f Gb) " % ( args.ref_pck, os.path.getsize( args.ref_pck ) / 1024 / 1024 / 1024))
    with open( args.ref_pck, "rb") as f:
        rmodel = pickle.load( f )
    f.close()
except:
    print( red("[PYEQ ERROR] could not read: %s" % (args.ref_pck )))
    sys.exit()

print("-- Reading model pck file to be inserted")

try:
    print("-- Loading %s (%.2f Gb) " % ( args.i_pck, os.path.getsize( args.i_pck ) / 1024 / 1024 / 1024))
    with open( args.i_pck, "rb") as f:
        imodel = pickle.load( f )
    f.close()
except:
    print( red("[PYEQ ERROR] could not read: %s" % (args.i_pck )))
    sys.exit()

###################################################################
# CHECK MODELS ARE CONSISTENT
###################################################################

# check geometry
# check green

###################################################################
# HANDLING DATES
###################################################################
# Here, we define the list of np_mid_model_date_s (dates in seconds)
# for the merged model. Taking np_mid_model_date_s should enable
# direct interpolation of slip rate value

# get the list of model dates to be updated
mmodel.np_mid_mmodel_date_s = np.array( sorted(list(set( rmodel.np_mid_model_date_s.tolist() + imodel.np_mid_model_date_s.tolist()) )))

###################################################################
# INTERPOLATES MODELS AT NEW DATES AND INSERT IN MERGED MODEL
###################################################################

mmodel.nfaults = rmodel.nfaults
mmodel.slip = np.zeros( ( mmodel.np_mid_mmodel_date_s.shape[0] , mmodel.nfaults ) )

# loop on faults
print("###############################################################################")
print("MERGING IMPUT MODELS")
print("###############################################################################")

R_RATE_SLIP_PER_TIME_STEP = rmodel.slip.reshape(-1, rmodel.nfaults, 1)
I_RATE_SLIP_PER_TIME_STEP = imodel.slip.reshape(-1, imodel.nfaults, 1)

print("--- Making time interpolation and merging")
for idx_fault in np.arange( mmodel.nfaults ):
    if verbose:
        print("  %04d / %04d " % (idx_fault,mmodel.nfaults.shape[0]))


    rslip = np.interp( mmodel.np_mid_mmodel_date_s, rmodel.np_mid_model_date_s, R_RATE_SLIP_PER_TIME_STEP[:,idx_fault,0], left=np.nan, right=np.nan )
    islip = np.interp( mmodel.np_mid_mmodel_date_s, imodel.np_mid_model_date_s, I_RATE_SLIP_PER_TIME_STEP[:,idx_fault,0], left=np.nan, right=np.nan )
    mmodel.slip[:,idx_fault] = np.where(np.isfinite( islip ),islip,rslip)

    print(idx_fault)
    print(rslip)


###################################################################
# PRINT RESULTS
###################################################################

rmodel.slip = mmodel.slip.flatten()
rmodel.name = 'merged'
rmodel.conf = 'merged'
print("###############################################################################")
print("PRINTING RESULTS")
print("###############################################################################")

# Added memory usage
# This command gives the maximum memory used by the process
# value in Gb for Linux, Mb for Mac OS X

rmodel.memory_usage = pyeq.log.get_process_memory_usage()

model = copy.deepcopy( rmodel )

"""
Print kinematic inversion results for the new pyeq >= 0.50.8
"""

###################################################################
# IMPORT
###################################################################

import numpy as np
import pickle

import pyeq.log.make_dir_pyeq_output


###################################################################
# OUTPUT PRINT_RESULTS DIRECTORY
###################################################################
model.odir = model.name
odir = model.odir
print("-- printing inversion results in: %s" % model.odir)
pyeq.log.make_dir_pyeq_output.make_dir_pyeq_output(model.odir)

###################################################################
# PRINT CONF FILE
###################################################################
#print("-- printing conf files: %s" % (model.odir + '/conf'))
#pyeq.lib.log.print_conf(model)

###################################################################
# SAVE MODEL AND OBSERVATION DATES
###################################################################
print("-- printing model and observation dates in: %s" % model.odir)
model = pyeq.log.print_dates(model)

###################################################################
# SAVE PYEQ_KINEMATICS COMMAND LINE IN INFO DIR
###################################################################

print('-- saving pyeq command line in ', model.odir + '/info/command_line.dat')

fcmd = open(model.odir + '/info/command_line.dat', 'w')
fcmd.write(model.cmd_line)
fcmd.close()

###########################################################################
# SAVE GEOMETRY IN INFO DIR
###########################################################################

print("-- saving geometry file in %s/info/geometry.dat" % model.odir)
model = pyeq.log.print_geometry(model)

###########################################################################
# SAVE WARNING FILE
###########################################################################

fwarning = open(model.odir + '/info/warning.dat', 'w')
fwarning.write(model.warning)
fwarning.close()

###########################################################################
# NPY TENSORS
# RATE_SLIP, DELTA_SLIP & SLIP
# HANDLES ORIGIN TIME CONSTANT
# HANDLES VARIABLES RAKE CASE
# SAVE SOME NPY FILES in NPY DIR
###########################################################################

print("-- saving slip, input_npz, t_obs, green tensors as npy in dir: %s" % (model.odir + '/npy'))

# SOLUTION
np.save(model.odir + '/npy/slip.npy', model.slip)

# COPY INPUT NPZ
#copyfile(model.input_npz, model.odir + '/npy/input.npz')

# OBS TENSOR
###########################################################################
np.save(model.odir + '/npy/t_obs.npy', model.t_obs)

# GREEN TENSOR
###########################################################################
np.save(model.odir + '/npy/green.npy', model.green)

# INCREMENTAL SLIP, SLIP RATE
###########################################################################
model = pyeq.log.print_sol_to_slip(model)

# Geometry
###########################################################################
np.save(model.odir + '/npy/geometry.npy', model.geometry)

###########################################################################
# SAVE SLIP TIME SERIES AS TEXT FILE
###########################################################################

print("-- saving text files: rate, incremental and cumulative slip time series in %s/slip_time_series" % (model.odir))
model = pyeq.log.print_slip_time_series(model)

###########################################################################
# STF, INC_STF, CSTF as npy and text file
###########################################################################

print("-- print stf files in %s and %s" % (model.odir + '/npy', model.odir + '/stf'))
model = pyeq.log.print_stf(model)

###########################################################################
# PRINT SLIP MODEL
###########################################################################

print("-- print slip models in %s" % (model.odir + '/slip'))
# no model here
pyeq.log.print_slip_model(model)

###########################################################################
# PRINT MODEL PREDICTED TIME SERIES
###########################################################################

model = pyeq.log.print_modeled_time_series(model)

###########################################################################
# OBS TIME SERIES REALLY USED IN THE INVERSION
###########################################################################

model = pyeq.log.print_observed_time_series(model)

###########################################################################
# RESIDUAL TIME SERIES
###########################################################################
# no return here

model = pyeq.log.print_residual_time_series(model)

###########################################################################
# STATS
###########################################################################

model = pyeq.log.print_stats(model)

###########################################################################
# MODEL/OBS/RESIDUAL DISPLACEMENT FIELDS FOR EACH OBSERVATION DATE
###########################################################################

pyeq.log.print_displacement_fields(model)

###########################################################################
# WRITE SUMMARY
###########################################################################

model.M0 = model.CSTF[-1]
model.magnitude = 2. / 3. * (np.log10(model.M0) - 9.1)

model = pyeq.log.print_summary(model)

print("-- all results written in %s " % model.odir)

###################################################################
# SAVE MODEL PCK (MODEL AS A PICKLE)
###################################################################
print("-- writting model as pickle in %s" % (odir + '/npy/model.pck'))
ofile = open(odir + '/npy/model.pck', 'wb')
pickle.dump(model, ofile, pickle.HIGHEST_PROTOCOL)
ofile.close()





