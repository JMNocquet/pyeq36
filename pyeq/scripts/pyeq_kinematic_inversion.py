#!/usr/bin/env python
'''
###################################################################
# SCRIPT    : pyeq_kinematic_inversion.py
# AUTHOR    : 
# DATE      : March 2013 - February 2018 - March/April 2020
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
# INIT MODEL
# MODEL WILL INCLUDE ALL INFORMATION
###################################################################

model = pyeq.lib.objects.pyeq_model()

print("###############################################################################")
print("STARTING PYEQ")
print("###############################################################################")


# SET STARTING TIME
model.start_time = time()

###################################################################
# PARSE COMMAND LINE
###################################################################

args = pyeq.lib.conf.parse_command_line()

###################################################################
# ARGS PARSING & INIT MODEL
###################################################################
model = pyeq.lib.conf.args_to_model( model , args )

###################################################################
# COLLECTING SYSTEM INFORMATION
###################################################################

model = pyeq.lib.conf.get_resources_info( model )

###################################################################
# NOT IMPLEMENTED YET OR NOT TESTED
###################################################################

if model.rake_constraint >0:
    print('!!!ERROR: variable rake option not implemented yet.')
    sys.exit()

###################################################################
# IS THAT AN INVERSION ONLY RUN ?
###################################################################

if model.pck is not None:
    model = pyeq.lib.conf.run_from_pck( model )

else:

    ###################################################################
    # READING INPUT NPZ FILES
    ###################################################################
    
    print("###############################################################################")
    print("READING GREEN TENSOR")
    print("###############################################################################")
    if model.debug:model_tmp = copy.deepcopy( model )
    
    print("-- Reading input file: %s " % ( model.input_npz ) )
    
    model.sgeometry , model.geometry, model.dm, model.green, _GREEN_UP, model.obs, model.name_obs, _OBS_UP, _NAME_OBS_UP = \
            pyeq.lib.elastic_tensor.read_pyeq_input_npz( model.input_npz )
    if model.verbose:
        pyeq.lib.log.print_model_tensors_shape( model )
    
    if model.debug:pyeq.lib.log.print_diff_model_attributes( model_tmp , model )

    print("-- Shape of green tensor: " , model.green.shape )

    
    ###################################################################
    # REMOVE SOME DISLOCATIONS IF USER PROVIDED ARGUMENT
    ###################################################################
    
    if ( model.geometry_remove_idx is not None ) or ( [model.geometry_range_lon , model.geometry_range_lat , model.geometry_range_depth ] != [None,None,None,None]):
        
        print("-- Excluding user-requested fault elements")
        
        model.green, model.geometry, model.sgeometry , model.dm = \
                pyeq.lib.elastic_tensor.exclude_dislocation(  model.green , model.geometry, model.dm, \
                                                              exclude_idx = model.geometry_remove_idx, \
                                                              range_lon   = model.geometry_range_lon , \
                                                              range_lat   = model.geometry_range_lat , \
                                                              range_depth = model.geometry_range_depth , \
                                                              verbose     = model.verbose)
    
    
    print("-- Shape after pyeq.lib.elastic_tensor.exclude_dislocation: " , model.green.shape )
    #if model.verbose:
    #    print("-- Shape after pyeq.lib.elastic_tensor.exclude_dislocation")
    #    pyeq.lib.log.print_model_tensors_shape( model )
    
    
    ##############################################################################################################################
    # RE-ORDER THE GREEN TENSOR ACCORDING TO THE MAIN AND CONJUGATE RAKE
    ##############################################################################################################################
    if model.debug:model_tmp = copy.deepcopy( model )
    
    
    print("-- Reformatting the GREEN tensor to account for the main and conjugate rake")
    model.green = pyeq.lib.green_tensor.GREEN_ENU_TO_GREEN_RAKE_MAIN_RAKE_CONJUGATE(\
                                        model.green, model.geometry, model.sgeometry, model.rake_type, model.rake_value)
    if model.verbose:
        pyeq.lib.log.print_model_tensors_shape( model )
    
    if model.debug:pyeq.lib.log.print_diff_model_attributes( model_tmp , model )
    
    ###################################################################
    # READ INPUT GPS TIME SERIES
    ###################################################################
    if model.debug:model_tmp = copy.deepcopy( model )
    
    print("###############################################################################")
    print("READING GPS TIME SERIES")
    print("###############################################################################")
    
    # ALL PYACS SGTS FORMAT ARE ALLOWED
    # THE MOST EFFICIENT FORMAT IS PCK (PYTHON PICKELS)
    # OR TSR (OBS_TENSOR) SPECIFIED (STILL NEEDS TO BE IMPLEMENTED)
    ###################################################################
    model = pyeq.lib.gps_time_series.read( model )
    
    if model.debug:pyeq.lib.log.print_diff_model_attributes( model_tmp , model )
    
    ###########################################################################
    # RESCALE UNCERTAINTY
    ###########################################################################
    
    if model.s_h ==0:
        print("-- setting data uncertainty to 1 mm")
        model.t_obs_raw[:,:,3:] = 1.
        model.h_uncertainty = 'set to 1'
    else:
        print("-- rescaling data uncertainty by: %.2lf" % model.s_h )
        model.t_obs_raw[:,:,3:5] = model.t_obs_raw[:,:,3:5] * model.s_h 
        model.h_uncertainty = ("rescaled by %.1f" % ( model.s_h) )
    
    ###########################################################################
    # RESCALE SIGMA UP
    ###########################################################################
    
    print("-- rescaling up sigma by: %.2lf" % model.s_up )
    model.t_obs_raw[:,:,5] = model.t_obs_raw[:,:,5] * model.s_up 
    
    
    ###################################################################
    # DEALING WITH DATES
    ###################################################################
    
    print("###############################################################################")
    print("GET MODEL DATES")
    print("###############################################################################")
    model.np_model_date_s = pyeq.lib.date.get_np_dates_from_arg( model.dates , model.np_obs_date_s, rounding=model.rounding , verbose=model.debug)
    
    ##############################################################################################################################
    # MAKING OBSERVATION AND GREEN TENSORS MUTUALLY CONSISTENT - NEW PYEQ >= 0.50.3
    ##############################################################################################################################
    
    if model.debug:model_tmp = copy.deepcopy( model )
    
    print("###############################################################################")
    print("CHECKING THAT OBSERVATIONS AND GREEN TENSOR ARE MUTUALLY CONSISTENT")
    print("###############################################################################")
    
    model = pyeq.lib.elastic_tensor.check_obs_vs_green( model )
    if model.verbose:
        pyeq.lib.log.print_model_tensors_shape( model )
    
    print("###############################################################################")
    print("CHECKING THAT OBSERVATION AND MODEL DATES ARE MUTUALLY CONSISTENT")
    print("###############################################################################")
    
    model = pyeq.lib.date.check_date_obs_vs_model( model )
    if model.verbose:
        pyeq.lib.log.print_model_tensors_shape( model )
    
    print("###############################################################################")
    print("SET FIRST OBSERVATION AS ZERO")
    print("###############################################################################")
    
    #### data so far is not yet set to zeros for the first obs, correct that
    print("-- Setting the first observation as the reference using pyeq.lib.obs_tensor.set_zero_at_first_obs")
    model.t_obs = pyeq.lib.obs_tensor.set_zero_at_first_obs.set_zero_at_first_obs( model.t_obs , verbose=model.verbose )
    
    if model.debug:pyeq.lib.log.print_diff_model_attributes( model_tmp , model )
    
    ##############################################################################################################################
    # NORMAL SYSTEM
    ##############################################################################################################################
    
    if model.debug:model_tmp = copy.deepcopy( model )
    
    print("###############################################################################")
    print("BUILDING OBSERVATION NORMAL SYSTEM")
    print("###############################################################################")
    
    model.memory_before_ons = pyeq.lib.log.get_process_memory_usage()
    
    # now these are firmly defined
    
    model.nfaults     = model.green.shape[1]
    model.nstep       = model.np_model_date_s.shape[0]-1
    model.nparameters = model.nfaults * model.nstep 
    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Change that !!!!
    model.ncomponent = 2
    
    print('-- model.up ' , model.up )

    print('-- number of model parameters ', model.nparameters ) 

    
    T0 = time()
    
    if model.build == 0:
        print("-- Building normal system using pyeq.lib.forward_model.build0")
        model.normal_system_build  = 'pyeq.lib.build.G_d'
        (G1,d1,diag_Cd1) = pyeq.lib.build.G_d( model , tol_date=5, window=0.1 )
        # still needs to build ATPA / ATPB from G1 & d1
    
    if model.build == 1:
        print("-- Building normal system using pyeq.lib.forward_model.build1")
        model.normal_system_build = 'pyeq.lib.forward_model.build1'
        (ATPA_NEW,ATPB_NEW) = pyeq.lib.forward_model.build1( model , date_tol=0.00000001 )
        
    if model.build == 2:
        print("-- Building normal system using pyeq.lib.build2.build2")
        model.normal_system_build = 'pyeq.lib.build2.build2'
        (ATPA_NEW,ATPB_NEW) = pyeq.lib.forward_model.build2( model )
    
    if model.build == 3:
        print("-- Building normal system using pyeq.lib.build3.build3")
        ( model.N ,  model.Nd ) = pyeq.lib.forward_model.build3( model )
        model.nconstant = 0
    
    if model.build == 4:
        print("-- Building normal system using pyeq.lib.build4.build4")
        RAKE=False
        model.normal_system_build = 'pyeq.lib.build4.build4'
        ( model.N ,  model.Nd ) = pyeq.lib.forward_model.build4( model )
        #N = model.N
        #Nd = model.Nd
    
    if model.build == 5:
        print("-- Building normal system using pyeq.lib.build5.build5")
        RAKE=False
        model.normal_system_build = 'pyeq.lib.build4.build5'
        ( model.N ,  model.Nd ) = pyeq.lib.forward_model.build5( model )
    model.time_build = time() - T0
    print(("-- time building the observation normal system: %.1lf s" % (model.time_build))) 

    if model.debug:pyeq.lib.log.print_diff_model_attributes( model_tmp , model )
    if model.debug:pyeq.lib.log.display_array( model.N *1.E6 )
    
    print('-- model.green  shape ' , model.green.shape )
    print('-- model.N  shape ' , model.N.shape )
    print('-- model.N without constants ' ,  model.nfaults * model.nstep )
    print('-- model.Nd shape ' , model.Nd.shape )
    print('-- model.nfaults ' ,  model.nfaults)
    print('-- model.nstep ' ,  model.nstep )
    print('-- n constant ' ,  model.N.shape[0] - ( model.nstep * model.nfaults ) )
    print('-- n constant vs n gps ' , ( model.N.shape[0] - model.nstep * model.nfaults ) / 2 , model.green.shape[0] )
    
    
    ##############################################################################################################################
    # SAVE OPTION
    ##############################################################################################################################
    
    if model.save:
        print("###############################################################################")
        print("SAVING OBSERVATION NORMAL SYSTEM IN MODEL_NOBS.PCK")
        print("###############################################################################")
        ofile = open( 'model_Nobs.pck', 'wb') 
        pickle.dump( model , ofile , pickle.HIGHEST_PROTOCOL)
        ofile.close()


##############################################################################################################################
# PART FOR REGULAR & RUN FROM MODEL.PCK 
##############################################################################################################################


##############################################################################################################################
# SLIP BOUNDS & CONSTRAINTS
# NOT (RE)-IMPLEMENTED YET
##############################################################################################################################
if model.debug:model_tmp = copy.deepcopy( model )

#model = pyeq.lib.regularization.bounds()

###########################################################################
# GET SIGMA INFO
###########################################################################
print('-- deciphering sigma arg and preparing input constraints')
model = pyeq.lib.regularization.decipher_sigma_arg( model )

print("###############################################################################")
print("REGULARIZATION %s - %s - TEMPORAL CORRELATION %d DAYS" % ( model.regularization.upper() , model.sigma_type.upper() , model.tau ) )
print("###############################################################################")

T0 = time()

model.memory_before_rns = pyeq.lib.log.get_process_memory_usage()

##############################################################################################################################
# LAPLACIAN REGULARIZATION
##############################################################################################################################

if model.regularization == 'laplacian':
    #model = pyeq.lib.regularization.make_regularization_laplacian( model )
    model = pyeq.lib.regularization.make_regularization_valette_01( model )

##############################################################################################################################
# COVARIANCE REGULARIZATION
##############################################################################################################################
if model.regularization == 'covariance':
    model = pyeq.lib.regularization.make_model_covariance_03( model )

model.time_regularization = time() - T0
print(("-- time building regularization normal system: %.1lf s" % (model.time_regularization))) 
print("-- memory usage: %.2lf Gb" % pyeq.lib.log.get_process_memory_usage() )

if model.debug:pyeq.lib.log.print_diff_model_attributes( model_tmp , model )


# print("###############################################################################")
# print("MERGING OBSERVATION AND REGULARIZATION NORMAL SYSTEMS")
# print("###############################################################################")
model.memory_before_mns = pyeq.lib.log.get_process_memory_usage()
model.time_merging_obs_regularization = time() - T0

# no prior value for model
print("-- Adding observation and normalization RHS vector")
#ATB = ATPB_NEW
print("-- memory usage: %.2lf Gb" % pyeq.lib.log.get_process_memory_usage() )

if model.debug:pyeq.lib.log.print_diff_model_attributes( model_tmp , model )

model.N_size = model.N.nbytes / 1024 / 1024 / 1024


###################################################################
# MODIFIES THE NORMAL SYSTEM SO THAT NEGATIVE VALUES FOR SOME
# PARAMETERS ARE ALLOWED, OR CHANGE THE LOWER BOUND.
# HANDLES ORIGIN TIME OFFSET CASE
###################################################################

if model.offset > 0:
    model.shift_constant = -10 # 10 mm should be enough
    print("-- modifying the normal system to handle offset at the origin time")
    model.Nd +=  np.dot ( model.N[ : , -model.nconstant: ] , np.ones( ( model.nconstant , 1 ) ) * -model.shift_constant ).flatten() 
else:
    model.shift_constant = 0.

###################################################################
# MODIFIES THE NORMAL SYSTEM SO THAT NEGATIVE VALUES FOR SOME
# PARAMETERS, OR CHANGE THE LOWER BOUND.
# HANDLES THE CASE OF NEGATIVE (BACK-SLIP BOUNDED) SLIP  
# FOR INTERSEISMIC MODELLING
###################################################################

if model.interseismic != 0:
    print("###############################################################################")
    print("INTERSEISMIC CASE")
    print("###############################################################################")

    print("-- modifying the normal system to allow negative bounded slip for the interseismic case")
    if model.nconstant == 0:
        model.Nd +=  np.dot ( model.N , np.ones( ( ( model.nfaults * model.nstep ) , 1 ) ) * -model.interseismic/365.25 ).flatten() 
    else:
        model.Nd +=  np.dot ( model.N[ : , :-model.nconstant ] , np.ones( ( ( model.nfaults * model.nstep ) , 1 ) ) * -model.interseismic/365.25 ).flatten() 
        
else:
    print("###############################################################################")
    print("TRANSIENT SLIP CASE")
    print("###############################################################################")

###################################################################
# INVERSION
###################################################################

if model.no_opt:
    print("-- no_opt option. User requested end.")
    sys.exit()

print("###############################################################################")
print("MAKING INVERSION")
print("###############################################################################")

model.memory_before_inv = pyeq.lib.log.get_process_memory_usage()

print("-- Number of parameters: %d" % ( model.nparameters ) )
T0 = time()
model.parameters , time_inversion = pyeq.lib.make_inversion.pyeq_nnls( model.N, model.Nd , model.nnls, verbose=model.verbose) 
model.time_inversion = time() - T0

if model.interseismic != 0:
    model.parameters[ : model.nfaults*model.nstep ] = model.parameters[ : model.nfaults*model.nstep ] + model.interseismic/365.25

###################################################################
# REMOVE N & Nd ATTRIBUTES FROM MODEL
###################################################################

print("-- Deleting Normal System from model")

try:
    delattr( model,  'N' )
    delattr( model,  'Nd' )
except:
    print("!!!ERROR deleting normal system attributes")

###################################################################
# NAME OF THE RUN - THINK ABOUT BETTER NAMING
###################################################################

model.name = pyeq.lib.conf.make_model_name( model )
model.odir = model.name

print("-- run name: %s" % model.name )

###################################################################
# PRINT RESULTS
###################################################################
print("###############################################################################")
print("PRINTING RESULTS")
print("###############################################################################")

# Added memory usage
# This command gives the maximum memory used by the process
# value in Gb for Linux, Mb for Mac OS X

model.memory_usage = pyeq.lib.log.get_process_memory_usage()
pyeq.lib.log.print_results( model )

###################################################################
# MAKE PLOTS
###################################################################
if model.plot:
    print("###############################################################################")
    print("MAKING PLOTS")
    print("###############################################################################")
    
    import pyeq.lib.plot

    pyeq.lib.plot.make_plot( model )

###################################################################
# MAKE A COMPRESSED TAR FILE
###################################################################

if model.tar:

    print("###############################################################################")
    print("MAKING A TAR ARCHIVE")
    print("###############################################################################")
    
    import traceback
    import subprocess
    from datetime import datetime
    
    print("-- creating " , model.odir+'.tar.bz2' )
    
    try:
        str_time = s1 = datetime.now().strftime("-%Y-%m-%d-%H-%M-%S")
        cmd = ['tar', 'cfj', model.odir+str_time+'.tar.bz2', model.odir ]
        output = subprocess.check_output(cmd).decode("utf-8").strip() 
        print(output)
        print("-- removing directory: " , model.odir )
        cmd = ['rm', '-Rf', model.odir ]
        output = subprocess.check_output(cmd).decode("utf-8").strip() 
        print(output)   
               
    except Exception:       
        print(f"E: {traceback.format_exc()}")

