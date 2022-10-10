#!/usr/local/geodesy/anaconda38/bin/python
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
#TODO prior value of slip rate
#TODO save regularization constraints?
#TODO put vertical component working again
#TODO offer more option to set the reference position rather than the simplest 0 at t=0
#TODO save sigma if provided as a file
#TODO change make_model_name

###################################################################
# MODULES IMPORT
###################################################################

# GENERAL

import sys
import numpy as np
from time import time
import copy
import logging

# PYACS

import pyeq.date
import pyeq.regularization
import pyeq.regularization.damping
import pyacs.lib.astrotime as at

# PYEQ
import pyeq.forward_model
import pyeq.lib.objects
import pyeq.log
import pyeq.conf
import pyeq.elastic_tensor
import pyeq.lib.green_tensor
import pyeq.gps_time_series
import pyeq.obs_tensor.set_zero_at_first_obs
import pyeq.optimization.wrapper.make_inversion

import pyeq.message.message as MESSAGE
import pyeq.message.verbose_message as VERBOSE
import pyeq.message.error as ERROR
import pyeq.message.warning as WARNING
import pyeq.message.debug_message as DEBUG

###################################################################
# INIT MODEL
# MODEL OBJECT WILL INCLUDE ALL INFORMATION
###################################################################

logging.getLogger("my_logger").setLevel(logging.INFO)


model = pyeq.lib.objects.pyeq_model()

MESSAGE("STARTING PYEQ", level=1)

# SET STARTING TIME
model.start_time = time()

###################################################################
# PARSE COMMAND LINE
###################################################################

args = pyeq.conf.parse_command_line()

###################################################################
# ARGS PARSING & INIT MODEL FROM CONF FILE AND COMMAND LINE
###################################################################
model = pyeq.conf.args_to_model(model, args)

###################################################################
# SET VERBOSE MODE THROUGH LOGGING MODULE
###################################################################

logging.getLogger("my_logger").setLevel(logging.WARNING)

if not hasattr(model, 'verbose'):
    # WARNING LEVEL: ONLY MAJOR STEPS AND WARNINGS WILL BE PRINTED
    logging.getLogger("my_logger").setLevel(logging.WARNING)
else:
    if model.verbose:
        # VERBOSE MODE
        logging.getLogger("my_logger").setLevel(logging.INFO)
    else:
        # WARNING MODE
        logging.getLogger("my_logger").setLevel(logging.WARNING)

if model.debug:
    # DEBUG MODE
    logging.getLogger("my_logger").setLevel(logging.DEBUG)

if logging.getLogger("my_logger").getEffectiveLevel() == logging.WARNING:
    MESSAGE("verbose mode: WARNING")

if logging.getLogger("my_logger").getEffectiveLevel() == logging.INFO:
    MESSAGE("verbose mode: VERBOSE")

if logging.getLogger("my_logger").getEffectiveLevel() == logging.DEBUG:
    MESSAGE("verbose mode: DEBUG")

###################################################################
# COLLECTING SYSTEM INFORMATION
###################################################################
VERBOSE("Collecting system information")
model = pyeq.conf.get_resources_info(model)

###################################################################
# NOT IMPLEMENTED YET OR NOT TESTED
###################################################################

if model.rake_constraint >0:
    ERROR("variable rake option not implemented yet",exit=True)
    sys.exit()

###################################################################
# IS THAT AN INVERSION ONLY RUN ?
###################################################################

if model.mpck is not None:
    VERBOSE("run from mpck option")
    model = pyeq.conf.run_from_pck(model)

else:

    ###################################################################
    # READING INPUT NPZ FILES
    ###################################################################
    
    MESSAGE("READING INPUT GEOMETRY AND GREEN TENSOR FROM INPUT_NPZ",level=1)

    if model.debug:model_tmp = copy.deepcopy( model )
    
    VERBOSE(("Reading input file: %s " % ( model.input_npz ) ))
    
    model.sgeometry , model.geometry, model.dm, model.green, _GREEN_UP, model.obs, model.name_obs, _OBS_UP, _NAME_OBS_UP = \
            pyeq.elastic_tensor.read_pyeq_input_npz(model.input_npz)

    VERBOSE(("green    tensor shape: %d,%d,%d,%d" % model.green.shape))
    VERBOSE(("obs      tensor shape: %d,%d" % model.obs.shape))
    VERBOSE(("dm       tensor shape: %d,%d" % model.dm.shape))
    VERBOSE(("name_obs tensor shape: %d" % model.name_obs.shape))


    if model.debug: pyeq.log.print_diff_model_attributes(model_tmp, model)


    ###################################################################
    # REMOVE SOME DISLOCATIONS IF USER PROVIDED ARGUMENT
    ###################################################################
    
    if ( model.geometry_remove_idx is not None ) or ( [model.geometry_range_lon , model.geometry_range_lat , model.geometry_range_depth ] != [None,None,None,None]):
        
        MESSAGE("Excluding user-requested fault elements")
        
        model.green, model.geometry, model.sgeometry , model.dm = \
                pyeq.elastic_tensor.exclude_dislocation(model.green, model.geometry, model.dm, \
                                                        exclude_idx = model.geometry_remove_idx, \
                                                        range_lon   = model.geometry_range_lon, \
                                                        range_lat   = model.geometry_range_lat, \
                                                        range_depth = model.geometry_range_depth, \
                                                        verbose     = model.verbose)
    
        MESSAGE(("new green tensor shape: %d,%d,%d,%d" % model.green.shape))


    ##############################################################################################################################
    # RE-ORDER THE GREEN TENSOR ACCORDING TO THE MAIN AND CONJUGATE RAKE
    ##############################################################################################################################
    if model.debug:model_tmp = copy.deepcopy( model )

    MESSAGE("Reformatting GREEN tensor to account for the main and conjugate rake")
    model.green = pyeq.lib.green_tensor.GREEN_ENU_TO_GREEN_RAKE_MAIN_RAKE_CONJUGATE(\
                                        model.green, model.geometry, model.sgeometry, model.rake_type, model.rake_value)

    if model.debug: pyeq.log.print_diff_model_attributes(model_tmp, model)
    
    ###################################################################
    # READ INPUT GPS TIME SERIES
    ###################################################################
    if model.debug:model_tmp = copy.deepcopy( model )
    
    MESSAGE("READING GPS TIME SERIES",level=1)

    # ALL PYACS SGTS FORMAT ARE ALLOWED
    # THE MOST EFFICIENT FORMAT IS PCK (PYTHON PICKELS)
    # OR TSR (OBS_TENSOR) SPECIFIED (STILL NEEDS TO BE IMPLEMENTED)
    ###################################################################
    model = pyeq.gps_time_series.read(model)
    
    if model.debug: pyeq.log.print_diff_model_attributes(model_tmp, model)
    
    ###########################################################################
    # RESCALE UNCERTAINTY
    ###########################################################################

    MESSAGE("DATA UNCERTAINTY",level=1)

    if model.s_h ==0:
        MESSAGE("setting data uncertainty to 1 mm")
        model.t_obs_raw[:,:,3:] = 1.
        model.h_uncertainty = 'set to 1'
    else:
        MESSAGE("rescaling data uncertainty by: %.2lf" % model.s_h )
        model.t_obs_raw[:,:,3:5] = model.t_obs_raw[:,:,3:5] * model.s_h 
        model.h_uncertainty = ("rescaled by %.1f" % ( model.s_h) )
    
    ###########################################################################
    # RESCALE SIGMA UP
    ###########################################################################
    
    MESSAGE("rescaling up sigma by: %.2lf" % model.s_up )
    model.t_obs_raw[:,:,5] = model.t_obs_raw[:,:,5] * model.s_up 
    
    #TODO unweight data

    ###################################################################
    # DEALING WITH DATES
    ###################################################################
    
    MESSAGE("DEFINE MODEL DATES",level=1)

    model.np_model_date_s = pyeq.date.get_np_dates_from_arg(model.dates, model.np_obs_date_s, rounding=model.rounding, verbose=model.debug)
    str_sd = at.seconds2datetime(model.np_model_date_s[0]).isoformat(' ')
    str_ed = at.seconds2datetime(model.np_model_date_s[-1]).isoformat(' ')

    MESSAGE("User requested %d model dates from %s to %s" % (model.np_model_date_s.shape[0], str_sd,str_ed))

    ##############################################################################################################################
    # MAKING OBSERVATION AND GREEN TENSORS MUTUALLY CONSISTENT - NEW PYEQ >= 0.50.3
    ##############################################################################################################################
    
    if model.debug:model_tmp = copy.deepcopy( model )
    
    DEBUG("Checking again that observation tensor and green tensor are mutually consistent")

    model = pyeq.elastic_tensor.check_obs_vs_green(model)

    DEBUG(("green    tensor shape: %d,%d,%d,%d" % model.green.shape))
    DEBUG(("obs      tensor shape: %d,%d" % model.obs.shape))
    DEBUG(("dm       tensor shape: %d,%d" % model.dm.shape))
    DEBUG(("name_obs tensor shape: %d" % model.name_obs.shape))


    VERBOSE("Making dates from observation and model consistent")

    model = pyeq.date.check_date_obs_vs_model(model)

    DEBUG(("green    tensor shape: %d,%d,%d,%d" % model.green.shape))
    DEBUG(("obs      tensor shape: %d,%d" % model.obs.shape))
    DEBUG(("dm       tensor shape: %d,%d" % model.dm.shape))
    DEBUG(("name_obs tensor shape: %d" % model.name_obs.shape))


    str_sd = at.seconds2datetime(model.np_model_date_s[0]).isoformat(' ')
    str_ed = at.seconds2datetime(model.np_model_date_s[-1]).isoformat(' ')
    MESSAGE("Model       will include %04d dates from %s to %s" % (model.np_model_date_s.shape[0], str_sd,str_ed))

    str_sd = at.seconds2datetime(model.np_obs_date_s[0]).isoformat(' ')
    str_ed = at.seconds2datetime(model.np_obs_date_s[-1]).isoformat(' ')
    MESSAGE("Observation will include %04d dates from %s to %s" % (model.np_model_date_s.shape[0], str_sd,str_ed))

    #### data so far is not yet set to zeros for the first obs, correct that
    # TODO IMPROVE THIS BY ADDING OPTIONS
    WARNING("SET FIRST OBSERVATION AS ZERO")
    model.warning += "SET FIRST OBSERVATION AS ZERO\n"
    MESSAGE("Setting the first observation as the reference")
    model.t_obs = pyeq.obs_tensor.set_zero_at_first_obs.set_zero_at_first_obs(model.t_obs, verbose=model.verbose)
    
    if model.debug: pyeq.log.print_diff_model_attributes(model_tmp, model)
    
    ##############################################################################################################################
    # OBSERVATION NORMAL SYSTEM
    ##############################################################################################################################
    
    if model.debug:model_tmp = copy.deepcopy( model )
    
    MESSAGE("BUILDING OBSERVATION NORMAL SYSTEM",level=1)

    model.memory_before_ons = pyeq.log.get_process_memory_usage()
    
    # now these are firmly defined
    
    model.nfaults     = model.green.shape[1]
    model.nstep       = model.np_model_date_s.shape[0]-1
    model.nparameters = model.nfaults * model.nstep 

    if model.up:
        model.ncomponent = 3
    else:
        model.ncomponent = 2


    MESSAGE("Use up component: %s" % model.up )
    MESSAGE("number of faults: %d" % model.nfaults )
    MESSAGE("number of model time steps: %d" % model.nstep )
    MESSAGE("number of model parameters: %d" % model.nparameters )

    
    T0 = time()
    MESSAGE("Building observation normal system. Can be long...")
    if model.build == 0:
        WARNING("This is obsolete routine. Building normal system using pyeq.forward_model.build0")
        model.normal_system_build  = 'pyeq.forward_model.build.G_d'
        (G1,d1,diag_Cd1) = pyeq.lib.build.G_d( model , tol_date=5, window=0.1 )
        # still needs to build ATPA / ATPB from G1 & d1
    
    if model.build == 1:
        WARNING("This is obsolete routine. Building normal system using pyeq.forward_model.build1")
        model.normal_system_build = 'pyeq.forward_model.build1'
        (ATPA_NEW,ATPB_NEW) = pyeq.forward_model.build1(model, date_tol=0.00000001)
        
    if model.build == 2:
        WARNING("This is obsolete routine. Building normal system using pyeq.forward_model.build2")
        model.normal_system_build = 'pyeq.forward_model.build2'
        (ATPA_NEW,ATPB_NEW) = pyeq.forward_model.build2(model)
    
    if model.build == 3:
        WARNING("This is obsolete routine. Building normal system using pyeq.forward_model.build3")
        ( model.N ,  model.Nd ) = pyeq.forward_model.build3(model)
        model.nconstant = 0
    
    if model.build == 4:
        WARNING("This is obsolete routine. Building normal system using pyeq.lib.build4.build4")
        RAKE=False
        model.normal_system_build = 'pyeq.forward_model.build4'
        ( model.N ,  model.Nd ) = pyeq.forward_model.build4(model)
        #N = model.N
        #Nd = model.Nd
    
    if model.build == 5:
        DEBUG("Building normal system using pyeq.forward_model.build5")
        RAKE=False
        model.normal_system_build = 'pyeq.forward_model.build5'
        ( model.N ,  model.Nd ) = pyeq.forward_model.build5(model)

    model.time_build = time() - T0
    MESSAGE(("time building the observation normal system: %.1lf s" % (model.time_build)))
    MESSAGE(("normal system shape: %s" % (model.N.shape,)))
    MESSAGE(("current used memory size (Gb): %.2lf" % (pyeq.log.get_process_memory_usage())))


    if model.debug: pyeq.log.print_diff_model_attributes(model_tmp, model)
    if model.debug: pyeq.log.display_array(model.N * 1.E6)
    
    DEBUG("model.green  shape: %s" % (model.green.shape,) )
    DEBUG("model.N  shape: %s " % (model.N.shape,) )
    DEBUG("model.N without constants: %d " %  model.nfaults * model.nstep )
    DEBUG("model.Nd shape: %s " % (model.Nd.shape,) )
    DEBUG("model.nfaults: %d " %   model.nfaults)
    DEBUG("model.nstep: %d " %  model.nstep )
    DEBUG("n constant: %d " %  (model.N.shape[0] - ( model.nstep * model.nfaults )) )
    DEBUG("n constant vs n gps: %d vs %d " %  (( model.N.shape[0] - model.nstep * model.nfaults ) / 2 , model.green.shape[0]) )
    
    
    ##############################################################################################################################
    # SAVE OPTION
    ##############################################################################################################################
    
    if model.save:
        MESSAGE("SAVING OBSERVATION NORMAL SYSTEM AS PICKLE (.MPCK)" , level=1)
        import pyeq.log

        pyeq.log.save_model_N(model)


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
# START REGULARIZATION PART
###########################################################################

MESSAGE(("BUILDING REGULARIZATION NORMAL SYSTEM: %s  " % ( model.regularization.upper() ) ),level=1)

###########################################################################
# GET SIGMA INFO
###########################################################################
#if model.regularization == 'covariance':
VERBOSE("deciphering sigma argument")
model = pyeq.regularization.damping.decipher_sigma_arg( model )

T0 = time()

model.memory_before_rns = pyeq.log.get_process_memory_usage()

##############################################################################################################################
# LAPLACIAN REGULARIZATION
##############################################################################################################################

if model.regularization == 'laplacian':
    import pyeq.regularization.laplace
    model = pyeq.regularization.laplace.add_laplace_cons( model )
    model = pyeq.regularization.damping.add_damping_cons( model )


##############################################################################################################################
# COVARIANCE REGULARIZATION
##############################################################################################################################

# COVARIANCE
if model.regularization == 'covariance':
    model = pyeq.lib.regularization.make_model_covariance_03( model )

##############################################################################################################################
# OBSOLETE OR NOT FINISHED IMPLEMENTATIONS
##############################################################################################################################

# This is the old Laplacian
#if model.regularization == 'laplacian':
#    #model = pyeq.lib.regularization.make_regularization_laplacian( model )
#    model = pyeq.lib.regularization.make_regularization_valette_01( model )
# LAPLACIAN_LIKE
#if model.regularization == 'laplacian_like':
#    model = pyeq.lib.regularization.make_regularization_laplacian_like( model )
# CVXOPT
#if model.regularization == 'cvxopt':
#    model = pyeq.lib.regularization.get_regularization_operators( model )


model.time_regularization = time() - T0
VERBOSE(("time building regularization normal system: %.1lf s" % (model.time_regularization)))
VERBOSE("memory usage: %.2lf Gb" % pyeq.log.get_process_memory_usage())

if model.debug: pyeq.log.print_diff_model_attributes(model_tmp, model)

model.memory_before_mns = pyeq.log.get_process_memory_usage()
model.time_merging_obs_regularization = time() - T0

# no prior value for model
#print("-- Adding observation and normalization RHS vector")
#ATB = ATPB_NEW
VERBOSE("memory usage: %.2lf Gb" % pyeq.log.get_process_memory_usage())

if model.debug: pyeq.log.print_diff_model_attributes(model_tmp, model)

model.N_size = model.N.nbytes / 1024 / 1024 / 1024
VERBOSE("Normal matrix memory size: %.2lf Gb" % model.N_size)


###################################################################
# MODIFIES THE NORMAL SYSTEM SO THAT NEGATIVE VALUES FOR SOME
# PARAMETERS ARE ALLOWED, OR CHANGE THE LOWER BOUND.
# HANDLES ORIGIN TIME OFFSET CASE
###################################################################

if model.offset > 0:
    model.shift_constant = -10 # 10 mm should be enough
    MESSAGE("modifying the normal system to handle offset at the origin time")
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
    MESSAGE(("INTERSEISMIC CASE"),level=1)
    ERROR("NOT IMPLEMENTED YET",exit=True)
    print("-- modifying the normal system to allow negative bounded slip for the interseismic case")
    if model.nconstant == 0:
        model.Nd +=  np.dot ( model.N , np.ones( ( ( model.nfaults * model.nstep ) , 1 ) ) * -model.interseismic/365.25 ).flatten() 
    else:
        model.Nd +=  np.dot ( model.N[ : , :-model.nconstant ] , np.ones( ( ( model.nfaults * model.nstep ) , 1 ) ) * -model.interseismic/365.25 ).flatten() 
        
else:
    MESSAGE("TRANSIENT SLIP CASE",level=1)

###################################################################
# NO_OPT CASE
###################################################################

if model.no_opt:
    MESSAGE("no_opt option set to True. Stopping before optimization")
    sys.exit()

###################################################################
# RUNNING INVERSION
###################################################################
MESSAGE("MAKING INVERSION",level=1)

model.memory_before_inv = pyeq.log.get_process_memory_usage()
MESSAGE("memory usage before inversion: %.2lf Gb" % model.memory_before_inv)
MESSAGE("Number of estimated parameters: %d" % ( model.nparameters ) )

T0 = time()
if model.regularization in ['covariance','laplacian']:
    model.parameters , time_inversion = pyeq.optimization.wrapper.make_inversion.pyeq_nnls(model.N, model.Nd, model.nnls, verbose=model.verbose)

# NOT IMPLEMENTED YET
if model.regularization in ['cvxopt']:
    model.parameters , time_inversion = pyeq.optimization.wrapper.make_inversion.pyeq_cvxopt(model)

model.time_inversion = time() - T0
VERBOSE(("time inversion: %.1lf s" % (model.time_inversion)))

# NOT IMPLEMENTED YET
if model.interseismic != 0:
    model.parameters[ : model.nfaults*model.nstep ] = model.parameters[ : model.nfaults*model.nstep ] + model.interseismic/365.25


###################################################################
# RECORDS MEMORY USAGE
###################################################################
model.memory_usage = pyeq.log.get_process_memory_usage()
VERBOSE("memory usage after inversion: %.2lf Gb" % model.memory_usage)

###################################################################
# REMOVE N & Nd ATTRIBUTES FROM MODEL
###################################################################

VERBOSE("Deleting Normal System from model")

try:
    delattr( model,  'N' )
    delattr( model,  'Nd' )
except:
    ERROR("deleting normal system attributes")

###################################################################
# NAME OF THE RUN - THINK ABOUT BETTER NAMING
###################################################################

model.name = pyeq.conf.make_model_name(model)
model.odir = model.name
# TODO improve name with 0 instead of E
MESSAGE("run name: %s" % model.name )

###################################################################
# PRINT RESULTS OR JUST SAVE MODEL.PCK
###################################################################

if model.print_result:

    MESSAGE("PRINTING RESULTS",level=1)
    pyeq.log.print_results(model)

else:
    model_pck = 'model_'+ model.name + '.mpck'
    MESSAGE(("SAVING MODEL AS PICKLE (.mpck): %s" % model_pck),level=1)
    model.write_pickle( model_pck )
    sys.exit()

###################################################################
# MAKE PLOTS
###################################################################
if model.plot:
    MESSAGE("MAKING PLOTS",level=1)

    import pyeq.plot
    pyeq.plot.make_plot(model)

###################################################################
# MAKE A COMPRESSED TAR FILE
###################################################################

if model.tar:
    MESSAGE(("MAKING A COMPRESSED TAR ARCHIVE"),level=1)

    import traceback
    import subprocess
    from datetime import datetime
    
    MESSAGE("creating %s.tar.bz2" % model.odir )
    
    try:
        str_time = s1 = datetime.now().strftime("-%Y-%m-%d-%H-%M-%S")
        cmd = ['tar', 'cfj', model.odir+str_time+'.tar.bz2', model.odir ]
        output = subprocess.check_output(cmd).decode("utf-8").strip() 
        DEBUG(output)
        VERBOSE("removing directory: %s" % model.odir )
        cmd = ['rm', '-Rf', model.odir ]
        output = subprocess.check_output(cmd).decode("utf-8").strip() 
        print(output)   
               
    except Exception:       
        print(f"E: {traceback.format_exc()}")

