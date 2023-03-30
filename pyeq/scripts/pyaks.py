#!/usr/bin/env python
'''
###################################################################
# SCRIPT    : pyaks.py
# AUTHOR    :
# DATE      : March 2013 - February 2018 - March/April 2020
# INPUT     :
# OUTPUT    :
# NOTE      : pyaks is the test version of pyeq_kinematic_inversion.py
#           : adding variable rake 26/03/2021
###################################################################
'''
# TODO prior value of slip rate
# TODO save regularization constraints?
# TODO offer more option to set the reference position rather than the simplest 0 at t=0
# TODO save sigma if provided as a file
# TODO change make_model_name
# TODO reactivate printing resolution.dat ; there is a compatibility issue when model.geometry_type = 'TDV'
# TODO implement model.geometry_type = 'RDE' (rectangular dislocation element)
# TODO check incompatible options: mpck option is not compatible with backwards and cross_validation

###################################################################
# MODULES IMPORT
###################################################################

# GENERAL

import sys
import numpy as np
from time import time
import copy
import logging
import pkg_resources

# PYACS

import pyacs
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
# change 27/03/2021
import pyeq.elastic_tensor.green_ds_ss_to_main_rake_conjugate
import pyeq.gps_time_series
import pyeq.obs_tensor.set_zero_at_first_obs
import pyeq.obs_tensor.set_zero_obs
import pyeq.optimization.wrapper.make_inversion

# VERBOSE
import pyeq.message.message as MESSAGE
import pyeq.message.verbose_message as VERBOSE
import pyeq.message.error as ERROR
import pyeq.message.warning as WARNING
import pyeq.message.debug_message as DEBUG

from art import tprint
tprint('PYAKS')
MESSAGE("PYEQ  version: %s" % pkg_resources.get_distribution("pyeq").version)
MESSAGE("PYACS version: %s" % pkg_resources.get_distribution("pyacs").version)


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

MESSAGE("regularization parameters after args: %s %s %s %s %s " % (model.sigma, model.lambda_spatial_smoothing,
                                                        model.lambda_temporal_smoothing,
                                                        model.lambda_final_spatial_smoothing, model.lambda_stf))

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

#if model.rake_constraint > 0:
#    ERROR("variable rake option not implemented yet", exit=True)
#    sys.exit()

###################################################################
# RUN CROSS VALIDATION
###################################################################

if model.cross_validation == 'run':
    # residuals
    def make_res_model(obs, model):
        # takes displacement tensors and computes basic stats
        # corrects model so that 0 is at the first epoch of obs
        import numpy as np
        # find and correct bias for each site i
        mmodel = np.copy(model)
        for i in np.arange(obs.shape[1]):
            #print('site i',i)
            #print(obs[:, i, 0])
            #print('Is Nan ? ' , np.isnan(obs[:, i, 0]).all())
            if np.isnan(obs[:, i, 0]).all():continue
            idx = np.where(~np.isnan(obs[:, i, 0]))[0][0]
            mmodel[:, i, :] = mmodel[:, i, :] - model[idx, i, :]
        # residuals
        res = obs[:, :, :3] - mmodel
        return res

    # rms routine
    def make_rms_model(obs, model):
        # takes displacement tensors and computes basic stats
        # corrects model so that 0 is at the first epoch of obs
        import numpy as np
        # find and correct bias for each site i
        mmodel = np.copy(model)
        for i in np.arange(obs.shape[1]):
            #print('site i',i)
            #print(obs[:, i, 0])
            #print('Is Nan ? ' , np.isnan(obs[:, i, 0]).all())
            if np.isnan(obs[:, i, 0]).all():continue
            idx = np.where(~np.isnan(obs[:, i, 0]))[0][0]
            mmodel[:, i, :] = mmodel[:, i, :] - model[idx, i, :]
        # residuals
        res = obs[:, :, :3] - mmodel
        # computes statistics
        rms_h = np.sqrt(np.nanmean(res[:, :, :2] ** 2))
        rms_v = np.sqrt(np.nanmean(res[:, :, 2] ** 2))
        chi2 = np.nansum((res[:, :, :3] / obs[:, :, 3:]) ** 2)
        num = np.nansum((1. / obs[:, :, 3:]) ** 2)
        wrms_3D = np.sqrt(chi2 / num)

        return rms_h, rms_v, chi2, wrms_3D

    # load reference model

    # import
    import pickle
    import numpy as np
    # load model
    # checking regularization parameters
    MESSAGE("regularization parameters: %s %s %s %s %s " % ( model.sigma, model.lambda_spatial_smoothing,
                model.lambda_temporal_smoothing, model.lambda_final_spatial_smoothing, model.lambda_stf ))
    MESSAGE("Loading reference model model_CROSS_VALIDATION.mpck")
    with open( 'model_CROSS_VALIDATION.mpck', "rb") as f:
        model_ref = pickle.load( f )
    f.close()

    # checking regularization parameters
    MESSAGE("regularization parameters: %s %s %s %s %s " % ( model.sigma, model.lambda_spatial_smoothing,
                model.lambda_temporal_smoothing, model.lambda_final_spatial_smoothing, model.lambda_stf ))

    model_ref.sigma, model_ref.lambda_spatial_smoothing, \
    model_ref.lambda_temporal_smoothing, model_ref.lambda_final_spatial_smoothing, model_ref.lambda_stf = \
    model.sigma, model.lambda_spatial_smoothing, \
    model.lambda_temporal_smoothing, model.lambda_final_spatial_smoothing, model.lambda_stf

    model = model_ref

    MESSAGE("regularization parameters used for cross-validation: %s %s %s %s %s " % ( model.sigma, model.lambda_spatial_smoothing,
                model.lambda_temporal_smoothing, model.lambda_final_spatial_smoothing, model.lambda_stf ))


    # run K0 with all observations
    model.N = np.load('CROSS_VALIDATION_N.npy')
    model.Nd = np.load('CROSS_VALIDATION_Nd.npy')

    import pyeq.regularization.laplace

    model = pyeq.regularization.laplace.add_laplace_cons(model)
    model = pyeq.regularization.damping.add_damping_cons(model)

    if model.regularization in ['covariance', 'laplacian']:
        model.parameters, time_inversion = pyeq.optimization.wrapper.make_inversion.pyeq_nnls(model.N, model.Nd, model.nnls,
                                                                                              verbose=model.verbose)
    # print results

    if model.geometry_type == 'TDV':
        MESSAGE("Reformatting output slip from triangular vertices (TDV) to triangular dislocation elements (TDE)")
        SRN = model.parameters.reshape(-1, model.green_node.shape[1]).T
        model.green = model.green_subfault
        SR = np.zeros((model.matrix_subfault_to_node.shape[0], SRN.shape[1]))

        for i in np.arange(SRN.shape[1]):
            SR[:, i] = np.dot(model.matrix_subfault_to_node, SRN[:, i])
        model.parameters = SR.T.flatten()
        model.nfaults = model.green.shape[1]


    ###################################################################
    # REMOVE N & Nd ATTRIBUTES FROM MODEL
    ###################################################################

    VERBOSE("Deleting Normal System from model")

    try:
        delattr(model, 'N')
        delattr(model, 'Nd')
    except:
        ERROR("deleting normal system attributes")

    ###################################################################
    # INVERSION SUMMARY
    ###################################################################
    MESSAGE(("INVERSION RESULTS SUMMARY ALL DATA"), level=1)
    if model.nconstant > 0:
        model.slip = model.parameters[:-model.nconstant]
        model.estimated_offsets = model.parameters[-model.nconstant:]
    else:
        model.slip = model.parameters
        model.estimated_offsets = None

    import pyeq.log

    model = pyeq.log.print_dates(model, save=False)
    model = pyeq.log.print_sol_to_slip(model, save=False)
    model = pyeq.log.print_stf(model, save=False)
    model.M0 = model.CSTF[-1]
    model = pyeq.log.print_offset(model, save=False)
    model = pyeq.log.print_modeled_time_series(model, save=False)

    if model.M0 > 0:
        model.magnitude = 2. / 3. * (np.log10(model.M0) - 9.1)
    else:
        model.magnitude = -99
    # moment
    MESSAGE("Moment (N.n): %8.2E (Mw%.1lf)" % (model.M0, model.magnitude))
    # max cumulated slip

    # maximum moment rate
    mmr = np.max(model.STF)
    if mmr > 0:
        MESSAGE("Maximum moment rate per day: %8.2E (Mw%.1lf)" % (mmr, 2. / 3. * (np.log10(mmr) - 9.1)))
    else:
        MESSAGE("Maximum moment rate per day: %8.2E (Mw%.1lf)" % (mmr, -99))

    # max slip rate
    i_cs_max = np.unravel_index(np.argmax(model.NORM_CUMULATIVE_SLIP_PER_TIME_STEP),
                                model.NORM_CUMULATIVE_SLIP_PER_TIME_STEP.shape)
    slip_max = model.NORM_CUMULATIVE_SLIP_PER_TIME_STEP[i_cs_max]
    model.idx_cumulative_slip_max = i_cs_max[1]

    MESSAGE("cumulative slip max  (mm)        : %8.1f fault #%d (%.2lf,%.2lf)" % (slip_max, i_cs_max[1],
                                                                                  model.sgeometry.centroid_long[
                                                                                      i_cs_max[1]],
                                                                                  model.sgeometry.centroid_lat[
                                                                                      i_cs_max[1]]))


    rms_h, rms_v, chi2, wrms_3d = make_rms_model(model.t_obs, model.t_mod)

    MESSAGE(("rms_h (mm):%.2lf" % rms_h))
    MESSAGE(("rms_v (mm):%.2lf" % rms_v))
    MESSAGE(("wrms_3d (mm):%.2lf" % wrms_3d))
    MESSAGE(("chi2 (mm**2):%.2lf" % chi2))


    with open('L_CURVE.dat','a') as flcurve:
        flcurve.write("%8.3lf %8.3lf %8.3lf %8.3lf %8.3lf %8.3lf %8.3lf %8.3lf %8.3lf\n" % \
                tuple(map(float,(model.sigma, model.lambda_spatial_smoothing, \
                model.lambda_temporal_smoothing, model.lambda_final_spatial_smoothing, model.lambda_stf, \
                 rms_h, rms_v, chi2, wrms_3d))))


    ###########################################################################
    # Loop on K-fold
    ###########################################################################

    # load index info
    np_idxcv_dates = np.load('CROSS_VALIDATION_idxcv_dates.npy')
    np_idxcv_gps   = np.load('CROSS_VALIDATION_idxcv_gps.npy')

    MSE_CV = 0.

    for k in np.arange(1,10):

        MESSAGE(("Running K-fold: %d with regularization parameters : %s/%s/%s/%s/%s" %
                 (k,model.sigma, model.lambda_spatial_smoothing,
        model.lambda_temporal_smoothing, model.lambda_final_spatial_smoothing, model.lambda_stf)) , level=1)

        MESSAGE("GPS sites with 1/3 data deleted")
        str_message = " ".join(model.np_gps_site[np_idxcv_gps[k,0]:np_idxcv_gps[k,1]])
        MESSAGE(str_message)

        # load N & Nd
        Nfile = ("CROSS_VALIDATION_NK%d.npy" % k )
        Ndfile = ("CROSS_VALIDATION_NdK%d.npy" % k )
        model.N = np.load( Nfile )
        model.Nd = np.load( Ndfile )

        MESSAGE("regularization parameters used for cross-validation: %s %s %s %s %s " % (
        model.sigma, model.lambda_spatial_smoothing,
        model.lambda_temporal_smoothing, model.lambda_final_spatial_smoothing, model.lambda_stf))

        import pyeq.regularization.laplace

        model = pyeq.regularization.laplace.add_laplace_cons(model)
        model = pyeq.regularization.damping.add_damping_cons(model)

        if model.regularization in ['covariance', 'laplacian']:
            model.parameters, time_inversion = pyeq.optimization.wrapper.make_inversion.pyeq_nnls(model.N, model.Nd, model.nnls,
                                                                                                  verbose=model.verbose)
        # print results

        if model.geometry_type == 'TDV':
            MESSAGE("Reformatting output slip from triangular vertices (TDV) to triangular dislocation elements (TDE)")
            SRN = model.parameters.reshape(-1, model.green_node.shape[1]).T
            model.green = model.green_subfault
            SR = np.zeros((model.matrix_subfault_to_node.shape[0], SRN.shape[1]))

            for i in np.arange(SRN.shape[1]):
                SR[:, i] = np.dot(model.matrix_subfault_to_node, SRN[:, i])
            model.parameters = SR.T.flatten()
            model.nfaults = model.green.shape[1]

        ###################################################################
        # REMOVE N & Nd ATTRIBUTES FROM MODEL
        ###################################################################

        VERBOSE("Deleting Normal System from model")

        try:
            delattr(model, 'N')
            delattr(model, 'Nd')
        except:
            ERROR("deleting normal system attributes")

        ###################################################################
        # INVERSION SUMMARY
        ###################################################################
        MESSAGE(("INVERSION RESULTS SUMMARY K=%d" % k), level=1)
        if model.nconstant > 0:
            model.slip = model.parameters[:-model.nconstant]
            model.estimated_offsets = model.parameters[-model.nconstant:]
        else:
            model.slip = model.parameters
            model.estimated_offsets = None

        import pyeq.log

        model = pyeq.log.print_dates(model, save=False)
        model = pyeq.log.print_sol_to_slip(model, save=False)
        model = pyeq.log.print_stf(model, save=False)
        model.M0 = model.CSTF[-1]
        model = pyeq.log.print_offset(model, save=False)
        model = pyeq.log.print_modeled_time_series(model, save=False)

        if model.M0 > 0:
            model.magnitude = 2. / 3. * (np.log10(model.M0) - 9.1)
        else:
            model.magnitude = -99
        # moment
        MESSAGE("Moment (N.n): %8.2E (Mw%.1lf)" % (model.M0, model.magnitude))
        # max cumulated slip

        # maximum moment rate
        mmr = np.max(model.STF)
        if mmr > 0:
            MESSAGE("Maximum moment rate per day: %8.2E (Mw%.1lf)" % (mmr, 2. / 3. * (np.log10(mmr) - 9.1)))
        else:
            MESSAGE("Maximum moment rate per day: %8.2E (Mw%.1lf)" % (mmr, -99))

        # max slip rate
        i_cs_max = np.unravel_index(np.argmax(model.NORM_CUMULATIVE_SLIP_PER_TIME_STEP),
                                    model.NORM_CUMULATIVE_SLIP_PER_TIME_STEP.shape)
        slip_max = model.NORM_CUMULATIVE_SLIP_PER_TIME_STEP[i_cs_max]
        model.idx_cumulative_slip_max = i_cs_max[1]

        MESSAGE("cumulative slip max  (mm)        : %8.1f fault #%d (%.2lf,%.2lf)" % (slip_max, i_cs_max[1],
                                                                                      model.sgeometry.centroid_long[
                                                                                          i_cs_max[1]],
                                                                                      model.sgeometry.centroid_lat[
                                                                                          i_cs_max[1]]))

        rms_h, rms_v, chi2, wrms_3d = make_rms_model(model.t_obs, model.t_mod)

        MESSAGE(("rms_h (mm):%.2lf" % rms_h))
        MESSAGE(("rms_v (mm):%.2lf" % rms_v))
        MESSAGE(("wrms_3d (mm):%.2lf" % wrms_3d))
        MESSAGE(("chi2 (mm**2):%.2lf" % chi2))


        #t_obs = model.t_obs[np_idxcv_dates[k,0]:np_idxcv_dates[k,1],np_idxcv_gps[k,0]:np_idxcv_gps[k,1],:]
        #t_mod = model.t_mod[np_idxcv_dates[k,0]:np_idxcv_dates[k,1],np_idxcv_gps[k,0]:np_idxcv_gps[k,1],:]
        #rms_h, rms_v, chi2, wrms_3d = make_rms_model(t_obs, t_mod)

        res = make_res_model(model.t_obs, model.t_mod)

        # computes statistics
        rms_h = np.sqrt(np.nanmean(res[np_idxcv_dates[k,0]:np_idxcv_dates[k,1],np_idxcv_gps[k,0]:np_idxcv_gps[k,1], :2] ** 2))
        rms_v = np.sqrt(np.nanmean(res[np_idxcv_dates[k,0]:np_idxcv_dates[k,1],np_idxcv_gps[k,0]:np_idxcv_gps[k,1], 2] ** 2))
        chi2 = np.nansum((res[np_idxcv_dates[k,0]:np_idxcv_dates[k,1],np_idxcv_gps[k,0]:np_idxcv_gps[k,1], :3] / model.t_obs[np_idxcv_dates[k,0]:np_idxcv_dates[k,1],np_idxcv_gps[k,0]:np_idxcv_gps[k,1], 3:]) ** 2)
        num = np.nansum((1. / model.t_obs[np_idxcv_dates[k,0]:np_idxcv_dates[k,1],np_idxcv_gps[k,0]:np_idxcv_gps[k,1], 3:]) ** 2)
        wrms_3D = np.sqrt(chi2 / num)


        MESSAGE(("MSE K=%d rms_h (mm):%.2lf" % (k,rms_h)))
        MESSAGE(("MSE K=%d rms_v (mm):%.2lf" % (k,rms_v)))
        MESSAGE(("MSE K=%d wrms_3d (mm):%.2lf" % (k,wrms_3d)))
        MESSAGE(("MSE K=%d chi2 (mm**2):%.2lf" % (k,chi2)))

        with open('MSE_K.dat','a') as flcurve:
            flcurve.write("%8.3lf %8.3lf %8.3lf %8.3lf %8.3lf %d %8.3lf %8.3lf %8.3lf %8.3lf\n" % \
                    tuple(map(float,(model.sigma, model.lambda_spatial_smoothing, \
                    model.lambda_temporal_smoothing, model.lambda_final_spatial_smoothing, model.lambda_stf, \
                     k, rms_h, rms_v, chi2, wrms_3d))))

        MSE_CV = MSE_CV + wrms_3d**2

    with open('MSE_CV.dat','a') as flcurve:
        flcurve.write("%8.3lf %8.3lf %8.3lf %8.3lf %8.3lf %8.3lf\n" % \
                tuple(map(float,(model.sigma, model.lambda_spatial_smoothing, \
                model.lambda_temporal_smoothing, model.lambda_final_spatial_smoothing, model.lambda_stf, \
                 MSE_CV))))


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

    MESSAGE("READING INPUT GEOMETRY AND GREEN TENSOR FROM INPUT_NPZ", level=1)

    if model.debug: model_tmp = copy.deepcopy(model)

    VERBOSE(("Reading input file: %s " % (model.input_npz)))

    model.sgeometry, model.geometry, model.dm, model.green, _GREEN_UP, model.obs, model.name_obs, _OBS_UP, _NAME_OBS_UP = \
        pyeq.elastic_tensor.read_pyeq_input_npz(model.input_npz)

    VERBOSE(("green    tensor shape: %d,%d,%d,%d" % model.green.shape))
    VERBOSE(("obs      tensor shape: %d,%d" % model.obs.shape))
    VERBOSE(("dm       tensor shape: %d,%d" % model.dm.shape))
    VERBOSE(("name_obs tensor shape: %d" % model.name_obs.shape))

    if model.debug: pyeq.log.print_diff_model_attributes(model_tmp, model)


    ###################################################################
    # REMOVE SOME DISLOCATIONS IF USER PROVIDED ARGUMENT
    ###################################################################

    if (model.geometry_remove_idx is not None) or (
            [model.geometry_range_lon, model.geometry_range_lat, model.geometry_range_depth] != [None, None, None,
                                                                                                 None]):
        MESSAGE("Excluding user-requested fault elements")

        model.green, model.geometry, model.sgeometry, model.dm = \
            pyeq.elastic_tensor.exclude_dislocation(model.green, model.geometry, model.dm, \
                                                    exclude_idx=model.geometry_remove_idx, \
                                                    range_lon=model.geometry_range_lon, \
                                                    range_lat=model.geometry_range_lat, \
                                                    range_depth=model.geometry_range_depth, \
                                                    verbose=model.verbose)

        MESSAGE(("new green tensor shape: %d,%d,%d,%d" % model.green.shape))

    ###################################################################
    # SWAP GREEN TENSOR AXIS
    # 28/03/2021
    # for compatibility because
    # original pyeq.lib.green_tensor.GREEN_ENU_TO_GREEN_RAKE_MAIN_RAKE_CONJUGATE
    # swap dislocation (axis=0) and OBS (axis=1)
    ###################################################################

    MESSAGE("Swaping dislocation and observation in Green tensor ordering")
    model.green = np.swapaxes(model.green , 0, 1)

    ###################################################################
    # READ INPUT GPS TIME SERIES
    ###################################################################
    if model.debug: model_tmp = copy.deepcopy(model)

    MESSAGE("READING GPS TIME SERIES", level=1)

    # ALL PYACS SGTS FORMAT ARE ALLOWED
    # THE MOST EFFICIENT FORMAT IS PCK (PYTHON PICKLES)
    # OR TSR (OBS_TENSOR) SPECIFIED (STILL NEEDS TO BE IMPLEMENTED)
    ###################################################################
    model = pyeq.gps_time_series.read(model)

    if model.debug: pyeq.log.print_diff_model_attributes(model_tmp, model)

    ###########################################################################
    # RESCALE UNCERTAINTY
    ###########################################################################

    MESSAGE("DATA UNCERTAINTY", level=1)

    if model.s_h == 0:
        MESSAGE("setting data uncertainty to 1 mm")
        model.t_obs_raw[:, :, 3:] = 1.
        model.h_uncertainty = 'set to 1'
    else:
        MESSAGE("rescaling data uncertainty by: %.2lf" % model.s_h)
        model.t_obs_raw[:, :, 3:5] = model.t_obs_raw[:, :, 3:5] * model.s_h
        model.h_uncertainty = ("rescaled by %.1f" % (model.s_h))

    ###########################################################################
    # RESCALE SIGMA UP
    ###########################################################################

    MESSAGE("rescaling up sigma by: %.2lf" % model.s_up)
    model.t_obs_raw[:, :, 5] = model.t_obs_raw[:, :, 5] * model.s_up

    ###########################################################################
    # RE-WEIGHT
    ###########################################################################

    unweight_factor = 1000.

    if model.lunweight_gps != []:
        for i in np.arange( model.np_names_t_obs.shape[0] ):
             if model.np_names_t_obs[i] in model.lunweight_gps:
                VERBOSE("Re-weighting uncertainty for %s by a factor of 1000. " % model.np_names_t_obs[i])
                model.t_obs_raw[:, i, 3:] = model.t_obs_raw[:, i, 3:] * unweight_factor


    ###################################################################
    # DEALING WITH DATES
    ###################################################################

    MESSAGE("DEFINE MODEL DATES", level=1)

    model.np_model_date_s = pyeq.date.get_np_dates_from_arg(model.dates, model.np_obs_date_s, rounding=model.rounding,
                                                            verbose=model.debug)
    str_sd = at.seconds2datetime(model.np_model_date_s[0]).isoformat(' ')
    str_ed = at.seconds2datetime(model.np_model_date_s[-1]).isoformat(' ')

    MESSAGE("User requested %d model dates from %s to %s" % (model.np_model_date_s.shape[0], str_sd, str_ed))

    ##############################################################################################################################
    # MAKING OBSERVATION AND GREEN TENSORS MUTUALLY CONSISTENT - NEW PYEQ >= 0.50.3
    ##############################################################################################################################

    if model.debug: model_tmp = copy.deepcopy(model)

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
    MESSAGE("Model       will include %04d dates from %s to %s" % (model.np_model_date_s.shape[0], str_sd, str_ed))

    str_sd = at.seconds2datetime(model.np_obs_date_s[0]).isoformat(' ')
    str_ed = at.seconds2datetime(model.np_obs_date_s[-1]).isoformat(' ')
    MESSAGE("Observation will include %04d dates from %s to %s" % (model.np_model_date_s.shape[0], str_sd, str_ed))


    ##############################################################################################################################
    # BACKWARD CASE
    ##############################################################################################################################
    if model.backward:
        WARNING("Backward mode")
        WARNING("model.t_obs.shape %s" % str(model.t_obs.shape))
        model.t_obs[:,:, 0:3] = -np.flip( model.t_obs[:,:, 0:3] , axis=0 )
        model.t_obs[:,:,3:] = np.flip( model.t_obs[:,:,3:] , axis=0 )
        model.np_obs_date_s = np.append(model.np_obs_date_s[0], model.np_obs_date_s[0] - np.cumsum(np.diff(np.flip(model.np_obs_date_s))))
        model.np_model_date_s = np.append(model.np_model_date_s[0], model.np_model_date_s[0] - np.cumsum(np.diff(np.flip(model.np_model_date_s))))


    #### data so far is not yet set to zeros for the first obs, correct that
    # TODO IMPROVE THIS BY ADDING OPTIONS
    #WARNING("SET FIRST OBSERVATION AS ZERO")
    #model.warning += "SET FIRST OBSERVATION AS ZERO\n"
    #MESSAGE("Setting the first observation as the reference")
    #model.t_obs = pyeq.obs_tensor.set_zero_at_first_obs.set_zero_at_first_obs(model.t_obs, verbose=model.verbose)
    model.t_obs = pyeq.obs_tensor.set_zero_obs.set_zero_obs(model.t_obs, model.ts_zero)

    if model.debug: pyeq.log.print_diff_model_attributes(model_tmp, model)

    ##############################################################################################################################
    # RE-ORDER THE GREEN TENSOR ACCORDING TO THE MAIN AND CONJUGATE RAKE
    ##############################################################################################################################
    if model.debug: model_tmp = copy.deepcopy(model)

    # save original green tensor
    model.green_ds_ss = np.copy( model.green )

    MESSAGE("Reformatting GREEN tensor to account for the main and conjugate rake")

    # computing model.rake_subfault, model.slip_az, slip_dir_en
    import pyacs.lib.faultslip

    if model.rake_type.lower() == 'euler':
        model.rake_subfault = pyacs.lib.faultslip.rake_from_euler( model.sgeometry.centroid_long, model.sgeometry.centroid_lat, model.sgeometry.strike, model.sgeometry.dip, model.rake_value )

    if model.rake_type in ['constant','fixed']:
        model.rake_subfault = np.zeros( model.geometry.shape[0] ) + float(model.rake_value)

    if model.rake_type == 'file':
        MESSAGE( ("Reading rake from file: %s" % ( model.rake_value )))
        model.rake_subfault = np.genfromtxt(model.rake_value )[:,-1]
        # if model.rake_subfault.shape[0] != model.sgeometry.shape[0] then try to interpolate
        if model.rake_subfault.shape[0] != model.sgeometry.shape[0]:
            WARNING("rake file %s and geometry have different length. Will try to perform interpolation" % model.rake_value)
            rake_data = np.genfromtxt( model.rake_value )
            from scipy.interpolate import griddata

            model.rake_subfault = griddata( rake_data[:,1:3], rake_data[:,-1],
                                    (model.sgeometry.centroid_long, model.sgeometry.centroid_lat), method='nearest')

            if np.isnan( model.rake_subfault ).any():
                ERROR("Could not properly make interpolation. Check grid dimenensions:")
                ERROR("Bounds of model geometry : %8.2lf - %8.2lf / %8.2lf - %8.2lf" % \
                      ( np.min( model.sgeometry.centroid_long ) , \
                        np.max( model.sgeometry.centroid_long ) , \
                        np.min(model.sgeometry.centroid_lat), \
                        np.max(model.sgeometry.centroid_lat)
                        ))

                ERROR("Bounds of rake file      : %8.2lf - %8.2lf / %8.2lf - %8.2lf" % \
                      ( np.min( rake_data[:,1] ) , \
                        np.max( rake_data[:,1] ) , \
                        np.min( rake_data[:,2] ), \
                        np.max( rake_data[:,2] )
                        ),exit=True)
            else:
                WARNING("Successfully interpolated rake file on model geometry")
                WARNING("Rake range is : %8.2lf - %8.2lf  " % (np.min(model.rake_subfault) , np.max(model.rake_subfault)))

    model.slip_az = pyacs.lib.faultslip.strike_dip_rake_to_dir( model.sgeometry.strike, model.sgeometry.dip, model.rake_subfault)
    model.slip_dir_en = np.array( [np.sin( np.radians( model.slip_az )) ,  np.cos( np.radians( model.slip_az ) ) ] ).T


    model.green , model.rake = pyeq.elastic_tensor.green_ds_ss_to_main_rake_conjugate( \
        model.green, model.geometry, model.sgeometry, 'vector', model.rake_subfault)

    MESSAGE("Green tensor shape is now: %s" % str(model.green.shape))

    if model.debug: pyeq.log.print_diff_model_attributes(model_tmp, model)


    ##############################################################################################################################
    # GETTING MESH TOPOLOGY
    ##############################################################################################################################

    if model.geometry_type in ['TDE','TDV']:
        MESSAGE("Building mesh topology")
        import pyeq.lib.geometry.make_topology_trimesh
        model = pyeq.lib.geometry.make_topology_trimesh( model )

        MESSAGE("Mesh used for inversion has %d triangular subfaults and %d triangle vertices" % (model.nfaults,model.nvertices))

    ##############################################################################################################################
    # SOLVE AT NODES RATHER AT SUBFAULTS
    ##############################################################################################################################


    if model.geometry_type == 'TDV':
        MESSAGE("Triangular dislocation vertices (TDV) option")
        MESSAGE("Slip rate will be estimated on vertices")
        VERBOSE("Building subfault to node transformation matrix")
        model.matrix_subfault_to_node = np.zeros(( ( model.green.shape[1], model.topology.vertex_coor.shape[0] ) ))
        VERBOSE("model.matrix_subfault_to_node.shape %s" % str(model.matrix_subfault_to_node.shape))
        for i in np.arange( model.matrix_subfault_to_node.shape[0] ):
            model.matrix_subfault_to_node[i, model.topology.cell_vertex_idx[i,0]] = 1./3.
            model.matrix_subfault_to_node[i, model.topology.cell_vertex_idx[i,1]] = 1./3.
            model.matrix_subfault_to_node[i, model.topology.cell_vertex_idx[i,2]] = 1./3.
        VERBOSE("Applying subfault to vertice transformation matrix")
        model.green_subfault = np.copy( model.green )
        model.green_node = np.zeros(( model.green.shape[0], model.topology.vertex_coor.shape[0], 3, 2 ))
        model.green_node[:,:,0,0] = np.dot( model.green[:,:,0,0], model.matrix_subfault_to_node)
        model.green_node[:,:,0,1] = np.dot( model.green[:,:,0,1], model.matrix_subfault_to_node)
        model.green_node[:,:,1,0] = np.dot( model.green[:,:,1,0], model.matrix_subfault_to_node)
        model.green_node[:,:,1,1] = np.dot( model.green[:,:,1,1], model.matrix_subfault_to_node)
        model.green_node[:,:,2,0] = np.dot( model.green[:,:,2,0], model.matrix_subfault_to_node)
        model.green_node[:,:,2,1] = np.dot( model.green[:,:,2,1], model.matrix_subfault_to_node)
        VERBOSE("model.green_node has shape %s" % str(model.green_node.shape))
        VERBOSE("copying model.green_node to model.green")
        model.green = model.green_node

    ##############################################################################################################################
    # CROSS-VALIDATION CASE
    ##############################################################################################################################

    if model.cross_validation == 'build':
        MESSAGE("Building matrices for Cross-validation" )
        # divides obs in K folds K=9

        date_n3 = int( model.t_obs[:,0,0].size / 3  )
        r_n3 = model.t_obs[:,0,0].size % 3
        if r_n3 == 0:
            len_dates = [date_n3,date_n3,date_n3]
        if r_n3 == 1:
            len_dates = [date_n3,date_n3+1,date_n3]
        if r_n3 == 2:
            len_dates = [date_n3,date_n3+1,date_n3+1]

        gps_n3 = int( model.t_obs[0,:,0].size / 3  )
        r_n3 = model.t_obs[0,:,0].size % 3
        if r_n3 == 0:
            len_gps = [gps_n3,gps_n3,gps_n3]
        if r_n3 == 1:
            len_gps = [gps_n3, gps_n3+1, gps_n3]
        if r_n3 == 2:
            len_gps = [gps_n3,gps_n3+1,gps_n3+1]

        # array of indices of obs used for validation
        np_idxcv_dates = np.zeros(( 10,2 ),dtype=int)
        np_idxcv_gps   = np.zeros(( 10,2 ),dtype=int)

        # K1
        np_idxcv_dates[1,:] = [0,len_dates[0]]
        np_idxcv_gps[1,:]   = [0,len_gps[0]]
        # K2
        np_idxcv_dates[2,:] = [0,len_dates[0]]
        np_idxcv_gps[2,:]   = [len_gps[0],len_gps[0]+len_gps[1]]
        # etc K 3-> 9
        np_idxcv_dates[3,:] = [0,len_dates[0]]
        np_idxcv_gps[3,:]   = [len_gps[0]+len_gps[1],len_gps[0]+len_gps[1]+len_gps[2]]

        np_idxcv_dates[4,:] = [len_dates[0],len_dates[0]+len_dates[1]]
        np_idxcv_gps[4,:]   = [0,len_gps[0]]
        np_idxcv_dates[5,:] = [len_dates[0],len_dates[0]+len_dates[1]]
        np_idxcv_gps[5,:]   = [len_gps[0],len_gps[0]+len_gps[1]]
        np_idxcv_dates[6,:] = [len_dates[0],len_dates[0]+len_dates[1]]
        np_idxcv_gps[6,:]   = [len_gps[0]+len_gps[1],len_gps[0]+len_gps[1]+len_gps[2]]

        np_idxcv_dates[7,:] = [len_dates[0]+len_dates[1],len_dates[0]+len_dates[1]+len_dates[2]]
        np_idxcv_gps[7,:]   = [0,len_gps[0]]
        np_idxcv_dates[8,:] = [len_dates[0]+len_dates[1],len_dates[0]+len_dates[1]+len_dates[2]]
        np_idxcv_gps[8,:]   = [len_gps[0],len_gps[0]+len_gps[1]]
        np_idxcv_dates[9,:] = [len_dates[0]+len_dates[1],len_dates[0]+len_dates[1]+len_dates[2]]
        np_idxcv_gps[9,:]   = [len_gps[0]+len_gps[1],len_gps[0]+len_gps[1]+len_gps[2]]

        MESSAGE("Saving indices of K-folds")
        np.save( 'CROSS_VALIDATION_idxcv_dates.npy',np_idxcv_dates )
        np.save( 'CROSS_VALIDATION_idxcv_gps.npy',np_idxcv_gps )

        MESSAGE("Building and Saving Normal Systems for K-fold cross-validation")
        # now these are firmly defined

        model.nfaults = model.green.shape[1]
        model.nstep = model.np_model_date_s.shape[0] - 1
        model.nparameters = model.nfaults * model.nstep

        if model.up:
            model.ncomponent = 3
        else:
            model.ncomponent = 2

        MESSAGE("Use up component: %s" % model.up)
        MESSAGE("number of faults/vertices: %d" % model.nfaults)
        MESSAGE("number of model time steps: %d" % model.nstep)
        MESSAGE("number of model parameters: %d" % model.nparameters)

        MESSAGE("Building full observation system")
        (model.N, model.Nd) = pyeq.forward_model.build5(model)
        MESSAGE("Saving CROSS_VALIDATION_N.npy & CROSS_VALIDATION_Nd.npy")
        np.save('CROSS_VALIDATION_N.npy',model.N)
        np.save('CROSS_VALIDATION_Nd.npy',model.Nd)

        # save model.t_obs
        save_t_obs = np.copy( model.t_obs )

        for i in np.arange(1,10):
            model.t_obs = np.copy( save_t_obs )
            MESSAGE("Building observation system for K-fold : %d" % i)

            # change model.t_obs
            model.t_obs[np_idxcv_dates[i,0]:np_idxcv_dates[i,1],np_idxcv_gps[i,0]:np_idxcv_gps[i,1],:] = np.nan
            MESSAGE("Building normal system for K-fold: %d" % i)
            (model.N, model.Nd) = pyeq.forward_model.build5(model)
            npyN = ("CROSS_VALIDATION_NK%d.npy" % i)
            npyNd = ("CROSS_VALIDATION_NdK%d.npy" % i)
            MESSAGE("Saving %s & %s" % (npyN,npyNd))
            np.save(npyN,model.N)
            np.save(npyNd,model.Nd)

        MESSAGE("Renaming output")
        model.name = 'CROSS_VALIDATION'
        model.t_obs = np.copy(save_t_obs)

    ##############################################################################################################################
    # OBSERVATION NORMAL SYSTEM
    ##############################################################################################################################

    if model.debug: model_tmp = copy.deepcopy(model)

    MESSAGE("BUILDING OBSERVATION NORMAL SYSTEM", level=1)

    model.memory_before_ons = pyeq.log.get_process_memory_usage()

    # now these are firmly defined

    model.nfaults = model.green.shape[1]
    model.nstep = model.np_model_date_s.shape[0] - 1
    model.nparameters = model.nfaults * model.nstep

    if model.up:
        model.ncomponent = 3
    else:
        model.ncomponent = 2

    MESSAGE("Use up component: %s" % model.up)
    MESSAGE("number of faults/vertices: %d" % model.nfaults)
    MESSAGE("number of model time steps: %d" % model.nstep)
    MESSAGE("number of model parameters: %d" % model.nparameters)

    T0 = time()
    MESSAGE("Building observation normal system. Can be long...")
    if model.build == 0:
        WARNING("This is obsolete routine. Building normal system using pyeq.forward_model.build0")
        model.normal_system_build = 'pyeq.forward_model.build.G_d'
        (G1, d1, diag_Cd1) = pyeq.lib.build.G_d(model, tol_date=5, window=0.1)
        # still needs to build ATPA / ATPB from G1 & d1

    if model.build == 1:
        WARNING("This is obsolete routine. Building normal system using pyeq.forward_model.build1")
        model.normal_system_build = 'pyeq.forward_model.build1'
        (ATPA_NEW, ATPB_NEW) = pyeq.forward_model.build1(model, date_tol=0.00000001)

    if model.build == 2:
        WARNING("This is obsolete routine. Building normal system using pyeq.forward_model.build2")
        model.normal_system_build = 'pyeq.forward_model.build2'
        (ATPA_NEW, ATPB_NEW) = pyeq.forward_model.build2(model)

    if model.build == 3:
        WARNING("This is obsolete routine. Building normal system using pyeq.forward_model.build3")
        (model.N, model.Nd) = pyeq.forward_model.build3(model)
        model.nconstant = 0

    if model.build == 4:
        WARNING("This is obsolete routine. Building normal system using pyeq.lib.build4.build4")
        RAKE = False
        model.normal_system_build = 'pyeq.forward_model.build4'
        (model.N, model.Nd) = pyeq.forward_model.build4(model)
        # N = model.N
        # Nd = model.Nd

    if model.build == 5:
        DEBUG("Building normal system using pyeq.forward_model.build5")
        RAKE = False
        model.normal_system_build = 'pyeq.forward_model.build5'
        (model.N, model.Nd) = pyeq.forward_model.build5(model)

    model.time_build = time() - T0
    MESSAGE(("time building the observation normal system: %.1lf s" % (model.time_build)))
    MESSAGE(("normal system shape: %s" % (model.N.shape,)))
    MESSAGE(("current used memory size (Gb): %.2lf" % (pyeq.log.get_process_memory_usage())))

    if model.debug: pyeq.log.print_diff_model_attributes(model_tmp, model)
    if model.debug: pyeq.log.display_array(model.N * 1.E6)

    DEBUG("model.green  shape: %s" % (model.green.shape,))
    DEBUG("model.N  shape: %s " % (model.N.shape,))
    DEBUG("model.N without constants: %d " % model.nfaults * model.nstep)
    DEBUG("model.Nd shape: %s " % (model.Nd.shape,))
    DEBUG("model.nfaults: %d " % model.nfaults)
    DEBUG("model.nstep: %d " % model.nstep)
    DEBUG("n constant: %d " % (model.N.shape[0] - (model.nstep * model.nfaults)))
    DEBUG(
        "n constant vs n gps: %d vs %d " % ((model.N.shape[0] - model.nstep * model.nfaults) / 2, model.green.shape[0]))

    ##############################################################################################################################
    # SAVE OPTION
    ##############################################################################################################################

    if model.save:
        MESSAGE("SAVING OBSERVATION NORMAL SYSTEM AS PICKLE (.MPCK)", level=1)
        import pyeq.log

        pyeq.log.save_model_N(model)

##############################################################################################################################
# CROSS_VALIDATION_RUN
##############################################################################################################################

if model.cross_validation == 'run':
    MESSAGE("Cross validation run", level=1)
    MESSAGE("regularization parameters: %8.3lf %8.3lf %8.3lf %8.3lf %8.3lf " % ( model.sigma, model.lambda_spatial_smoothing,
                model.lambda_spatial_smoothing, model.lambda_final_spatial_smoothing, model.lambda_stf ))

    # load model


    # full data run
    model.N = np.load( 'CROSS_VALIDATION_N.npy' )
    model.Nd = np.load( 'CROSS_VALIDATION_Nd.npy' )

    # adds regularization
    model = pyeq.regularization.laplace.add_laplace_cons(model)
    model = pyeq.regularization.damping.add_damping_cons(model)

    if model.offset > 0:
        model.shift_constant = -10  # 10 mm should be enough
        MESSAGE("modifying the normal system to handle offset at the origin time")
        model.Nd += np.dot(model.N[:, -model.nconstant:], np.ones((model.nconstant, 1)) * -model.shift_constant).flatten()
    else:
        model.shift_constant = 0.

    ###################################################################
    # MODIFIES THE NORMAL SYSTEM SO THAT NEGATIVE VALUES FOR SOME
    # PARAMETERS, OR CHANGE THE LOWER BOUND.
    # HANDLES THE CASE OF NEGATIVE (BACK-SLIP BOUNDED) SLIP
    # FOR INTERSEISMIC MODELLING
    ###################################################################

    if model.interseismic != 0:
        MESSAGE(("INTERSEISMIC CASE"), level=1)
        ERROR("NOT IMPLEMENTED YET", exit=True)
        print("-- modifying the normal system to allow negative bounded slip for the interseismic case")
        if model.nconstant == 0:
            model.Nd += np.dot(model.N,
                               np.ones(((model.nfaults * model.nstep), 1)) * -model.interseismic / 365.25).flatten()
        else:
            model.Nd += np.dot(model.N[:, :-model.nconstant],
                               np.ones(((model.nfaults * model.nstep), 1)) * -model.interseismic / 365.25).flatten()

    else:
        MESSAGE("TRANSIENT SLIP CASE", level=1)

    ###################################################################
    # RUNNING INVERSION - NO RAKE
    ###################################################################
    MESSAGE("MAKING INVERSION", level=1)

    model.memory_before_inv = pyeq.log.get_process_memory_usage()
    MESSAGE("memory usage before inversion: %.2lf Gb" % model.memory_before_inv)
    MESSAGE("Number of estimated parameters: %d" % (model.nparameters))

    T0 = time()
    if model.regularization in ['covariance', 'laplacian']:
        model.parameters, time_inversion = pyeq.optimization.wrapper.make_inversion.pyeq_nnls(model.N, model.Nd, model.nnls,
                                                                                              verbose=model.verbose)

    # NOT IMPLEMENTED YET
    if model.regularization in ['cvxopt']:
        model.parameters, time_inversion = pyeq.optimization.wrapper.make_inversion.pyeq_cvxopt(model)

    model.time_inversion = time() - T0
    VERBOSE(("time inversion: %.1lf s" % (model.time_inversion)))




##############################################################################################################################
# PART FOR REGULAR & RUN FROM MODEL.PCK
##############################################################################################################################


##############################################################################################################################
# SLIP BOUNDS & CONSTRAINTS
# NOT (RE)-IMPLEMENTED YET
##############################################################################################################################
if model.debug: model_tmp = copy.deepcopy(model)
# model = pyeq.lib.regularization.bounds()


###########################################################################
# START REGULARIZATION PART
###########################################################################

MESSAGE(("BUILDING REGULARIZATION NORMAL SYSTEM: %s  " % (model.regularization.upper())), level=1)

###########################################################################
# GET SIGMA INFO
###########################################################################
# if model.regularization == 'covariance':
VERBOSE("deciphering sigma argument")
model = pyeq.regularization.damping.decipher_sigma_arg(model)

T0 = time()

model.memory_before_rns = pyeq.log.get_process_memory_usage()

##############################################################################################################################
# LAPLACIAN REGULARIZATION
##############################################################################################################################

if model.regularization == 'laplacian':
    import pyeq.regularization.laplace

    model = pyeq.regularization.laplace.add_laplace_cons(model)
    model = pyeq.regularization.damping.add_damping_cons(model)

##############################################################################################################################
# COVARIANCE REGULARIZATION
##############################################################################################################################

# COVARIANCE
if model.regularization == 'covariance':
    model = pyeq.lib.regularization.make_model_covariance_03(model)

##############################################################################################################################
# OBSOLETE OR NOT FINISHED IMPLEMENTATIONS
##############################################################################################################################

# This is the old Laplacian
# if model.regularization == 'laplacian':
#    #model = pyeq.lib.regularization.make_regularization_laplacian( model )
#    model = pyeq.lib.regularization.make_regularization_valette_01( model )
# LAPLACIAN_LIKE
# if model.regularization == 'laplacian_like':
#    model = pyeq.lib.regularization.make_regularization_laplacian_like( model )
# CVXOPT
# if model.regularization == 'cvxopt':
#    model = pyeq.lib.regularization.get_regularization_operators( model )


model.time_regularization = time() - T0
VERBOSE(("time building regularization normal system: %.1lf s" % (model.time_regularization)))
VERBOSE("memory usage: %.2lf Gb" % pyeq.log.get_process_memory_usage())

if model.debug: pyeq.log.print_diff_model_attributes(model_tmp, model)

model.memory_before_mns = pyeq.log.get_process_memory_usage()
model.time_merging_obs_regularization = time() - T0

# no prior value for model
# print("-- Adding observation and normalization RHS vector")
# ATB = ATPB_NEW

if model.debug: pyeq.log.print_diff_model_attributes(model_tmp, model)

model.N_size = model.N.nbytes / 1024 / 1024 / 1024
VERBOSE("Estimated normal matrix memory size: %.2lf Gb" % model.N_size)

###################################################################
# MODIFIES THE NORMAL SYSTEM SO THAT NEGATIVE VALUES FOR SOME
# PARAMETERS ARE ALLOWED, OR CHANGE THE LOWER BOUND.
# HANDLES ORIGIN TIME OFFSET CASE
###################################################################

if model.offset > 0:
    model.shift_constant = -10  # 10 mm should be enough
    MESSAGE("modifying the normal system to handle offset at the origin time")
    model.Nd += np.dot(model.N[:, -model.nconstant:], np.ones((model.nconstant, 1)) * -model.shift_constant).flatten()
else:
    model.shift_constant = 0.

###################################################################
# MODIFIES THE NORMAL SYSTEM SO THAT NEGATIVE VALUES FOR SOME
# PARAMETERS, OR CHANGE THE LOWER BOUND.
# HANDLES THE CASE OF NEGATIVE (BACK-SLIP BOUNDED) SLIP
# FOR INTERSEISMIC MODELLING
###################################################################

if model.interseismic != 0:
    MESSAGE(("INTERSEISMIC CASE"), level=1)
    ERROR("NOT IMPLEMENTED YET", exit=True)
    print("-- modifying the normal system to allow negative bounded slip for the interseismic case")
    if model.nconstant == 0:
        model.Nd += np.dot(model.N,
                           np.ones(((model.nfaults * model.nstep), 1)) * -model.interseismic / 365.25).flatten()
    else:
        model.Nd += np.dot(model.N[:, :-model.nconstant],
                           np.ones(((model.nfaults * model.nstep), 1)) * -model.interseismic / 365.25).flatten()

else:
    MESSAGE("TRANSIENT SLIP CASE", level=1)

###################################################################
# NO_OPT CASE
###################################################################

if model.no_opt:
    MESSAGE("no_opt option set to True. Stopping before optimization")
    sys.exit()

###################################################################
# THIS A TEST FOR THE RAKE VARIABLE INVERSION
###################################################################

if model.rake_constraint[0] != '0':

    MESSAGE("SEARCHING FOR BEST SUBFAULT VARIABLE SLIP DIRECTION" , level=1 )

    model.memory_before_inv = pyeq.log.get_process_memory_usage()
    MESSAGE("memory usage before inversion: %.2lf Gb" % model.memory_before_inv)
    MESSAGE("Number of estimated parameters: %d" % (model.nparameters))

    T0 = time()

    import pyeq.optimization.rake.subfault_variable
    model = pyeq.optimization.rake.subfault_variable( model )

    model.time_inversion = time() - T0
    VERBOSE(("time inversion: %.1lf s" % (model.time_inversion)))

else:

    # no log for variable rake
    model.log_variable_rake = None

    ###################################################################
    # RUNNING INVERSION - NO RAKE
    ###################################################################
    MESSAGE("MAKING INVERSION", level=1)

    model.memory_before_inv = pyeq.log.get_process_memory_usage()
    MESSAGE("memory usage before inversion: %.2lf Gb" % model.memory_before_inv)
    MESSAGE("Number of estimated parameters: %d" % (model.nparameters))

    T0 = time()
    if model.regularization in ['covariance', 'laplacian']:
        model.parameters, time_inversion = pyeq.optimization.wrapper.make_inversion.pyeq_nnls(model.N, model.Nd, model.nnls,
                                                                                              verbose=model.verbose)

    # NOT IMPLEMENTED YET
    if model.regularization in ['cvxopt']:
        model.parameters, time_inversion = pyeq.optimization.wrapper.make_inversion.pyeq_cvxopt(model)

    model.time_inversion = time() - T0
    VERBOSE(("time inversion: %.1lf s" % (model.time_inversion)))

    # NOT IMPLEMENTED YET
    if model.interseismic != 0:
        model.parameters[: model.nfaults * model.nstep] = model.parameters[
                                                          : model.nfaults * model.nstep] + model.interseismic / 365.25

    ###################################################################
    # BACKWARD - SET IT BACKWARD AGAIN TO GET FORWARD
    ###################################################################
    if model.backward:
        WARNING("Backward mode - Setting back to forward")
        WARNING("model.t_obs.shape %s" % str(model.t_obs.shape))
        model.t_obs[:,:, 0:3] = -np.flip( model.t_obs[:,:, 0:3] , axis=0 )
        model.t_obs[:,:,3:] = np.flip( model.t_obs[:,:,3:] , axis=0 )
        model.np_obs_date_s = np.append(model.np_obs_date_s[0], model.np_obs_date_s[0] - np.cumsum(np.diff(np.flip(model.np_obs_date_s))))
        model.np_model_date_s = np.append(model.np_model_date_s[0], model.np_model_date_s[0] - np.cumsum(np.diff(np.flip(model.np_model_date_s))))

        model.parameters = np.flip(model.parameters.reshape( -1,model.green.shape[1] ).T, axis=1).T.flatten()


###################################################################
# REFORMATING IF TDV
###################################################################

if model.geometry_type == 'TDV':
    MESSAGE("Reformatting output slip from triangular vertices (TDV) to triangular dislocation elements (TDE)")
    SRN =  model.parameters.reshape( -1,model.green.shape[1] ).T
    model.green = model.green_subfault
    SR = np.zeros((model.matrix_subfault_to_node.shape[0],SRN.shape[1]))

    for i in np.arange( SRN.shape[1] ):
        SR[:,i] = np.dot( model.matrix_subfault_to_node,SRN[:,i])
    model.parameters = SR.T.flatten()
    model.nfaults = model.green.shape[1]

if pyacs.debug():
    pyeq.log.print_model_attributes( model )


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
    delattr(model, 'N')
    delattr(model, 'Nd')
except:
    ERROR("deleting normal system attributes")

###################################################################
# NAME OF THE RUN - THINK ABOUT BETTER NAMING
###################################################################
if model.name is None:
    model.name = pyeq.conf.make_model_name(model)

model.odir = model.name
# TODO improve name with 0 instead of E
MESSAGE("run name: %s" % model.name)


###################################################################
# INVERSION SUMMARY
###################################################################
MESSAGE(("INVERSION RESULTS SUMMARY"), level=1)
if model.nconstant > 0:
    model.slip = model.parameters[:-model.nconstant]
    model.estimated_offsets = model.parameters[-model.nconstant:]
else:
    model.slip = model.parameters
    model.estimated_offsets = None

import pyeq.log
model = pyeq.log.print_dates( model, save=False )
model = pyeq.log.print_sol_to_slip( model, save=False  )
model = pyeq.log.print_stf( model, save=False  )
model.M0 = model.CSTF[-1]
model = pyeq.log.print_offset( model, save=False  )
model = pyeq.log.print_modeled_time_series( model, save=False  )

if model.M0 > 0:
    model.magnitude = 2./3.*(np.log10( model.M0 )-9.1)
else:
    model.magnitude = -99
# moment
MESSAGE("Moment (N.n): %8.2E (Mw%.1lf)" % (model.M0,model.magnitude))
# max cumulated slip

# maximum moment rate
mmr = np.max(model.STF)
if mmr>0:
    MESSAGE("Maximum moment rate per day: %8.2E (Mw%.1lf)" % (mmr,2./3.*(np.log10( mmr )-9.1)))
else:
    MESSAGE("Maximum moment rate per day: %8.2E (Mw%.1lf)" % (mmr,-99))

# max slip rate
i_cs_max = np.unravel_index(np.argmax(model.NORM_CUMULATIVE_SLIP_PER_TIME_STEP),
                            model.NORM_CUMULATIVE_SLIP_PER_TIME_STEP.shape)
slip_max = model.NORM_CUMULATIVE_SLIP_PER_TIME_STEP[i_cs_max]
model.idx_cumulative_slip_max = i_cs_max[1]

MESSAGE("cumulative slip max  (mm)        : %8.1f fault #%d (%.2lf,%.2lf)" % (slip_max, i_cs_max[1],
                                                                                model.sgeometry.centroid_long[i_cs_max[1]],
                                                                                model.sgeometry.centroid_lat[i_cs_max[1]]))

# rms
def make_rms_model( obs , model ):
    # takes displacement tensors and computes basic stats
    # corrects model so that 0 is at the first epoch of obs
    import numpy as np
    # find and correct bias for each site i
    mmodel = np.copy(model)
    for i in np.arange(obs.shape[1]):
        idx = np.where(~np.isnan(obs[:,i,0]))[0][0]
        mmodel[:,i,:] = mmodel[:,i,:] - model[idx,i,:]
    # residuals
    res = obs[:,:,:3] - mmodel
    # computes statistics
    rms_h = np.sqrt( np.nanmean(res[:,:,:2]**2) )
    rms_v = np.sqrt( np.nanmean(res[:,:,2]**2) )
    chi2 = np.nansum( (res[:,:,:3]/obs[:,:,3:])**2 )
    num =  np.nansum( (1./obs[:,:,3:])**2 )
    wrms_3D = np.sqrt( chi2/num )

    return rms_h, rms_v, chi2, wrms_3D

rms_h, rms_v, chi2, wrms_3d = make_rms_model(model.t_obs, model.t_mod)

MESSAGE(("rms_h (mm):%.2lf"% rms_h))
MESSAGE(("rms_v (mm):%.2lf"% rms_v))
MESSAGE(("wrms_3d (mm):%.2lf"% wrms_3d))

# cross validation

# if model.cross_validation is not None:
#     # compute MSE
#     SE = np.nansum((model.t_obs[idxcv_dates[0]:idxcv_dates[1],idxcv_gps[0]:idxcv_gps[1], :3]
#           - model.t_mod[idxcv_dates[0]:idxcv_dates[1],idxcv_gps[0]:idxcv_gps[1], :3]) ** 2
#             / model.t_obs[idxcv_dates[0]:idxcv_dates[1],idxcv_gps[0]:idxcv_gps[1], 3:] ** 2)
#     MSE = SE / ( model.t_obs[idxcv_dates[0]:idxcv_dates[1],idxcv_gps[0]:idxcv_gps[1], :3].size - np.count_nonzero(np.isnan(model.t_obs[idxcv_dates[0]:idxcv_dates[1],idxcv_gps[0]:idxcv_gps[1], :3]))  )
#
#     SQRT_MSE = np.sqrt( MSE )
#     MESSAGE("MSE for %s: %.4lf %.3lf mm" % ( model.cross_validation,MSE, SQRT_MSE ))

###################################################################
# PRINT RESULTS OR JUST SAVE MODEL.PCK
###################################################################

if model.print_result:

    MESSAGE("PRINTING RESULTS", level=1)
    pyeq.log.print_results(model)

else:
    model_pck = 'model_' + model.name + '.mpck'
    MESSAGE(("SAVING MODEL AS PICKLE (.mpck): %s" % model_pck), level=1)
    model.write_pickle(model_pck)

###################################################################
# MAKE PLOTS
###################################################################
if model.plot:
    MESSAGE("MAKING PLOTS", level=1)

    import pyeq.plot

    pyeq.plot.make_plot(model)

###################################################################
# MAKE A COMPRESSED TAR FILE
###################################################################

if model.tar:
    MESSAGE(("MAKING A COMPRESSED TAR ARCHIVE"), level=1)

    import traceback
    import subprocess
    from datetime import datetime

    MESSAGE("creating %s.tar.bz2" % model.odir)

    try:
        str_time = s1 = datetime.now().strftime("-%Y-%m-%d-%H-%M-%S")
        cmd = ['tar', 'cfj', model.odir + str_time + '.tar.bz2', model.odir]
        output = subprocess.check_output(cmd).decode("utf-8").strip()
        DEBUG(output)
        VERBOSE("removing directory: %s" % model.odir)
        cmd = ['rm', '-Rf', model.odir]
        output = subprocess.check_output(cmd).decode("utf-8").strip()
        print(output)

    except Exception:
        print(f"E: {traceback.format_exc()}")
