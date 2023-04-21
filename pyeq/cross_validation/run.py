def run( model ):
    """
    runs K-fold cross_validation
    requires to have run pyeq.cross_validation.build first

    pyeq.cross_validation.build should have created the following files:
    CROSS_VALIDATION_N.npy
    CROSS_VALIDATION_NK1.npy
    CROSS_VALIDATION_NK2.npy
    CROSS_VALIDATION_NK3.npy
    CROSS_VALIDATION_NK4.npy
    CROSS_VALIDATION_NK5.npy
    CROSS_VALIDATION_NK6.npy
    CROSS_VALIDATION_NK7.npy
    CROSS_VALIDATION_NK8.npy
    CROSS_VALIDATION_NK9.npy
    CROSS_VALIDATION_Nd.npy
    CROSS_VALIDATION_NdK1.npy
    CROSS_VALIDATION_NdK2.npy
    CROSS_VALIDATION_NdK3.npy
    CROSS_VALIDATION_NdK4.npy
    CROSS_VALIDATION_NdK5.npy
    CROSS_VALIDATION_NdK6.npy
    CROSS_VALIDATION_NdK7.npy
    CROSS_VALIDATION_NdK8.npy
    CROSS_VALIDATION_NdK9.npy
    CROSS_VALIDATION_idxcv_dates.npy
    CROSS_VALIDATION_idxcv_gps.npy
    model_CROSS_VALIDATION.mpck
    """

    # import
    import pyeq.message.message as MESSAGE
    import pyeq.message.verbose_message as VERBOSE
    import pyeq.message.error as ERROR
    import pyeq.message.warning as WARNING
    import pyeq.message.debug_message as DEBUG

    # residuals routine
    def make_res_model(obs, model):
        # takes displacement tensors and computes basic stats
        # corrects model so that 0 is at the first epoch of obs
        import numpy as np
        # find and correct bias for each site i
        mmodel = np.copy(model)
        for i in np.arange(obs.shape[1]):
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

    # load index info
    np_idxcv_dates = np.load('CROSS_VALIDATION_idxcv_dates.npy')
    np_idxcv_gps   = np.load('CROSS_VALIDATION_idxcv_gps.npy')


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

    model = pyeq.regularization.damping.decipher_sigma_arg(model)

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

    # computes rms
    rms_h, rms_v, chi2, wrms_3d = make_rms_model(model.t_obs, model.t_mod)

    # computes reference statistics
    res = make_res_model(model.t_obs, model.t_mod)

    # computes statistics
    kstd_h = np.zeros(10)
    kstd_v = np.zeros(10)
    krms_h = np.zeros(10)
    krms_v = np.zeros(10)
    kchi2 = np.zeros(10)
    knum = np.zeros(10)
    kwrms_3D = np.zeros(10)
    for k in np.arange(1,10):
        kstd_h[k] = np.nanstd(res[np_idxcv_dates[k, 0]:np_idxcv_dates[k, 1], np_idxcv_gps[k, 0]:np_idxcv_gps[k, 1], :2])
        kstd_v[k] = np.nanstd(res[np_idxcv_dates[k, 0]:np_idxcv_dates[k, 1], np_idxcv_gps[k, 0]:np_idxcv_gps[k, 1], 2])
        krms_h[k] = np.sqrt(
            np.nanmean(res[np_idxcv_dates[k, 0]:np_idxcv_dates[k, 1], np_idxcv_gps[k, 0]:np_idxcv_gps[k, 1], :2]**2))
        krms_v[k] = np.sqrt(
            np.nanmean(res[np_idxcv_dates[k, 0]:np_idxcv_dates[k, 1], np_idxcv_gps[k, 0]:np_idxcv_gps[k, 1], 2]**2))
        kchi2[k] = np.nansum((res[np_idxcv_dates[k, 0]:np_idxcv_dates[k, 1], np_idxcv_gps[k, 0]:np_idxcv_gps[k, 1],
                          :3] / model.t_obs[np_idxcv_dates[k, 0]:np_idxcv_dates[k, 1],
                                np_idxcv_gps[k, 0]:np_idxcv_gps[k, 1], 3:])**2)
        knum[k] = np.nansum(
            (1. / model.t_obs[np_idxcv_dates[k, 0]:np_idxcv_dates[k, 1], np_idxcv_gps[k, 0]:np_idxcv_gps[k, 1], 3:])**2)
        kwrms_3D[k] = np.sqrt(kchi2[k] / knum[k])

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


    MSE_CV = 0.
    MSE_CV_STDH = 0.
    MSE_CV_RMSH = 0.

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

        model = pyeq.regularization.damping.decipher_sigma_arg(model)
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

        # computes residuals
        res = make_res_model(model.t_obs, model.t_mod)
        # extract relevant residuals
        rres = res[np_idxcv_dates[k,0]:np_idxcv_dates[k,1],np_idxcv_gps[k,0]:np_idxcv_gps[k,1],:]
        # extract relevant obs component std
        robs = model.t_obs[np_idxcv_dates[k,0]:np_idxcv_dates[k,1],np_idxcv_gps[k,0]:np_idxcv_gps[k,1], 3:]
        # set NaN where rres is NaN
        robs = rres*robs*0. + robs

        # computes statistics
        std_h = np.nanstd(rres[:,:, :2] )
        std_v = np.nanstd(rres[:,:, 2] )
        rms_h = np.sqrt(np.nanmean(rres[:,:, :2]**2))
        rms_v = np.sqrt(np.nanmean(rres[:,:, 2]**2))
        chi2 = np.nansum    ( ( rres[:,:, :3] / robs )**2)
        num = np.nansum((1. / robs )**2)
        wrms_3D = np.sqrt(chi2 / num)

        MESSAGE(("    MSE K=%d (mm) rms_h rms_v std_h std_v wrms_3d chi2" % k))
        MESSAGE("                  %.2lf  %.2lf  %.2lf  %.2lf    %.2lf  %.2lf" % (rms_h, rms_v, std_h, std_v, wrms_3D, chi2))
        MESSAGE(("REF MSE K=%d (mm) rms_h rms_v std_h std_v wrms_3d chi2" % k))
        MESSAGE("                  %.2lf  %.2lf  %.2lf  %.2lf    %.2lf  %.2lf" % (krms_h[k], krms_v[k], kstd_h[k], kstd_v[k], kwrms_3D[k], kchi2[k]))

        with open('MSE_K.dat','a') as flcurve:
            flcurve.write("%8.3lf %8.3lf %8.3lf %8.3lf %8.3lf %d %8.3lf %8.3lf %8.3lf %8.3lf %8.3lf %8.3lf\n" % \
                    tuple(map(float,(model.sigma, model.lambda_spatial_smoothing, \
                    model.lambda_temporal_smoothing, model.lambda_final_spatial_smoothing, model.lambda_stf, \
                     k, rms_h, rms_v, chi2, wrms_3d, std_h, std_v))))

        MSE_CV = MSE_CV + wrms_3d**2
        MSE_CV_STDH = MSE_CV_STDH + std_h**2
        MSE_CV_RMSH = MSE_CV_RMSH + rms_h**2

    with open('MSE_CV.dat','a') as flcurve:
        flcurve.write("%8.3lf %8.3lf %8.3lf %8.3lf %8.3lf %8.3lf %8.3lf %8.3lf\n" % \
                tuple(map(float,(model.sigma, model.lambda_spatial_smoothing, \
                model.lambda_temporal_smoothing, model.lambda_final_spatial_smoothing, model.lambda_stf, \
                 MSE_CV, MSE_CV_STDH,MSE_CV_RMSH))))


