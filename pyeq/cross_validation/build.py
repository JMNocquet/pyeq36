def build( model ):
    """
    Creates the files required for subsequent cross validation

    """

    # import
    import numpy as np
    import pyeq.forward_model
    import pyeq.message.message as MESSAGE
    import pyeq.message.verbose_message as VERBOSE
    import pyeq.message.error as ERROR
    import pyeq.message.warning as WARNING
    import pyeq.message.debug_message as DEBUG


    MESSAGE("Building matrices for Cross-validation")
    # divides obs in K folds K=9

    date_n3 = int(model.t_obs[:, 0, 0].size / 3)
    r_n3 = model.t_obs[:, 0, 0].size % 3
    if r_n3 == 0:
        len_dates = [date_n3, date_n3, date_n3]
    if r_n3 == 1:
        len_dates = [date_n3, date_n3 + 1, date_n3]
    if r_n3 == 2:
        len_dates = [date_n3, date_n3 + 1, date_n3 + 1]

    gps_n3 = int(model.t_obs[0, :, 0].size / 3)
    r_n3 = model.t_obs[0, :, 0].size % 3
    if r_n3 == 0:
        len_gps = [gps_n3, gps_n3, gps_n3]
    if r_n3 == 1:
        len_gps = [gps_n3, gps_n3 + 1, gps_n3]
    if r_n3 == 2:
        len_gps = [gps_n3, gps_n3 + 1, gps_n3 + 1]

    # array of indices of obs used for validation
    np_idxcv_dates = np.zeros((10, 2), dtype=int)
    np_idxcv_gps = np.zeros((10, 2), dtype=int)

    # K1
    np_idxcv_dates[1, :] = [0, len_dates[0]]
    np_idxcv_gps[1, :] = [0, len_gps[0]]
    # K2
    np_idxcv_dates[2, :] = [0, len_dates[0]]
    np_idxcv_gps[2, :] = [len_gps[0], len_gps[0] + len_gps[1]]
    # etc K 3-> 9
    np_idxcv_dates[3, :] = [0, len_dates[0]]
    np_idxcv_gps[3, :] = [len_gps[0] + len_gps[1], len_gps[0] + len_gps[1] + len_gps[2]]

    np_idxcv_dates[4, :] = [len_dates[0], len_dates[0] + len_dates[1]]
    np_idxcv_gps[4, :] = [0, len_gps[0]]
    np_idxcv_dates[5, :] = [len_dates[0], len_dates[0] + len_dates[1]]
    np_idxcv_gps[5, :] = [len_gps[0], len_gps[0] + len_gps[1]]
    np_idxcv_dates[6, :] = [len_dates[0], len_dates[0] + len_dates[1]]
    np_idxcv_gps[6, :] = [len_gps[0] + len_gps[1], len_gps[0] + len_gps[1] + len_gps[2]]

    np_idxcv_dates[7, :] = [len_dates[0] + len_dates[1], len_dates[0] + len_dates[1] + len_dates[2]]
    np_idxcv_gps[7, :] = [0, len_gps[0]]
    np_idxcv_dates[8, :] = [len_dates[0] + len_dates[1], len_dates[0] + len_dates[1] + len_dates[2]]
    np_idxcv_gps[8, :] = [len_gps[0], len_gps[0] + len_gps[1]]
    np_idxcv_dates[9, :] = [len_dates[0] + len_dates[1], len_dates[0] + len_dates[1] + len_dates[2]]
    np_idxcv_gps[9, :] = [len_gps[0] + len_gps[1], len_gps[0] + len_gps[1] + len_gps[2]]

    MESSAGE("Saving indices of K-folds")
    np.save('CROSS_VALIDATION_idxcv_dates.npy', np_idxcv_dates)
    np.save('CROSS_VALIDATION_idxcv_gps.npy', np_idxcv_gps)

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
    np.save('CROSS_VALIDATION_N.npy', model.N)
    np.save('CROSS_VALIDATION_Nd.npy', model.Nd)

    # save model.t_obs
    save_t_obs = np.copy(model.t_obs)

    for i in np.arange(1, 10):
        model.t_obs = np.copy(save_t_obs)
        MESSAGE("Building observation system for K-fold : %d" % i)

        # change model.t_obs
        model.t_obs[np_idxcv_dates[i, 0]:np_idxcv_dates[i, 1], np_idxcv_gps[i, 0]:np_idxcv_gps[i, 1], :] = np.nan
        MESSAGE("Building normal system for K-fold: %d" % i)
        (model.N, model.Nd) = pyeq.forward_model.build5(model)
        npyN = ("CROSS_VALIDATION_NK%d.npy" % i)
        npyNd = ("CROSS_VALIDATION_NdK%d.npy" % i)
        MESSAGE("Saving %s & %s" % (npyN, npyNd))
        np.save(npyN, model.N)
        np.save(npyNd, model.Nd)

    MESSAGE("Renaming output")
    model.name = 'CROSS_VALIDATION'
    model.odir = model.name
    model.t_obs = np.copy(save_t_obs)

    #return model

    # del N & Nd attributes
    VERBOSE("Deleting Normal System from model")

    try:
       delattr(model, 'N')
       delattr(model, 'Nd')
    except:
       ERROR("deleting normal system attributes")

    # save model_CROSS_VALIDATION.mpck
    model_pck = 'model_' + model.name + '.mpck'
    MESSAGE(("SAVING MODEL AS PICKLE (.mpck): %s" % model_pck), level=1)
    model.write_pickle(model_pck)
