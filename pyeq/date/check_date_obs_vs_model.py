def check_date_obs_vs_model( model ):


    import pyacs.lib.astrotime as at
    import numpy as np
    import textwrap

    import pyeq.message.message as MESSAGE
    import pyeq.message.verbose_message as VERBOSE
    import pyeq.message.error as ERROR
    import pyeq.message.warning as WARNING
    import pyeq.message.debug_message as DEBUG

    ###########################################################################
    # MAKING OBSERVATION AND MODEL DATE CONSISTENT
    ###########################################################################
    
    # we will need np_obs_dates including the observation dates in seconds since np_model_dates_s[0]
    
    DEBUG("Checking that model and observation dates are consistent")
    DEBUG("number of model dates requested by user: %d" % ( model.np_model_date_s.shape[0] ) )
    DEBUG("start     model date  requested by user: %s" % ( at.seconds2datetime( model.np_model_date_s[0] ).isoformat(' ')) )
    DEBUG("end       model date  requested by user: %s" % ( at.seconds2datetime( model.np_model_date_s[-1] ).isoformat(' ')) )

    # check compatibility of obs and model dates
    
    DEBUG("Checking model/obs date consistency")
    if model.np_model_date_s[0] < model.np_obs_date_s[0]:
        WARNING("requested model date starts before data. Changing model start date to the first available observation date.")
        model.np_model_date_s = np.append( model.np_model_date_s , model.np_obs_date_s[0] )
    
    if model.np_model_date_s[-1] > model.np_obs_date_s[-1]:
        WARNING("requested model date (%d) ends after data end (%d)" % ( model.np_model_date_s[-1] , model.np_obs_date_s[-1] ))
        WARNING("Changing model end date to the last available observation date.")
        model.np_model_date_s = np.append( model.np_model_date_s , model.np_obs_date_s[-1] )

    model.np_model_date_s = np.unique( np.sort( model.np_model_date_s ))
    
    lindex = np.where( ( model.np_model_date_s>= model.np_obs_date_s[0]) & ( model.np_model_date_s<= model.np_obs_date_s[-1]) )
    model.np_model_date_s = model.np_model_date_s[lindex]
    
    # extract the model period for the time series
    
    DEBUG("Extracting time series for the user defined model period"  )
    lindex = np.where( (model.np_obs_date_s >= model.np_model_date_s[0]) &  (model.np_obs_date_s <= model.np_model_date_s[-1]) )
    model.t_obs = model.t_obs_raw[lindex]
    model.np_obs_date_s = model.np_obs_date_s[ lindex ]
    DEBUG("Extracted observation tensor for modeling period: %s" % (model.t_obs.shape,) )
    DEBUG("Number of observation dates: %d" % (model.np_obs_date_s.shape[0]) )
    DEBUG("number of observation dates used for modeling: %d " % model.t_obs.shape[0])
    
    # remove GPS site with no data
    DEBUG("removing GPS sites with less than two data during the model period")
    np_index_gps_no_data = np.array([] , dtype=int )
    
    for i in np.arange( model.np_gps_site.shape[0] ):
        code = model.np_gps_site[i]
        ndata = np.count_nonzero(~np.isnan( model.t_obs[:,i][:,0]))
        if  ndata < 2:
            np_index_gps_no_data = np.append( np_index_gps_no_data , i ) 

    if model.np_gps_site[ np_index_gps_no_data ].shape[0] > 0:
        WARNING("some GPS site have less than two data within the model period. Check info/warning.dat")
        DEBUG("%s" % ( textwrap.fill(" ".join( model.np_gps_site[ np_index_gps_no_data ].tolist() ) , width=80 ) ))
        model.warning += ("%s" % ( textwrap.fill(" ".join( model.np_gps_site[ np_index_gps_no_data ].tolist() ) , width=80 ) ))

    model.np_gps_site = np.delete( model.np_gps_site, np_index_gps_no_data )
    np_common_index_obs_tensor = np.delete( np.arange( model.t_obs.shape[1] ), np_index_gps_no_data )
    np_common_index_elastic_tensor = np.delete( np.arange( model.green.shape[0] ), np_index_gps_no_data )

    # making extraction    
    model.t_obs = model.t_obs[ : , np_common_index_obs_tensor, :]
    model.green = model.green[np_common_index_elastic_tensor]
    
    model.n_site_no_data = np_index_gps_no_data.shape[0]
    DEBUG("Observation tensor new shape: %s" % (model.t_obs.shape,) )
    DEBUG("Elastic tensor new shape: %s" % (model.green.shape,) )
    
    VERBOSE("final list of GPS sites used in inversion: %d" % ( model.np_gps_site.shape[0]) )
    DEBUG("%s" % ( textwrap.fill(" ".join( model.np_gps_site.tolist() ) , width=80 ) ))

    return model