def check_date_obs_vs_model( model ):


    import pyacs.lib.astrotime as at
    import numpy as np
    import textwrap


    ###########################################################################
    # MAKING OBSERVATION AND MODEL DATE CONSISTENT
    ###########################################################################
    
    # we will need np_obs_dates including the observation dates in seconds since np_model_dates_s[0]
    
    print("-- Checking that model and observation dates are consistent")
    print("-- number of model dates requested by user: %d" % ( model.np_model_date_s.shape[0] ) )
    print("-- start     model date  requested by user: %s" % ( at.seconds2datetime( model.np_model_date_s[0] ).isoformat(' ')) )
    print("-- end       model date  requested by user: %s" % ( at.seconds2datetime( model.np_model_date_s[-1] ).isoformat(' ')) )

    # check compatibility of obs and model dates
    
    print("-- Checking model/obs date consistency")
    if model.np_model_date_s[0] < model.np_obs_date_s[0]:
        print("!!!WARNING: requested model date starts before data. Changing model start date to the first available observation date.")
        model.np_model_date_s = np.append( model.np_model_date_s , model.np_obs_date_s[0] )
    
    if model.np_model_date_s[-1] > model.np_obs_date_s[-1]:
        print("!!!WARNING: requested model date (%d) ends after data end (%d)" % ( model.np_model_date_s[-1] , model.np_obs_date_s[-1] ))
        print("!!! Changing model end date to the last available observation date.")
        model.np_model_date_s = np.append( model.np_model_date_s , model.np_obs_date_s[-1] )

    model.np_model_date_s = np.unique( np.sort( model.np_model_date_s ))
    
    lindex = np.where( ( model.np_model_date_s>= model.np_obs_date_s[0]) & ( model.np_model_date_s<= model.np_obs_date_s[-1]) )
    model.np_model_date_s = model.np_model_date_s[lindex]
    
    # extract the model period for the time series
    
    print("-- Extracting time series for the user defined model period"  )
    lindex = np.where( (model.np_obs_date_s >= model.np_model_date_s[0]) &  (model.np_obs_date_s <= model.np_model_date_s[-1]) )
    model.t_obs = model.t_obs_raw[lindex]
    model.np_obs_date_s = model.np_obs_date_s[ lindex ]
    print("-- Extracted observation tensor for modeling period: " , model.t_obs.shape )
    print("-- Number of observation dates: %d" % (model.np_obs_date_s.shape[0]) )
    print("-- number of observation dates used for modeling: %d " % model.t_obs.shape[0])
    
    # remove GPS site with no data
    print("-- removing GPS sites with less than two data during the model period")
    np_index_gps_no_data = np.array([] , dtype=int )
    
    for i in np.arange( model.np_gps_site.shape[0] ):
        code = model.np_gps_site[i]
        ndata = np.count_nonzero(~np.isnan( model.t_obs[:,i][:,0]))
        if  ndata < 2:
            np_index_gps_no_data = np.append( np_index_gps_no_data , i ) 
         
    print("-- GPS site with less than two data within the model period")
    print("%s" % ( textwrap.fill(" ".join( model.np_gps_site[ np_index_gps_no_data ].tolist() ) , width=80 ) ))
    model.warning = model.warning + ("%s" % ( textwrap.fill(" ".join( model.np_gps_site[ np_index_gps_no_data ].tolist() ) , width=80 ) ))
    model.np_gps_site = np.delete( model.np_gps_site, np_index_gps_no_data )
    np_common_index_obs_tensor = np.delete( np.arange( model.t_obs.shape[1] ), np_index_gps_no_data )
    np_common_index_elastic_tensor = np.delete( np.arange( model.green.shape[0] ), np_index_gps_no_data )

    # making extraction    
    model.t_obs = model.t_obs[ : , np_common_index_obs_tensor, :]
    model.green = model.green[np_common_index_elastic_tensor]
    
    model.n_site_no_data = np_index_gps_no_data.shape[0]
    print("-- Observation tensor new shape: " , model.t_obs.shape )
    print("-- Elastic tensor new shape: " , model.green.shape )
    
    print("-- final list of GPS sites used in inversion: %d" % ( model.np_gps_site.shape[0]) )
    print("%s" % ( textwrap.fill(" ".join( model.np_gps_site.tolist() ) , width=80 ) ))

    return model