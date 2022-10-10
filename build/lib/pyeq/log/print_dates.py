def print_dates( model ):
    """
    Print model and observation dates
    """
    
    # import
    import numpy as np
    import pyacs.lib.astrotime as at

    import pyeq.message.message as MESSAGE
    import pyeq.message.verbose_message as VERBOSE
    import pyeq.message.error as ERROR
    import pyeq.message.warning as WARNING
    import pyeq.message.debug_message as DEBUG

    ###################################################################
    # DATES
    # as input we have two 1D numpy array of dates
    # np_model_date_s: dates at which the model cumulative slip is estimated
    # np_obs_date_s  : dates where at least one observation is available
    # both np_date_model_s & np_date_obs_s are in integer seconds since 1980/1/1 at 00:00:00
    # for printing the results, both dates are converted:
    # delta_d: float being the time in decimal days since the first model date
    # isoformat: is a string YYYY/MM/DAY HR:MN:SC
    # slip rate, stf (moment rate), inc_stf (incremental moment every model time step) correspond to the value during one model time step
    # their date is therefore the middle of the model time step.
    # cumulative slip and and cstf are the estimated model values at the model dates
    ###################################################################
    
    # input dates in integer seconds since 1980/1/1 00:00:00

    np_model_date_s = model.np_model_date_s
    np_obs_date_s = model.np_obs_date_s

    # model step duration in seconds
    np_model_step_duration_s =np.diff( np_model_date_s)
    model.np_model_step_duration_s = np_model_step_duration_s
    # middle of the model dates in seconds
    np_mid_model_date_s = ( np_model_date_s[:-1] +  np.diff( np_model_date_s)/2. ).astype( int )
    model.np_mid_model_date_s = np_mid_model_date_s
    # middle of the model dates in datetime
    np_mid_model_datetime = at.seconds2datetime( np_mid_model_date_s )
    model.np_mid_model_datetime = np_mid_model_datetime
    # middle of the model dates in isoformat
    np_mid_model_isoformat = np.array([x.isoformat(' ') for x in np_mid_model_datetime ])
    model.np_mid_model_isoformat = np_mid_model_isoformat
    # middle of the model dates in delta_days since model first date
    np_mid_model_delta_d = ( np_mid_model_date_s - np_model_date_s[0] ) / 86400.
    model.np_mid_model_delta_d = np_mid_model_delta_d
  
    # model step duration in decimal days
    np_model_step_duration_days = np_model_step_duration_s / 86400.
    model.np_model_step_duration_days = np_model_step_duration_days
    
    # model dates in datetime
    np_model_datetime = at.seconds2datetime( np_model_date_s ) 
    model.np_model_datetime = np_model_datetime
    
    # model dates in isoformat
    np_model_date_isoformat = np.array([x.isoformat(' ') for x in np_model_datetime ])
    model.np_model_date_isoformat = np_model_date_isoformat
    # model dates in delta_days since model first date
    np_model_delta_d = ( np_model_date_s - np_model_date_s[0] ) / 86400.
    model.np_model_delta_d = np_model_delta_d
    # obs dates in datetime
    np_obs_datetime = at.seconds2datetime( np_obs_date_s ) 
    model.np_obs_datetime = np_obs_datetime

    # obs dates in isoformat
    np_obs_date_isoformat = np.array([x.isoformat(' ') for x in np_obs_datetime ])
    model.np_obs_date_isoformat = np_obs_date_isoformat

    # obs dates in delta_days since model first date
    np_obs_delta_d = ( np_obs_date_s - np_model_date_s[0] ) / 86400.
    model.np_obs_delta_d = np_obs_delta_d

    syear, sdoy, _ut  = at.datetime2dayno( np_model_datetime[0] )
    eyear, edoy, _ut  = at.datetime2dayno( np_model_datetime[-1] )
     
    median_time_step_days   = np.median( np_model_step_duration_days )
    shortest_time_step_days = np.min( np_model_step_duration_days )
    longest_time_step_days  = np.max( np_model_step_duration_days )

    model.median_time_step_days   = median_time_step_days
    model.shortest_time_step_days = shortest_time_step_days
    model.longest_time_step_days  = longest_time_step_days


    ###########################################################################
    # SAVE MODEL DATES IN DIR: INFO
    ###########################################################################

    
    wf = open( model.odir+'/info/model_dates.dat' , 'w' )
    VERBOSE("saving model dates file in %s" % (model.odir+'/info/model_dates.dat') )
    
    wf.write("#step  start_date            end_date             delta_d delta_t0 start_doy   end_doy   start_decyear     end_decyear\n")
    for i in np.arange( np_model_date_s.shape[0]-1 ):
        
        sdatetime = np_model_datetime[i]
        edatetime = np_model_datetime[i+1]
        iso_sdate = np_model_date_isoformat[i]
        iso_edate = np_model_date_isoformat[i+1]
        delta_d   = ( np_model_date_s[i+1] - np_model_date_s[i] ) / 86400.
        delta_t0  = ( np_model_date_s[i+1] - np_model_date_s[0] ) / 86400.
        syear, sdoy, _ut  = at.datetime2dayno( sdatetime )
        eyear, edoy, _ut  = at.datetime2dayno( edatetime )

        sdecyear = at.datetime2decyear( sdatetime )
        edecyear = at.datetime2decyear( edatetime )

        wf.write("%04d  %s -> %s %8.3lf %8.3lf %04d %03d -> %04d %03d %15.10lf -> %15.10lf \n" % ( i, iso_sdate, iso_edate, delta_d, delta_t0, syear, sdoy, eyear, edoy, sdecyear, edecyear ))

    wf.close()



    
    ###########################################################################
    # SAVE OBSERVATION DATES IN DIR: INFO
    ###########################################################################


    wf = open( model.odir+'/info/obs_dates.dat' , 'w' )
    VERBOSE("saving observation dates file in %s" % (model.odir+'/info/obs_dates.dat') )
    
    wf.write("#obs  date                 delta_t0  year doy      decyear\n")
    for i in np.arange( np_obs_date_s.shape[0] ):
        
        sdatetime = np_obs_datetime[i]
        iso_sdate = np_obs_date_isoformat[i]
        delta_t0  = np_obs_delta_d[i]
        syear, sdoy, _ut  = at.datetime2dayno( sdatetime )
        sdecyear = at.datetime2decyear( sdatetime )

        wf.write("%04d  %s  %8.3lf  %04d %03d %15.10lf\n" % ( i, iso_sdate, delta_t0, syear, sdoy, sdecyear ))

    wf.close()

    return model


    