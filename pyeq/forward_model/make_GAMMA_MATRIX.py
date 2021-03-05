###########################################################################
# GAMMA_MATRIX gives information about obs available in model_time_step
# GAMMA_MATRIX[ idx_gps_site , model_time_step ] = 0 if the first observation available is after model_time_step
# GAMMA_MATRIX[ idx_gps_site , model_time_step ] = model_time_step_duration if the first observation is before model_time_step
# GAMMA_MATRIX[ idx_gps_site , model_time_step ] = duration of the site in the current model_time_step duration if the first observation is available during model_time_step, and alpha is the fraction of the current time step

def make_GAMMA_MATRIX( T_OBS , np_obs_date_s , np_model_date_s ):

    import numpy as np
    
    np_model_step_duration_in_days = np.diff( np_model_date_s ) / 86400.
    
    GAMMA_MATRIX = np.zeros( ( T_OBS.shape[1], np_model_date_s.shape[0]-1 ) )
    
    
    # loop on site
    
    for idx_site in np.arange( T_OBS.shape[1] ):
    
        # get the first not nan observation for the current site
        idx_first_obs = np.where( np.isfinite(T_OBS[:,idx_site,0])  )[0][0]
        
        GAMMA_MATRIX[ idx_site , : ] = np.heaviside( ( np_model_date_s[1:] - np_obs_date_s[ idx_first_obs ] ) , 0 ) 
        
        # get the first model time step for site idx_site 
        idx_ts = np.where( GAMMA_MATRIX[ idx_site , : ] == 1 )[0][0] 

        # rescale by np_model_step_duration_in_days
        GAMMA_MATRIX[ idx_site , : ] = GAMMA_MATRIX[ idx_site , : ] * np_model_step_duration_in_days

        
        # renormalize by the fraction of the time step x time step duration
        frac = 1. - ( np_obs_date_s[ idx_first_obs ] - np_model_date_s[ idx_ts ] ) / (  np_model_date_s[ idx_ts +1 ] -  np_model_date_s[ idx_ts ] ) 
        time_step_duration = ( np_model_date_s[ idx_ts +1 ] -  np_model_date_s[ idx_ts ] ) / 86400.
        
        # change 19/04/2020
        #GAMMA_MATRIX[ idx_site , idx_ts ] = frac
        GAMMA_MATRIX[ idx_site , idx_ts ] = frac * time_step_duration
        
        #if model_step < np_model_date_s.shape[0] - 1: 
        #    GAMMA_MATRIX[ idx_site , model_step+1: ] = 1.
        
    return GAMMA_MATRIX
