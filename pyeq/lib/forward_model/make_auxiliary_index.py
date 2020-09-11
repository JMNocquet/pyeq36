
def make_auxiliary_index( model ):
    """
    Builds a dictionary of index for observation time included at each model_time_step
    """
    
    import numpy as np
    import pyacs.lib.astrotime as at

    ###########################################################################
    # AUXILIARY INDEX
    ###########################################################################

    if model.verbose:
        print("-- Building H_i_k  : dictionary providing the list of indexes k for observation time for each model time step i")
    # H_i_k  : dictionary providing the list of indexes k for observation time for each model time step i
    H_i_k = {}

    
    for model_time_step in np.arange( model.np_model_date_s.shape[0]-2 , -1, step = -1 ):
        sdate_s = model.np_model_date_s[ model_time_step ]
        edate_s = model.np_model_date_s[ model_time_step + 1 ]
        if model.verbose:
            print("--- model_time_step: %04d %s -> %s " % ( model_time_step, at.seconds2datetime(sdate_s).isoformat(' '), at.seconds2datetime(edate_s).isoformat(' ') ) )
        H_i_k[ model_time_step ] = []
        
        for k in np.arange(model.np_obs_date_s.shape[0]):
            if ( ( model.np_obs_date_s[k] > sdate_s ) and ( model.np_obs_date_s[k] <= edate_s ) ):
                #print( model.np_obs_date_s[k] )
                H_i_k[ model_time_step ].append( k )
                 
        H_i_k[ model_time_step ] = sorted( H_i_k[ model_time_step ] )
        
    if model.verbose:
        print("--- list of observation index per model time step")
        for key, value in H_i_k.items():
            
            print("--- model time step %d: %s " % (key,value) )

    return H_i_k