def C_B_at_k( \
            k, 
            np_date_obs_s, \
            np_date_model_s, \
            model_time_step, \
            GAMMA_MATRIX, \
            G_PER_SITE, \
            dk, \
            sdk, \
            verbose = False, \
            ):
    """
    Calculates the C columns corresponding to observation at index k obs time step tk and the associated RHS B
    """

    # import
    
    import numpy as np
    from progress.bar import Bar

    # deep copy because we do not want GAMMA_MATRIX to be changed outside this routine
    
    GGAMMA_MATRIX = np.copy( GAMMA_MATRIX )
    
    # step duration in days
    
    # model step duration in seconds
    np_model_step_duration_day = np.diff( np_date_model_s) / 86400.
    
    # dimension of model parameters per time step
    n_model = G_PER_SITE.shape[2]
    
    # number of time step involved in building C and B. It is equal to the index of the current time step for obs date k + 1
    n_model_time_step_obs_k = model_time_step + 1
    
    # initialize
    
    C = np.zeros( ( n_model_time_step_obs_k * n_model , n_model ) )
    B = np.zeros( ( n_model_time_step_obs_k * n_model ) )

    # delta for site with no observation at date k
    
    delta = np.zeros( ( G_PER_SITE.shape[0] ) ) + 1.
    lidx_missing = np.where( np.isnan( dk[:,0]) )[0]
    delta[ lidx_missing ] = 0.
    
    # Normalize observation

    W_Sdk = 1./ np.nan_to_num(sdk , nan=1.)
    DK    = np.nan_to_num( dk, nan=0. ) * W_Sdk

    # Update the GGAMMA_MATRIX for the current time step model_time_step
    frac_time_step = ( np_date_obs_s[k] - np_date_model_s[ model_time_step ]  ) / ( np_date_model_s[ model_time_step+1 ] - np_date_model_s[ model_time_step ] )

    GGAMMA_MATRIX[ : , model_time_step ] =  np.heaviside( GGAMMA_MATRIX[ : , model_time_step ] - (1. - frac_time_step) , 0  )  \
                                           * ( GGAMMA_MATRIX[ : , model_time_step ] - (1. - frac_time_step) )  \
    

    # loop on model time step

    # print progression bar
    bar = Bar('', max= n_model_time_step_obs_k , suffix = '%(percent).1f%% - %(eta)ds')

    for idx_ts in np.arange( n_model_time_step_obs_k ):
        bar.next()

        # duration of current time step in day
        step_duration_day = np_model_step_duration_day[ idx_ts ]

        # loop on sites
        
        for idx_site in np.arange( G_PER_SITE.shape[0] ):
            green_factor = delta[ idx_site ]  * GGAMMA_MATRIX[ idx_site , idx_ts ] * step_duration_day
            
            if ( green_factor != 0):
                
                # This line does the following
                # If the site has no obs a date k, then delta is zero and the Green tensor is 0
                # GAMMA_MATRIX tells whether the multiplication for the preceding model time step
                
                
                Gk = (G_PER_SITE[ idx_site ].T * W_Sdk[ idx_site ]).T
                #Gk = np.ascontiguousarray( (G_PER_SITE[ idx_site ].T * W_Sdk[ idx_site ]).T )
                
                # RHS
                # BUG CORRECTED 08/04/2020 - AFFECT INVERSIONS WITH STEP DURATION != 1 DAY
                # Gk should be multiplied by green_factor 
                B[ idx_ts * n_model : (idx_ts+1) * n_model ] = B[ idx_ts * n_model :  (idx_ts+1) * n_model ] \
                                                                + np.dot( Gk.T , DK[ idx_site ].reshape(-1,1) ).flatten() \
                                                                * green_factor
            
                # C
                
                C[ idx_ts * n_model : (idx_ts+1) * n_model , : ] = C[ idx_ts * n_model : (idx_ts+1) * n_model , : ] \
                                                            + green_factor * np.dot( Gk.T , Gk )

            else:
                if verbose:
                    print("!!! No observation to be added for site %d at model step %d during model step handling %d " % (idx_site,idx_ts,model_time_step) )
            
    bar.finish()
    
    return C , B
                                                            
                                                            
                                                            
                                                            
                                                            
                                                            
                                                            
                                                            
                                                            
                                                            
                                                            
                                                             