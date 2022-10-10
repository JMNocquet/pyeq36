def C_B_at_k_build4( \
            k, \
            nstep, \
            np_site_name, \
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

    # Update the GGAMMA_MATRIX for the current time step model_time_step
    frac_time_step = ( np_date_obs_s[k] - np_date_model_s[ model_time_step ]  ) / ( np_date_model_s[ model_time_step+1 ] - np_date_model_s[ model_time_step ] )

    GGAMMA_MATRIX[ : , model_time_step ] =  np.heaviside( GGAMMA_MATRIX[ : , model_time_step ] - (1. - frac_time_step) , 0  )  \
                                           * ( GGAMMA_MATRIX[ : , model_time_step ] - (1. - frac_time_step) )  \
    
    ###########################################################################
    # INITIALIZE
    ###########################################################################
    
    BB = np.array([])
    AA = None

    # delta for site with no observation at date k
    
    delta = np.zeros( ( G_PER_SITE.shape[0] ) ) + 1.
    lidx_missing = np.where( np.isnan( dk[:,0]) )[0]
    delta[ lidx_missing ] = 0.
    
    # Normalize observation

    W_Sdk = 1./ np.nan_to_num(sdk , nan=1.)
    DK    = np.nan_to_num( dk, nan=0. ) * W_Sdk
    CDK   = ( np.nan_to_num( dk*0.+1., nan=0. ) * W_Sdk ).flatten()

    # Column matrix for offset at the origin time - ATPA LHS
    CCST = np.zeros( ( CDK.shape[0]  + ( nstep * G_PER_SITE.shape[2] ) , CDK.shape[0] ) ) 

    # Vector of ATPB RHS corresponding to offsets
    ATPB_OFFSET = np.zeros( CDK.shape[0] )

    ###########################################################################
    # MAIN LOOP ON THE MODEL TIME STEP BEFORE CURRENT OBSERVATION TIME
    ###########################################################################

    # print progression bar
    bar = Bar('', max= n_model_time_step_obs_k , suffix = '%(percent).1f%% - %(elapsed)ds')

    # loop on model the model time steps before the current observation
    # this produces a column-matrix to be added to N
    # the length of the columns involves all model time steps before the current observation time
    
    for idx_ts in np.arange( n_model_time_step_obs_k ):
        #print('sub-step #',idx_ts)
        bar.next()
        
        # duration of current time step in day
        step_duration_day = np_model_step_duration_day[ idx_ts ]
        # The Green's function are multiplied by the vector green_factor whose element tell for each site
        # delta[site] = 0: if there is no observation for the current date to be included to the normal system
        # GGAMMA_MATRIX[ site , idx_ts ] = 0 if the site has observation starting later than the end date of the idx_ts model step
        #                                = 1 if the site has observation starting before the start date of the idx_ts model step
        #                                = frac if the current observation is during the current idx_ts model time step. frac is the fraction 
        #                                  of the model step for the current date
         
        green_factor = delta[ : ]  * GGAMMA_MATRIX[ : , idx_ts ] * step_duration_day 
        # for test
        #green_factor = delta[ : ]  * GGAMMA_MATRIX[ : , idx_ts ] 
        
        #print( 'green factor \n' , green_factor )
        
        # GK has shape (ngps,ncomponent,nfault) allowing multiplication by green_factor (green_factor includes the value different for each site)
        #print('G_PER_SITE\n',G_PER_SITE)
        #print('GGAMMA_MATRIX[ : , idx_ts ]' , GGAMMA_MATRIX[ : , idx_ts ] )
        GK = (G_PER_SITE.T * green_factor ).T 
        #print('GK\n',GK)
        # GK transform to 2D-array
        GK = GK.reshape( GK.shape[0]*GK.shape[1] , GK.shape[2] ) 

        # renormalize by Cd-1/2
        GK = ( GK.T * CDK ).T

        # Compute the column matrix for constants 
        # get the index where there is no data involved
        
        try:
            CST = np.vstack((CST, GK.T  ))
        except:
            CST = GK.T
        
        # AA is obtained by stacking for all idx_ts
        # AA is initialized as None 
        try:
            AA = np.vstack(( AA, np.dot( GK.T , GK )  ))
        except:
            AA = np.dot( GK.T , GK )
        
        # BB is obtained by stacking for all idx_ts
        BB = np.hstack( (BB, np.dot( GK.T , DK.reshape(-1,1).flatten() )  ))
        
    bar.finish()


    # DEAL WITH OFFSETS
    
    # OFFSET CONSTANT ATPA COLUMN
    
    # SIGMA GK.T part
    CCST[ :CST.shape[0] , : ] = CST
    # CCST has an extra block of Identity * Cd-1 at its bottom corresponding to its own observation (see paper)
    LAST_CST = ( (W_Sdk**2).T * delta ).T.flatten()
    CCST[-LAST_CST.shape[0]: , : ] = np.diag( LAST_CST )

    # OFFSET ATPB SEGMENT

    ATPB_OFFSET = DK.flatten()
    
    return AA , BB , CCST , ATPB_OFFSET
                                                            
                                                            
                                                            
                                                            
                                                            
                                                            
                                                            
                                                            
                                                            
                                                            
                                                            
                                                             