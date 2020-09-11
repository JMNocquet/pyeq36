"""
Routines to build the linear system
This is version 3 handling data not starting a the model start date
Fixed rake only
"""


def build3( model ):

    print("-- BUILDING THE NORMAL SYSTEM USING PYEQ.LIB.BUILD3 --")

    ### TO BE DONE
    ### renormalize Gk by Cd-1/2
    ### RAKE handling - check whether GREEN has been modified before build2 call
    ### ADD ORIGIN TIME CONSTANT AS ESTIMATED PARAMETER
    ### PARAMETER TRANSLATION FOR CONJUGATE RAKE AND/OR ORIGIN TIME CONSTANT
    
    n_constant = 0
    RAKE = False

    ###########################################################################
    # import
    ###########################################################################
    
    import numpy as np
    import resource
    
    import pyeq.lib.elastic_tensor.shrink
    import pyeq.lib.obs_tensor.dk_from_obs_tensor
    import pyeq.lib.forward_model.C_B_at_k
    import pyeq.lib.obs_tensor
    import pyacs.lib.glinalg
    import pyacs.lib.astrotime as at

    ###########################################################################
    # COMPONENTS & RAKE
    ###########################################################################
    
    # components
    if model.up:
        np_component = np.array( [0,1,2] )
    else:
        np_component = np.array( [0,1] )
    
    if model.debug:
        print("-- np_component:" , np_component )
    
    # rake
    if RAKE:
        np_rake = np.array( [0,1] )
    else:
        np_rake = np.array( [0] )

    model.np_rake = np_rake

    nrake = np_rake.shape[0] 
    
    # debug
    if model.debug:
        print("-- np_rake:" , np_rake )


    ###########################################################################
    # PRINT SOME INFO
    ###########################################################################

    print('-- list of observation dates')
    for i in np.arange( model.t_obs.shape[0] ):
        print("%s" % ( at.seconds2datetime( model.np_obs_date_s[i] ).isoformat(' ')   ) )
        
    print('-- list of model dates')
    for i in np.arange( model.np_model_date_s.shape[0] ):
        print("%s" % ( at.seconds2datetime( model.np_model_date_s[i] ).isoformat(' ')   ) )

    
    ###########################################################################
    # GAMMA_MATRIX gives information about obs available in model_time_step
    # GAMMA_MATRIX[ idx_gps_site , model_time_step ] = 0 if the first observation available is after model_time_step
    # GAMMA_MATRIX[ idx_gps_site , model_time_step ] = 1 if the first observation is before model_time_step
    # GAMMA_MATRIX[ idx_gps_site , model_time_step ] = alpha if the first observation is available during model_time_step, and alpha is the fraction of the current time step

    def make_GAMMA_MATRIX( T_OBS , np_obs_date_s , np_model_date_s ):
        
        GAMMA_MATRIX = np.zeros( ( T_OBS.shape[1], np_model_date_s.shape[0]-1 ) )
        
        # loop on site
        
        for idx_site in np.arange( T_OBS.shape[1] ):
        
            # get the first not nan observation for the current site
            idx_first_obs = np.where( np.isfinite(T_OBS[:,idx_site,0])  )[0][0]
            
            GAMMA_MATRIX[ idx_site , : ] = np.heaviside( ( np_model_date_s[1:] - np_obs_date_s[ idx_first_obs ] ) , 0 )
            
            # get the first model time step for site idx_site 
            idx_ts = np.where( GAMMA_MATRIX[ idx_site , : ] == 1 )[0][0] 
            
            # renormalize by the fraction of the time step
            
            frac = 1. - ( np_obs_date_s[ idx_first_obs ] - np_model_date_s[ idx_ts ] ) / (  np_model_date_s[ idx_ts +1 ] -  np_model_date_s[ idx_ts ] ) 
            
            GAMMA_MATRIX[ idx_site , idx_ts ] = frac
            
            #if model_step < np_model_date_s.shape[0] - 1: 
            #    GAMMA_MATRIX[ idx_site , model_step+1: ] = 1.
            
        return GAMMA_MATRIX

    ###########################################################################

    for idx_site in np.arange( model.t_obs.shape[1] ):
    
        # get the first not nan observation for the current site
        idx_first_obs = np.where( np.isfinite(model.t_obs[:,idx_site,0])  )[0][0]
        date_s = model.np_obs_date_s[ idx_first_obs ]
        if model.verbose:
            print("-- site %s starts at idx date %04d  %s "  % (model.np_gps_site[ idx_site ] , idx_first_obs , at.seconds2datetime(date_s).isoformat(' ')  ) )
        model.warning = model.warning + \
            ("[WARNING] site %s starts at idx date %04d  %s \n"  % (model.np_gps_site[ idx_site ] , idx_first_obs , at.seconds2datetime(date_s).isoformat(' ')  ) )
    
    GAMMA_MATRIX = make_GAMMA_MATRIX(model.t_obs, model.np_obs_date_s, model.np_model_date_s)

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


    ###########################################################################
    # ATPA ATBP and various variables initialization
    ###########################################################################

    n_model_time_step = model.np_model_date_s.size - 1
    nfaults = model.green.shape[1]
    
    print("-- Initializing ATPA & ATPB")
    
    # ATPA
    
    memory_usage = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024 / 1024
    print("-- memory usage before ATPA initialization in Gb Linux/Mac OS X : %.1lf / %.1lf" %  (memory_usage,memory_usage/1024) )

    n_atpa = nfaults * nrake * n_model_time_step + n_constant
    ATPA = np.zeros( (n_atpa,n_atpa) )

    memory_usage = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024 / 1024
    print("-- memory usage after ATPA initialization in Gb Linux/Mac OS X : %.1lf / %.1lf" %  (memory_usage,memory_usage/1024) )
    
    # ATPB
    
    ATPB = np.zeros( (n_atpa) )

    memory_usage = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024 / 1024
    print("-- memory usage after ATPB initialization in Gb Linux/Mac OS X : %.1lf / %.1lf" %  (memory_usage,memory_usage/1024) )

    # COLUMN OF ATPA CORRESPONDING TO A TIME STEP

    # number of columns for C
    ncc = nfaults * nrake
    print('-- nrake' , nrake)
    print('-- nfaults' , nfaults)
    print("-- Initializing the individual model time step matrices and RHS vector")
    C = np.zeros( ( ncc * n_model_time_step , ncc ) )
    # RHS sub-vector
    Bi = np.zeros( ( ncc * ( n_model_time_step ) ) )

    # GAMMA MODEL STEP DURATION VECTOR
    
    print("-- Initializing the vector of time duration")
    GAMMA = np.diff( model.np_model_date_s , 1 )
    
    # G_PER_SITE
    if model.debug:
        print("-- Initializing GT_PER_SITE")
    
    G_PER_SITE = model.green[ : ][:, :][:,:,np_component][:,:,:,np_rake]
    G_PER_SITE = G_PER_SITE.reshape( G_PER_SITE.shape[0] , G_PER_SITE.shape[1]* G_PER_SITE.shape[3] , G_PER_SITE.shape[2] )
    G_PER_SITE = np.moveaxis(G_PER_SITE,[0,1,2],[0,2,1]) 

    if model.debug:
        print("-- G_PER_SITE " , G_PER_SITE.shape )
        print("GAMMA_MATRIX " , GAMMA_MATRIX.shape )
    
    ###########################################################################
    # loop on model time steps
    ###########################################################################

    
    # loop is made reverse in time
    for model_time_step in np.arange( n_model_time_step-1, -1, -1 ):

        ssdate_s = model.np_model_date_s[ model_time_step ]
        eedate_s = model.np_model_date_s[ model_time_step + 1 ]


        print('######################################################################')  
        #if model.verbose:
        print("-- Building equations for model time step: %04d %s -> %s" % \
              ( model_time_step , at.seconds2datetime(ssdate_s).isoformat(' '), at.seconds2datetime(eedate_s).isoformat(' ') ) )
     
        # Get the list of the k-index of observation time within the current model time step
        
        lk = H_i_k[ model_time_step ] 
        
        
        #######################################################################
        # loop on k observation time within the current model time step
        #######################################################################
        
        
        
        for k in lk:
            
            # get tk time for observation k
            
            tk = model.np_obs_date_s[ k ]
            
            print("-- model time step: %04d observation time step: %04d %s" % (model_time_step , k, at.seconds2datetime(tk).isoformat(' ') ) )
            
            # Get dk,sdk observation/sigma vector at date tk for available sites with the chosen components
            if model.debug:
                print("-- model time step: %04d observation time step: %04d get normalized observation using pyeq.lib.obs_tensor.dk_from_obs_tensor"  % (model_time_step , k ) )
            
            dk  = model.t_obs[k ,:,np_component].T
            sdk = model.t_obs[k ,:,np_component+3].T

            # normalize the Green matrix
            #print( GAMMA_MATRIX )
            Ctmp, Btmp = pyeq.lib.forward_model.C_B_at_k( k, model.np_obs_date_s, model.np_model_date_s, model_time_step, GAMMA_MATRIX, G_PER_SITE, dk, sdk)
            
                    
            # Update C (See "Algorithm to build the observation normal system in paper")
            if model.debug:
                print("-- model time step: %04d observation time step: %04d Individual column matrix for ATPA"  % (model_time_step , k ) )
            C = C + Ctmp
            
            # Update Bk
            if model.debug:
                print("-- model time step: %04d observation time step: %04d Individual vector for ATPB"  % (model_time_step , k ) )
            Bi = Bi + Btmp
            
            
        #######################################################################
        # end loop on k
        #######################################################################

        # C is ready to be inserted into ATPA
        
        # inserting C as a column
        if model.debug:
            print("-- model time step: %04d Inserting the normal observation column in ATPA"  % (model_time_step ) )
        ATPA[ :C.shape[0] , model_time_step*ncc:(model_time_step+1)*ncc ] = C


        # inserting C.T as a row
        if model.debug:
            print("-- model time step: %04d Inserting the normal observation row in ATPA"  % (model_time_step ) )
        ATPA[ model_time_step*ncc:(model_time_step+1)*ncc , :C.shape[0] ] = C.T

        # Updating ATPB !!!! THINK MORE
        if model.debug:
            print("-- model time step: %04d Updating ATPB"  % (model_time_step ) )
        
        ATPB[ model_time_step*ncc:(model_time_step+1)*ncc ] = Bi[ model_time_step*ncc:(model_time_step+1)*ncc ] 


        #######################################################################
        # end of the current model time step
        #######################################################################
        
        # Remove the last rows of C corresponding to the current time step
        if model.debug:
            print("-- model time step: %04d Resizing individual column matrix for next model time step"  % (model_time_step ) )
        
        C = C[:-ncc,:]
        Bi= Bi[ :-ncc ]
        
    return ATPA , ATPB


