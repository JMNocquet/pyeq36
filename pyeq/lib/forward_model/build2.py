"""
Routines to build the linear system
build 2 allows variable model time step but not site starting after the first observation date - use build 3
"""

###############################################################################
def build2( model ):
###############################################################################
    
    RAKE = False
    n_constant = 0
    

    ###########################################################################
    # import
    ###########################################################################
    
    import numpy as np
    import resource
    import textwrap
    
    # pyacs
    import pyacs.lib.glinalg
    import pyacs.lib.astrotime as at
    # pyeq
    import pyeq.lib.elastic_tensor
    import pyeq.lib.obs_tensor
    import pyeq.lib.obs_tensor.dk_from_obs_tensor
    

    terminal_width = 80
    
    ###########################################################################
    # CHECK WETHER ANY OBSERVATION STARTS LATER THAN THE FIRST OBS AVAILABLE
    ###########################################################################

    idx_site_no_obs_at_first_obs_date = np.where( np.isnan(model.t_obs[0,:,0])  )[0]

    if idx_site_no_obs_at_first_obs_date.shape[0] != 0:
        print("[WARNING] some sites do not start at the first observation date.")
        print("[WARNING] either use --build 3 option or remove those sites using --exclude_gps option")
        print("%s" % ( textwrap.fill(" ".join( model.name_obs[ idx_site_no_obs_at_first_obs_date  ] ) , width=terminal_width ) ))
        import sys
        sys.exit()


    ###########################################################################
    # COMPONENTS & RAKE
    ###########################################################################
    
    # components
    if model.UP:
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

    nrake = np_rake.shape[0] 
    
    # debug
    if model.debug:
        print("-- np_rake:" , np_rake )
        
    ###########################################################################
    # OBSERVATION TENSOR
    ###########################################################################

    # let's try to make a data tensor
    # D(i,j,k) would be the displacement divided by the uncertainty
    # observation time index i
    # for site j
    # for component k
    # no data would be Nan
    
    
    ###########################################################################
    # AUXILIARY INDEX
    ###########################################################################
    
    print('-- list of observation dates')
    for i in np.arange( model.t_obs.shape[0] ):
        print("%s" % ( at.seconds2datetime( model.np_obs_date_s[i] ).isoformat(' ')   ) )
    
    print("-- Building H_i_k  : dictionary providing the list of indexes k for observation time for each model time step i")
    # H_i_k  : dictionary providing the list of indexes k for observation time for each model time step i
    H_i_k = {}

    
    for model_time_step in np.arange( model.np_model_date_s.shape[0]-2 , -1, step = -1 ):
        sdate_s = model.np_model_date_s[ model_time_step ]
        edate_s = model.np_model_date_s[ model_time_step + 1 ]
        print("--- model_time_step: %04d %s -> %s " % ( model_time_step, at.seconds2datetime(sdate_s).isoformat(' '), at.seconds2datetime(edate_s).isoformat(' ') ) )
        H_i_k[ model_time_step ] = []
        
        for k in np.arange( model.np_obs_date_s.shape[0]):
            if ( ( model.np_obs_date_s[k] > sdate_s ) and ( model.np_obs_date_s[k] <= edate_s ) ):
                H_i_k[ model_time_step ].append( k )
                 
        H_i_k[ model_time_step ] = sorted( H_i_k[ model_time_step ] )
        
    if model.verbose:
        print("--- list of observation index per model time step")
        for key, value in H_i_k.items():
            
            print("--- model time step %d: %s " % (key,value) )
            

    # H_k_obs: dictionary providing the list of GNSS site indexes available at observation time t_k

    print("-- Building H_k_obs: dictionary providing the list of GNSS site indexes available at observation time t_k")
    
    H_k_obs = {}
    for k in np.arange( model.np_obs_date_s.shape[0]):
        H_k_obs[ k ] = np.where( np.isfinite( model.t_obs[ k,:,0 ]))

    if model.debug:
        print('H_k_obs')
        for key, value in H_k_obs.items():
            print(key,value)

        for klkl in np.arange( model.name_obs.shape[0] ):
            print("%04d %s" % ( klkl, name_obs[ klkl ] ))

    ## stupid?! H_k_tk : dictionary providing the observation date tk at the observation index k


    ###########################################################################
    # ATPA ATBP and various variables initialization
    ###########################################################################

    n_model_time_step = model.np_model_date_s.size - 1
    nfaults = model.green.shape[1]
    model.n_model_time_step = n_model_time_step
    
    print("-- Initializing ATPA & ATPB")
    
    # ATPA
    
    memory_usage = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024 / 1024
    print("-- memory usage before ATPA initialization in Gb Linux/Mac OS X : %.1lf / %.1lf" %  (memory_usage,memory_usage/1024) )

    n_atpa = nfaults * nrake * n_model_time_step + n_constant
    ATPA = np.zeros( (n_atpa,n_atpa) )

    memory_usage = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024 / 1024
    print("-- memory usage after ATPA initialization in Gb Linux/Mac OS X : %.1lf / %.1lf" %  (memory_usage,memory_usage/1024) )
    
    # ATPB
    
    ATPB = np.zeros( (n_atpa,1) )

    memory_usage = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024 / 1024
    print("-- memory usage after ATPB initialization in Gb Linux/Mac OS X : %.1lf / %.1lf" %  (memory_usage,memory_usage/1024) )

    # COLUMN OF ATPA CORRESPONDING TO A TIME STEP

    # number of columns for C
    ncc = nfaults * nrake
    print('-- nrake' , nrake)
    print('-- nfaults' , nfaults)
    print("-- Initializing the individual model time step matrices")
    C = np.zeros( ( ncc * n_model_time_step , ncc ) )

    # GAMMA MODEL STEP DURATION VECTOR
    
    print("-- Initializing the vector of time duration")
    GAMMA = np.diff( model.np_model_date_s , 1 )
    
    ###########################################################################
    # loop on model time steps
    ###########################################################################

    # RHS sub-vector
    Bi = np.zeros( ( nfaults*nrake , 1 ) )
    
    # loop is made reverse in time
    for model_time_step in np.arange( model.n_model_time_step-1, -1, -1 ):
#    for model_time_step in np.arange( np_model_date_s.shape[0]-2 , -1, step = -1 ):
  
        if model.verbose:
            print("-- Building equations for model time step: %04d" % model_time_step )
        
        # Get the list of the k-index of observation time within the current model time step
        
        lk = H_i_k[ model_time_step ]
        
        
        #######################################################################
        # loop on k observation time within the current model time step
        #######################################################################
        for k in lk:
            
            # get tk time for observation k
            
            tk = model.np_obs_date_s[ k ]
            print("-- model time step: %04d observation time step: %04d %s" % (model_time_step , k, at.seconds2datetime(tk).isoformat(' ') ) )
            
            # get Gkt and Gk2
            if model.debug:
                print("-- model time step: %04d observation time step: %04d Gkt from pyeq.lib.elastic_tensor.shrink"  % (model_time_step , k ) )
            
            # extract Gkt as a 2D matrix for available sites with the chosen components at date tk
            Gkt = pyeq.lib.elastic_tensor.shrink( model.green, lfaults=None, lobs=H_k_obs[k], lcomponent=np_component, lrake=np_rake).T

            # Get dk, observation vector at date tk for available sites with the chosen components
            if model.debug:
                print("-- model time step: %04d observation time step: %04d get normalized observation using pyeq.lib.obs_tensor.dk_from_obs_tensor"  % (model_time_step , k ) )
            
            Dk = model.t_obs[k][ H_k_obs[k] ][:,np_component].T.flatten()

            # normalize the observation vector
            Sdk = model.t_obs[k][ H_k_obs[k] ][:,np_component +3 ].T.flatten()
            Dk  = Dk / Sdk

            # normalize the Green matrix
            Gkt = (Gkt / Sdk ) 
            Gk2 = np.dot( Gkt , Gkt.T )
            
            # handle the Gamma functions
            # get GAMMA, change for the current time step, remove time steps after the current one to ensure size is compatible with C
            # we choose as a unit for parameters delta_m, slip rate expressed in mm/day
            if model.debug:
                print("-- model time step: %04d observation time step: %04d gamma_k: vector of duration"  % (model_time_step , k ) )

            gamma_k = np.copy( GAMMA )
            gamma_k[ model_time_step ] = ( tk - model.np_model_date_s[ model_time_step ] )
            gamma_k = gamma_k / 86400.

            # Update C (See "Algorithm to build the observation normal system in paper")
            if model.debug:
                print("-- model time step: %04d observation time step: %04d Individual column matrix for ATPA"  % (model_time_step , k ) )
            C = C + gamma_k[ model_time_step ] * pyacs.lib.glinalg.odot( gamma_k[:model_time_step+1] , pyacs.lib.glinalg.repeat_matrix_in_col(Gk2 , model_time_step+1 ) )
            
            # Update Bk
            if model.debug:
                print("-- model time step: %04d observation time step: %04d Individual vector for ATPB"  % (model_time_step , k ) )
            Bi = Bi + gamma_k[ model_time_step ] * np.dot( Gkt , Dk.reshape(-1,1) )
            
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

        # Updating ATPB
        if model.debug:
            print("-- model time step: %04d Updating ATPB"  % (model_time_step ) )
        
        ATPB[ model_time_step*ncc:(model_time_step+1)*ncc ] = Bi 

        #######################################################################
        # end of the current model time step
        #######################################################################
        
        # Remove the last rows of C corresponding to the current time step
        if model.debug:
            print("-- model time step: %04d Resizing individual column matrix for next model time step"  % (model_time_step ) )
        
        C = C[:-ncc,:]
        
    return ATPA , ATPB



