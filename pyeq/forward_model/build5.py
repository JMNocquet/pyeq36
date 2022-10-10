def build5( model ):
    """
    New scheme to build the linear system.
    Corrects a bug from build4 in the case of variable model time step
    """
    
    ###########################################################################
    # IMPORT
    ###########################################################################
    
    import numpy as np

    import pyeq.forward_model.C_B_at_k
    import pyacs.lib.astrotime as at
    from progress.bar import Bar

    import pyeq.message.message as MESSAGE
    import pyeq.message.verbose_message as VERBOSE
    import pyeq.message.error as ERROR
    import pyeq.message.warning as WARNING
    import pyeq.message.debug_message as DEBUG

    VERBOSE("Using build5")

    np.set_printoptions(precision=2, suppress=False, linewidth=150 )

    DEBUG("model.offset: %.2lf " % model.offset )

    ### TO BE DONE
    ### RAKE handling - check whether GREEN has been modified before build2 call
    ### ON GOING: PARAMETER TRANSLATION FOR CONJUGATE RAKE AND/OR ORIGIN TIME CONSTANT
    
    if model.offset == 0:
        n_constant = 0
        model.nconstant = n_constant
    else:
        n_constant = model.ncomponent * model.green.shape[0]
        model.nconstant = n_constant
    
    VERBOSE("number of offsets at origin time: %d" % n_constant)

    ###########################################################################
    # TEMPORARY - TO BE CHANGED WHEN VARIABLE RAKE IMPLEMENTED
    ###########################################################################
    RAKE = False

    ###########################################################################
    # IMPORT
    ###########################################################################

    import numpy as np
    import pyeq.forward_model
    #np.set_printoptions(precision=2, suppress=True)

    ###########################################################################
    # IMITIALIZE COMPONENTS, RAKE AND GREEN
    ###########################################################################
    
    # components
    if model.up:
        np_component = np.array( [0,1,2] )
    else:
        np_component = np.array( [0,1] )
    
    DEBUG("np_component: %s" % np_component )
    
    # rake
    if RAKE:
        np_rake = np.array( [0,1] )
    else:
        np_rake = np.array( [0] )

    model.np_rake = np_rake

    model.nrake = np_rake.shape[0] 
    
    # debug
    DEBUG("np_rake:%s" % np_rake )

    # Get the Green tensor and reshape it as a 2D array
    # Cd-1/2 still missing here, handled later
    G = model.green[ : ][:, :][:,:,np_component][:,:,:,np_rake]
    G = G.reshape( G.shape[0] , G.shape[1]* G.shape[3] , G.shape[2] )
    G = np.moveaxis(G,[0,1,2],[0,2,1]) 
    model.G = G.reshape( G.shape[0]*G.shape[1] , G.shape[2] ) 
    
    # GAMMA_MATRIX
    GAMMA_MATRIX = pyeq.forward_model.make_GAMMA_MATRIX(model.t_obs, model.np_obs_date_s, model.np_model_date_s)

    # DICTIONNARY OF OBS INDEX k AT EACH MODEL STEP i
    H_i_k = pyeq.forward_model.make_auxiliary_index(model)
    

    ###########################################################################
    # INITIALIZE LHS & RHS NORMAL SYSTEM
    ###########################################################################
    n_N = model.nfaults * model.nrake * model.nstep + model.nconstant
    model.N = np.zeros( ( n_N , n_N ) )
    model.Nd = np.zeros( n_N )
    
    ###########################################################################
    # BUILDS THE LINEAR SYSTEM GOING BACKWARDS IN MODEL TIME STEP
    ###########################################################################
                                                               
    for model_time_step in np.flip( np.arange( model.nstep ) ):
        
        ssdate_s = model.np_model_date_s[ model_time_step  ]
        eedate_s = model.np_model_date_s[ model_time_step + 1 ]
        
        VERBOSE("model #%04d %s to %s" % \
              ( model_time_step , at.seconds2datetime(ssdate_s).isoformat(' '), at.seconds2datetime(eedate_s).isoformat(' ') ) )
        
        # Get the list of the k-index of observation time within the current model time step
        lk = H_i_k[ model_time_step ] 
        
        for k in sorted( lk ):

            tk = model.np_obs_date_s[ k ]
            DEBUG("model #%04d observation #%04d %s" % (model_time_step , k, at.seconds2datetime(tk).isoformat(' ') ) )

            # inverse square root data covariance W_Sdk        
            dk  = model.t_obs[k ,:,np_component].T
            sdk = model.t_obs[k ,:,np_component+3].T
            
            W_Sdk = 1./ np.nan_to_num(sdk , nan=1.)
            DK    = np.nan_to_num( dk, nan=0. ) * W_Sdk
            CDK   = np.nan_to_num( dk*0.+1., nan=0. ) * W_Sdk 

            # Update the GGAMMA_MATRIX for the current time step model_time_step
            #frac_time_step = ( model.np_obs_date_s[k] - model.np_model_date_s[ model_time_step ]  ) / 86400.

            #frac_time_step = ( model.np_obs_date_s[k] - model.np_model_date_s[ model_time_step ]  ) \
            #                / ( model.np_model_date_s[ model_time_step+1 ] - model.np_model_date_s[ model_time_step ] )
        
            #GAMMA_MATRIX[ : , model_time_step ] =  np.heaviside( GAMMA_MATRIX[ : , model_time_step ] - (1. - frac_time_step) , 0  )  \
            #                                       * ( GAMMA_MATRIX[ : , model_time_step ] - (1. - frac_time_step) )  \

            for idx_site in np.arange( model.t_obs.shape[1] ):
                # get the first not nan observation for the current site
                try:
                    idx_first_obs = np.where( np.isfinite( model.t_obs[:,idx_site,0]) )[0][0]
                except:
                    continue
                # get the first available observation for the current site
                date_first_obs = model.np_obs_date_s[ idx_first_obs ]
                # update with the duration of site in the current time step
                if date_first_obs > model.np_obs_date_s[k]:
                    GAMMA_MATRIX[ idx_site , model_time_step ] = 0.
                else:
                    GAMMA_MATRIX[ idx_site , model_time_step ] = \
                    ( model.np_obs_date_s[k] - np.max( [ date_first_obs , model.np_model_date_s[ model_time_step ] ] ) ) / 86400.
            

            G = pyeq.forward_model.make_G_at_k(model, model_time_step, CDK, GAMMA_MATRIX)

            ###################################################################
            # loop on the time step for the current date
            ###################################################################
            
            # constant term
            if model.nconstant > 0:
                model.N[ -model.nconstant: , -model.nconstant: ] \
                    += np.diag( CDK.flatten()**2 )

                model.Nd[ -model.nconstant: ] += DK.flatten()
            
            bar = Bar('', max= (model_time_step + 2)* (model_time_step + 1) /2 , suffix = '%(percent).1f%% - %(elapsed)ds')

            for i in np.flip( np.arange( model_time_step + 1 ) ):
                
                # Nd
                model.Nd[ i*model.nfaults : (i+1)*model.nfaults ] += np.dot( G[i].T , DK.flatten().reshape(-1,1) ).flatten()

                # constant term
                if model.nconstant > 0:
                    model.N[ -model.nconstant: , i*model.nfaults : (i+1)*model.nfaults ] \
                        += np.dot( np.diag( CDK.flatten() ) , G[i] )     


                for j in np.flip( np.arange( i+1 ) ):
                    bar.next()
                    #print("-- Filling contribution to model.N at (%d,%d)" % (i,j))
                    
                    model.N[ i*model.nfaults : (i+1)*model.nfaults  , j*model.nfaults : (j+1)*model.nfaults ] \
                        += np.dot( G[i].T , G[j] )


                    #print(model.N)

            bar.finish()

    # Obs. eq. for t0 for constant
    if model.nconstant > 0:
        MESSAGE("Coordintes at t=0 will be adjusted.")
        dk  = model.t_obs[0 ,:,np_component].T
        sdk = model.t_obs[0 ,:,np_component+3].T
        
        # Normalize observation
    
        W_Sdk = 1./ np.nan_to_num(sdk , nan=1.)
        DK    = np.nan_to_num( dk, nan=0. ) * W_Sdk
        CDK   = ( np.nan_to_num( dk*0.+1., nan=0. ) * W_Sdk ).flatten()
    
        VERBOSE('adding constant column to the normal matrix')
        model.N[ -model.nconstant: , -model.nconstant: ] += np.diag( CDK.flatten()**2 )
    

    MESSAGE('Making the matrix symmetric. Could take some time...')

    bar = Bar('', max= ( model.nstep )* ( model.nstep - 1) /2 , suffix = '%(percent).1f%% - %(elapsed)ds')

    for i in np.flip( np.arange( model.nstep ) ):
        for j in np.flip( np.arange( i ) ):
            
            bar.next()

            model.N[ j*model.nfaults : (j+1)*model.nfaults  , i*model.nfaults : (i+1)*model.nfaults ] \
                =  model.N[ i*model.nfaults : (i+1)*model.nfaults  , j*model.nfaults : (j+1)*model.nfaults ]
    
    bar.finish()
    
    if model.nconstant > 0:
        for i in np.flip( np.arange( model.nstep + 1 ) ):
            model.N[ i*model.nfaults : (i+1)*model.nfaults  , -model.nconstant: ] = model.N[ -model.nconstant: , i*model.nfaults : (i+1)*model.nfaults ].T
    
    
    #model.N = pyacs.lib.glinalg.symmetrize(model.N,'tril')
    
    #model.nconstant = 0


    return model.N , model.Nd 
