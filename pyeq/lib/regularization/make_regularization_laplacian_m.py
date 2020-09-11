def make_regularization_laplacian( model ):
    """
    Build a Laplacian regularization matrix
    """
    
    # IMPORT
    
    import numpy as np
    import scipy.linalg    
    import pyacs.lib.glinalg
    import pyeq.lib.regularization
    import pyeq.lib.log
    from time import time
    from progress.bar import Bar
   
    # computes model step duration vector in days delta_time_in_days
    print("-- Computing model step duration vector" )
    delta_time_in_days = np.diff( model.np_model_date_s / 86400. )

###############################################################################
# COMPUTE THE SMOOTHING OPERATOR
###############################################################################

    # compute the individual matrix_d_over_dc
    matrix_d_over_dc = model.dm / model.dc

    # building the step regularization matrix
    #L = 1./ ( 1 + matrix_d_over_dc )
    
    # normalize by the sum of values per row
    # remove diagonal elements
    #np.fill_diagonal( L, 0. )
    # compute sum per row
    #S = np.sum( L , axis=1 )
    # normalize
    #L = ( L.T / S ).T
    # re-inject diagonal elements
    #np.fill_diagonal( L , -1./S )


###############################################################################
# SPATIAL SMOOTHING
###############################################################################
    print("-- Laplacian regularization - no temporal smoothing")
    
    # building the smoothing matrix including the individual weights
    #L = 1./ ( 1 + matrix_d_over_dc )**2
    L = np.exp( -matrix_d_over_dc )
    # building the Normal matrix for smoothing
    N_S = np.zeros( ( model.nfaults, model.nfaults ) )
    # non-diagonal elements
    for i in np.arange( model.nfaults ):
        for j in np.arange( i+1 , model.nfaults ):
        # non diagonal are directly the weights with a minus sign
            N_S[ i , j ] = N_S[ i , j ] - L[ i , j ]
            N_S[ j , i ] = N_S[ j , i ] - L[ j , i ]
            # diagonal
            N_S[ i,i ] = N_S[ i,i ] + L[ i , j ]
            N_S[ j,j ] = N_S[ j,j ] + L[ i , j ]
        
        N_S = N_S / np.mean( np.diag(N_S) ) * 1./model.sigma_time**2


        # damping
        np.fill_diagonal(N_S,np.diag(N_S)+1./model.sigma_time**2)
         
        # put in the model Normal matrix
        
        model.N =  model.N + scipy.linalg.block_diag(*([N_S] * model.nstep))

###############################################################################
# CASE TEMPORAL SMOOTHING - FOR NOW SMOOTHING IS INDEPENDANT FOR DAMPING
###############################################################################

    if model.tau != 0:
        print("-- Adding temporal smoothing. Correlation time is tau = %.2lf days" % model.tau )

        print("-- Computing the model step time distance matrix" )
        delta_day_2D = pyeq.lib.regularization.time_to_2D_diff_time_in_days(  (model.np_model_date_s / 86400.)[:-1] + delta_time_in_days / 2 )
        delta_day_2D = delta_day_2D / model.tau
        #Dm = model.dm / model.dc

        print("-- Computing the smoothing normal matrix" )
        bar = Bar('', max= model.nstep  , suffix = '%(percent).1f%% - %(eta)ds')
        # building the Normal matrix for smoothing
        N_S = np.zeros( ( model.nparameters, model.nparameters ) )
        # non-diagonal elements
        for i in np.arange( model.nstep):
            bar.next()
            LOCAL = pyacs.lib.glinalg.matrix_from_pattern( Dm , np.zeros(( 1, model.nstep-i ))+1. )
            LOCAL += pyacs.lib.glinalg.matrix_from_pattern( np.zeros(( model.nfaults , model.nfaults )) +1. , np.atleast_2d( delta_day_2D[i,i:] )  )
            
            LOCAL = np.exp( - LOCAL )

            N_S[ i*model.nfaults:(i+1)*model.nfaults , i *model.nfaults: ] = -LOCAL
            N_S.T[ i*model.nfaults:(i+1)*model.nfaults , i *model.nfaults: ] = -LOCAL

        bar.finish()
                
        # diagonal elements
        print("-- Computing diagonal terms" )
        np.fill_diagonal( N_S , 0 )
        norm = -np.sum( N_S , axis=0 )
        np.fill_diagonal( N_S , norm )
        
        # normalize                
        print("-- Normalizing diagonal terms by the mean of diagonal elements" )
        #trace_N_S = np.sum( np.diag(N_S) )

        # we want to equate damping and smoothing eight
        #print( trace_N_S )
        #trace_damping = 1./model.sigma_time**2 * model.nfaults
        
        #N_S = N_S / np.mean( np.diag(N_S) )
        N_S = N_S / np.mean( np.diag(N_S) ) * 1./model.sigma_time**2
        
        # damping
        print("-- Adding damping to smoothing" )
        np.fill_diagonal(N_S,np.diag(N_S)+1./model.sigma_time**2)
         
        # put in the model Normal matrix
        
        print("-- Adding regularization to the normal system" )
        model.N =  model.N + N_S
         
        return model
        
        

        for i in np.arange( model.nstep ):
            for j in np.arange( i, model.nstep ):
                
                bar.next()
                
                # block matrix 1./( 1+ d/dc + t/tau ) 
                R_step = 1./ ( matrix_d_over_dc + delta_day_2D[i,j]/model.tau + 1. )
                # normalize with sigma_time
                if 'time_constant' in model.sigma_type:
                    #print("-- sigma_time constant case")
                    R_step = 1./ ( model.sigma_time * delta_time_in_days[i] * model.sigma_time * delta_time_in_days[j] ) * R_step
                else:
                    #print("-- sigma_time variable case")
                    R_step = 1./ ( model.sigma_time[i] * delta_time_in_days[i] * model.sigma_time[j] * delta_time_in_days[j] ) * R_step
                # normalize with sigma_fault
                if 'fault_constant' in model.sigma_type:
                    pass
                else:
                    R_step = R_step * np.outer( model_sigma_fault , model_sigma_fault )
                
                model.N[i*model.nfaults:(i+1)*model.nfaults , j*model.nfaults:(j+1)*model.nfaults ] += R_step

                model.N[j*model.nfaults:(j+1)*model.nfaults , i*model.nfaults:(i+1)*model.nfaults ] += R_step
       
        bar.finish()

###############################################################################
# CASE NO TEMPORAL SMOOTHING
###############################################################################
    else:
        print("-- Computing the individual model_step inverse spatial correlation matrix (inv_corr_step)")
        

    
    return model
