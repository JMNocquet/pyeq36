def make_regularization_laplacian( model ):
    """
    Build a Laplacian regularization matrix
    """

    # NOTES
    # BETTER IMPLEMENT SMOOTHING KERNEL WITH CONF
    
    
###############################################################################
    # IMPORT
###############################################################################
    
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
# CASE NO TEMPORAL SMOOTHING
###############################################################################
#    if model.tau == 0:
        
    print("-- Laplacian regularization - no temporal smoothing")
    
    # make a selection of adjacent triangles
    # compute a critical length to tell that triangles are adjacent
    dsel = np.min( np.sort( model.dm[model.dm> 1.2*np.median( np.sort( model.dm ) , axis=0)[1]] ))
    
    # set all non-adjacent triangles to 0
    L = np.where(model.dm > desl , 0 , model.dm ) 
    
    # renormalize to the unit space
    L = L / model.dc
    
    # building the smoothing matrix including the individual weights
    # building the Normal matrix for smoothing
    N_S = np.zeros( ( model.nfaults, model.nfaults ) )
    # non-diagonal elements
    for i in np.arange( model.nfaults ):
        # diagonal
        N_S[ i,i ] = 1.
        for j in np.arange( i+1 , model.nfaults ):
        # non diagonal are directly the weights with a minus sign
            N_S[ i , j ] = N_S[ i , j ] - L[ i , j ]
            N_S[ j , i ] = N_S[ j , i ] - L[ j , i ]
    

    # renormalization
    
    STD = np.sqrt( np.diag(N_S) )
    N_S = N_S / np.outer(STD , STD )

    if 'fault_variable' in model.sigma_type:
        N_S = N_S / model.sigma_fault**2

    # 
    #N_S = N_S  * pyacs.lib.glinalg.matrix_from_pattern( np.ones( (1,model.nfaults )) * 1./model.sigma_time**2 , delta_time_in_days.reshape(1,-1) )
    # temporary solution - constant time step - constant everything

    if 'time_variable' in model.sigma_type:
        N_S = N_S / model.sigma_fault**2

        
        for i in np.arange( model.nstep ):
            model.N[ i*model.nfaults:(i+1)*model.nfaults , i*model.nfaults:(i+1)*model.nfaults ] += \
            N_S * 1./model.sigma_time[i]**2 * 1./delta_time_in_days[i] + \
            1./model.sigma_time[i]**2 * 1./delta_time_in_days[i]           

    else:
        for i in np.arange( model.nstep ):
            model.N[ i*model.nfaults:(i+1)*model.nfaults , i*model.nfaults:(i+1)*model.nfaults ] += \
            N_S * 1./model.sigma_time**2 * 1./delta_time_in_days[i] + \
            1./model.sigma_time**2 * 1./delta_time_in_days[i]
            
            
        #
        #N_S = N_S * 1./model.sigma_time**2

        

        # damping
        #np.fill_diagonal(N_S,np.diag(N_S)+1./model.sigma_time**2)
         
        # put in the model Normal matrix
        #model.N =  model.N + scipy.linalg.block_diag(*([N_S] * model.nstep))
         
        #return model
        

###############################################################################
# CASE TEMPORAL SMOOTHING - FOR NOW SMOOTHING IS INDEPENDANT FOR DAMPING
###############################################################################

    if model.tau != 0:
        print("-- Building the model regularization matrix with temporal smoothing. Correlation time is tau = %.2lf days" % model.tau )

        print("-- Computing the model step time distance matrix" )
        delta_day_2D = pyeq.lib.regularization.time_to_2D_diff_time_in_days(  (model.np_model_date_s / 86400.)[:-1] + delta_time_in_days / 2 )
        
        
        # set the elements at twice tau to 0
        lindex =np.where( delta_day_2D > 2 * model.tau )
        delta_day_2D[ lindex ] = 0.
        
        # apply kernel
        delta_day_2D = np.exp( -( delta_day_2D / model.tau )**2 / 2 )**2

        # normalize for the normal system
        np.fill_diagonal( delta_day_2D , 0 )
        np.fill_diagonal( delta_day_2D , np.sum( delta_day_2D , axis = 0 ) )
        STD = np.sqrt( np.diag( delta_day_2D ) )
        delta_day_2D = delta_day_2D / np.outer(STD , STD )

        # index for loop
        lindex = np.where( np.triu(delta_day_2D , 0 ) > 0 )
        
        # loop on these elements

        
        for model_step in np.arange( len(lindex[0]) ):
            # index of the block
            i = lindex[0][ model_step ]
            j = lindex[1][ model_step ]

            
            STDi = np.sqrt( np.diagonal( model.N[ i*model.nfaults:(i+1)*model.nfaults , i*model.nfaults:(i+1)*model.nfaults ]  ) )
            STDj = np.sqrt( np.diagonal( model.N[ j*model.nfaults:(j+1)*model.nfaults , j*model.nfaults:(j+1)*model.nfaults ]  ) )
            
            
            
            np.fill_diagonal( model.N[ i*model.nfaults:(i+1)*model.nfaults , j*model.nfaults:(j+1)*model.nfaults ] , \
                              np.diagonal( model.N[ i*model.nfaults:(i+1)*model.nfaults , j*model.nfaults:(j+1)*model.nfaults ] ) + delta_day_2D[i,j] * (STDi * STDj) )
        
            np.fill_diagonal( model.N[ j*model.nfaults:(j+1)*model.nfaults , i*model.nfaults:(i+1)*model.nfaults ] , \
                              np.diagonal( model.N[ j*model.nfaults:(j+1)*model.nfaults , i*model.nfaults:(i+1)*model.nfaults ] ) + delta_day_2D[i,j] * (STDi * STDj) )

    
###############################################################################
# HANDLE ORIGIN TIME CONSTANTS
###############################################################################
    if model.offset > 0:
        np.fill_diagonal( model.N[ -model.nconstant: , -model.nconstant: ] , np.diagonal( model.N[ -model.nconstant: , -model.nconstant: ] ) + 1./model.offset**2 ) 
    
    return model
        
