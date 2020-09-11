###############################################################################
def make_model_covariance_02( model ):
###############################################################################
    """
    Create a model covariance matrix Cm
    
    cov_02 assumes:
    - temporal smoothing: specified through critical duration model.tau. Could be 0 for no temporal smoothing.
    - model.np_sigma is a scalar controlling the a priori average slip rate per day (sigma**2 in mm**2 /day )
    - model.np_sigma is rescaled according to the model step duration
    
    - with this formulation, np_sigma can also be a vector of length model.nfaults
    - finally model.np_sigma can be a matrix whose columns are sigma vectors for each model time step
    """
    
    ###########################################################################
    # IMPORT
    ###########################################################################
    
    import numpy as np
    import scipy.linalg    
    import pyacs.lib.glinalg
    import pyeq.lib.regularization
    import pyeq.lib.log
    from time import time
    from progress.bar import Bar

    ###########################################################################
    # INITIALIZE SOME PARAMETERS
    ###########################################################################

    # first get the correlation matrix
    corr_step = pyeq.lib.regularization.make_model_spatial_correlation_matrix( model.dm , model.dc, cm_type =  model.cm_type )

    # computes normalization factor (not squared, d0 (min characteristic distance in model.sgeometry ) / dc )
    if model.cm_norm == 'd0':
        nf = pyeq.lib.regularization.normalization_factor( model.sgeometry , model.dc )
    else:
        nf = 1.
    
    print("-- Normalization factor accounting for fault discretization: %.3lf " % nf )

    # computes model step duration vector in days delta_time_in_days
    print("-- Computing model step duration vector" )
    delta_time_in_days = np.diff( model.np_model_date_s / 86400. )

    # handling the np_sigma
    
    print('-- Handling model covariance diagonal (sigma)')

    model.np_sigma = None
    if model.sigma_type == 'time_constant/fault_constant':
        model.np_sigma = np.zeros( (model.nparameters) ) + model.sigma_time

    if model.sigma_type == 'time_variable/fault_constant':
        model.np_sigma = pyacs.lib.glinalg.matrix_from_pattern(np.array([ model.sigma_time]), np.zeros((1, model.nstep ) ) + model.sigma_fault ).flatten() 
    
    if model.sigma_type == 'time_constant/fault_variable':
        model.np_sigma = pyacs.lib.glinalg.matrix_from_pattern(np.array([ model.sigma_fault]), np.zeros((1, model.nstep ) ) + model.sigma_time ).flatten() 

    if model.sigma_type == 'time_variable/fault_variable':
        model.np_sigma = pyacs.lib.glinalg.matrix_from_pattern(np.array([ model.sigma_fault]), np.array( [model.sigma_time] ) ).flatten()

    # fault discretization lengthscale normalization
    np_sigma = model.np_sigma * nf 

###############################################################################
# CASE TEMPORAL SMOOTHING
###############################################################################

    if model.tau != 0:
        print("-- Building the model covariance matrix CM with temporal smoothing. Correlation time is tau = %.2lf days" % model.tau )
        
    
        # spatial correlation at all time step
        ones = np.zeros( ( model.nstep , model.nstep ) ) + 1.
        print("-- Filling the model correlation matrix with spatial correlation blocks" )
        Corr_spatial =  pyacs.lib.glinalg.matrix_from_pattern( corr_step , ones )
        
        # compute the time distance among model time step matrix
        print("-- Computing the model step time distance matrix" )
        delta_day_2D = pyeq.lib.regularization.time_to_2D_diff_time_in_days(  (model.np_model_date_s / 86400.)[:-1] + delta_time_in_days / 2 )
    
        # temporal correlation at all time step
        print("-- Computing the model step time correlation matrix" )
        corr_time_2D = np.exp( -delta_day_2D / model.tau )

        # temporal correlation matrix            
        print("-- Filling the model correlation matrix with time correlation blocks" )
        ones = np.zeros( ( model.nfaults , model.nfaults ) ) + 1.
        Corr_time =  pyacs.lib.glinalg.matrix_from_pattern( ones , corr_time_2D )
        
        # spatial and time correlation matrix
        print("-- Merging spatial and time correlation matrices" )
        Corr = Corr_spatial * Corr_time
        
        # get the model covariance matrix
        print('-- Converting correlation to covariance with shape:' , Corr.shape )
        
        CM = pyacs.lib.glinalg.corr_to_cov(Corr, np_sigma )

        # get the inverse covariance matrix  
        print('-- Computing the inverse model covariance matrix of shape: ' , CM.shape )
        print('-- And build the normal matrix' )
        model.N =  model.N + pyacs.lib.glinalg.syminv( CM )

###############################################################################
# CASE NO TEMPORAL SMOOTHING
###############################################################################
    else:
        # compute the inversion correlation matrix
        # !!! THERE IS A SURPRISE HERE
        # THE COVARIANCE MATRIX APPEARS TO BE BETTER CONDITIONED
        # IF THE COVARIANCE IS INVERSED COMPARED TO THE CORRLATION
        # MATRIX INVERSED AND THEN MULTIPLIED BY 1/sigma**2
        # THIS LEADS TO SUB-OPTIMAL CODING
        print("-- Computing the individual model_step inverse spatial correlation matrix (inv_corr_step)")

        #######################################################################        
        if model.sigma_type == 'time_constant/fault_constant':
        #######################################################################        
            print("-- Computing the time step covariance matrix - sigma_type is %s " % ( model.sigma_type ) )

            if np.max( np.fabs( np.diff( delta_time_in_days ) ) ) == 0:
                print("-- Constant time step case ")
                normalization_factor = delta_time_in_days[0] * model.sigma_time**2 * nf
                cm_step = corr_step * normalization_factor
                inv_cm_step = pyacs.lib.glinalg.syminv( cm_step )
                model.N = model.N + scipy.linalg.block_diag(*([inv_cm_step] * model.nstep))
            else:
                print("-- Variable time step case ")
                # THIS COULD BE IMPROVED: HERE THIS REQUIRES NSTEP INVERSION
                # ALTHOUGH MANY OF THEM ARE PROBABLY THE SAME

                bar = Bar('', max= model.nstep , suffix = '%(percent).1f%% - %(eta)ds')

                for i in np.arange( model.nstep ):

                    bar.next()
                    
                    normalization_factor = delta_time_in_days[i] * model.sigma_time**2 * nf
                    cm_step = corr_step * normalization_factor
                    inv_cm_step = pyacs.lib.glinalg.syminv( cm_step )
                    model.N[ i * model.nfaults:(i+1) * model.nfaults , i * model.nfaults:(i+1) * model.nfaults ] = \
                        model.N[ i * model.nfaults:(i+1) * model.nfaults , i * model.nfaults:(i+1) * model.nfaults ] + inv_cm_step
                
                bar.finish()

        #######################################################################        
        if model.sigma_type == 'time_variable/fault_constant':
        #######################################################################        
            print("-- Computing the time step covariance matrix - sigma_type is %s " % ( model.sigma_type ) )

            bar = Bar('', max= model.nstep , suffix = '%(percent).1f%% - %(eta)ds')

            for i in np.arange( model.nstep ):

                bar.next()
                normalization_factor = delta_time_in_days[i] * model.sigma_time[i]**2 * nf
                cm_step = corr_step * normalization_factor
                inv_cm_step = pyacs.lib.glinalg.syminv( cm_step )
                model.N[ i * model.nfaults:(i+1) * model.nfaults , i * model.nfaults:(i+1) * model.nfaults ] = \
                    model.N[ i * model.nfaults:(i+1) * model.nfaults , i * model.nfaults:(i+1) * model.nfaults ] + inv_cm_step
            
            bar.finish()

        #######################################################################        
        if model.sigma_type == 'time_constant/fault_variable':
        #######################################################################        
            print("-- Computing the time step covariance matrix - sigma_type is %s " % ( model.sigma_type ) )

            if np.max( np.fabs( np.diff( delta_time_in_days ) ) ) == 0:
                print("-- Constant time step case ")
                std_vector = np.sqrt( delta_time_in_days[0]) * model.time_sigma * model.sigma_fault * nf
                inv_cm_step = pyacs.lib.glinalg.corr_to_cov( corr_step,  std_vector )
                model.N = model.N + scipy.linalg.block_diag(*([inv_cm_step] * model.nstep))
            else:
                print("-- Variable time step case ")
                # THIS COULD BE IMPROVED: HERE THIS REQUIRES NSTEP INVERSION
                # ALTHOUGH MANY OF THEM ARE PROBABLY THE SAME

                bar = Bar('', max= model.nstep , suffix = '%(percent).1f%% - %(eta)ds')

                for i in np.arange( model.nstep ):

                    bar.next()
                    std_vector = np.sqrt( delta_time_in_days[i] ) * model.fault_sigma * model.sigma_time * nf
                    cm_step = pyacs.lib.glinalg.corr_to_cov( corr_step,  std_vector )
                    inv_cm_step = pyacs.lib.glinalg.syminv( cm_step )
                    model.N[ i * model.nfaults:(i+1) * model.nfaults , i * model.nfaults:(i+1) * model.nfaults ] = \
                        model.N[ i * model.nfaults:(i+1) * model.nfaults , i * model.nfaults:(i+1) * model.nfaults ] + inv_cm_step
                
                bar.finish()

        #######################################################################        
        if model.sigma_type == 'time_variable/fault_variable':
        #######################################################################        
            print("-- Computing the time step covariance matrix - sigma_type is %s " % ( model.sigma_type ) )

            bar = Bar('', max= model.nstep , suffix = '%(percent).1f%% - %(eta)ds')

            for i in np.arange( model.nstep ):

                bar.next()
                std_vector = np.sqrt( delta_time_in_days[i] ) * model.fault_sigma * model.sigma_time[i]
                cm_step = pyacs.lib.glinalg.corr_to_cov( corr_step,  std_vector )
                inv_cm_step = pyacs.lib.glinalg.syminv( cm_step )
                model.N[ i * model.nfaults:(i+1) * model.nfaults , i * model.nfaults:(i+1) * model.nfaults ] = \
                    model.N[ i * model.nfaults:(i+1) * model.nfaults , i * model.nfaults:(i+1) * model.nfaults ] + inv_cm_step
            
            bar.finish()

    
    return model


