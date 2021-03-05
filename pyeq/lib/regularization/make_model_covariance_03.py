###############################################################################
def make_model_covariance_03( model ):
###############################################################################
    """
    Create a model covariance matrix Cm
    
    cov_03 is an update of cov_02 to handle offset estimates at the origin time 
    
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
    import pyeq.log
    from progress.bar import Bar

    ###########################################################################
    # INITIALIZE SOME PARAMETERS
    ###########################################################################

    # size of the block matrix associated with slip parameters
    ncm = model.nfaults * model.nstep

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

###############################################################################
# SAVE DIAGONAL BLOCKS BEFORE THEY GET MODIFIED BY THE SPATIAL REGULARIZATION
###############################################################################

    if model.tau !=0:
        print("-- saving Observation Normal diagonal for later use")
        DIAGONAL_N_OBS = np.copy( np.diag( model.N )[ :ncm ] )
        # since temporal smoothing is applied independently we need to rescale spatial and temporal variance 
        nsigma_tau = 1./np.sqrt(2.)
    else:
        nsigma_tau=1.

###############################################################################
# CASE NO TEMPORAL SMOOTHING
###############################################################################

    # compute the inversion correlation matrix
    # !!! THERE IS A SURPRISE HERE
    # THE COVARIANCE MATRIX APPEARS TO BE BETTER CONDITIONED
    # IF THE COVARIANCE IS INVERSED COMPARED TO THE CORRELATION
    # MATRIX INVERSED AND THEN MULTIPLIED BY 1/sigma**2
    # THIS LEADS TO SUB-OPTIMAL CODING
    print("-- Computing the individual model_step inverse spatial correlation matrix (inv_corr_step)")

    #######################################################################        
    if model.sigma_type == 'time_constant/fault_constant':
    #######################################################################        
        print("-- Computing the time step covariance matrix - sigma_type is %s " % ( model.sigma_type ) )

        # change on 02/11/2020
        # single model step case
        if delta_time_in_days.size == 1:
            boolean_constant_time_step = True

        # multiple model step case
        else:
            if np.max( np.fabs( np.diff( delta_time_in_days ) ) ) == 0:
                boolean_constant_time_step = True
            else:
                boolean_constant_time_step = False


        if boolean_constant_time_step:

            print("-- Constant time step case ")
            # change 16/04/2020
            # since parameters are expressed in mm/day (slip rate) and green's function have been
            # rescaled by the time step duration in the observation normal matrix,
            # regularization only applies on rates and do not need delta_day
            # normalization_factor = delta_time_in_days[0] * (model.sigma_time/nsigma_tau)**2 * nf
            normalization_factor = (model.sigma_time/nsigma_tau)**2 * nf
            cm_step = corr_step * normalization_factor
            inv_cm_step = pyacs.lib.glinalg.syminv( cm_step )
            print("-- Filling diagonal blocks")
            model.N[ :ncm , :ncm ] = model.N[ :ncm , :ncm ] + scipy.linalg.block_diag(*([inv_cm_step] * model.nstep))
        else:
            print("-- Variable time step case")
            # THIS COULD BE IMPROVED: HERE THIS REQUIRES NSTEP INVERSION
            # ALTHOUGH MANY OF THEM ARE PROBABLY THE SAME

            bar = Bar('', max= model.nstep , suffix = '%(percent).1f%% - %(eta)ds')

            for i in np.arange( model.nstep ):

                bar.next()
                
                # change 16/04/2020
                # since parameters are expressed in mm/day (slip rate) and green's function have been
                # rescaled, regularization only applies on rates and do not need delta_day
                #normalization_factor = delta_time_in_days[i] * (model.sigma_time/nsigma_tau)**2 * nf
                normalization_factor = (model.sigma_time/nsigma_tau)**2 * nf
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
#            normalization_factor = delta_time_in_days[i] * (model.sigma_time[i]/nsigma_tau)**2 * nf
            normalization_factor = (model.sigma_time[i]/nsigma_tau)**2 * nf
            cm_step = corr_step * normalization_factor
            inv_cm_step = pyacs.lib.glinalg.syminv( cm_step )
            model.N[ i * model.nfaults:(i+1) * model.nfaults , i * model.nfaults:(i+1) * model.nfaults ] = \
                model.N[ i * model.nfaults:(i+1) * model.nfaults , i * model.nfaults:(i+1) * model.nfaults ] + inv_cm_step
        
        bar.finish()

    #######################################################################        
    if model.sigma_type == 'time_constant/fault_variable':
    #######################################################################        
        print("-- Computing the time step covariance matrix - sigma_type is %s " % ( model.sigma_type ) )


        # change on 02/11/2020
        # single model step case
        if delta_time_in_days.size == 1:
            boolean_constant_time_step = True

        # multiple model step case
        else:
            if np.max(np.fabs(np.diff(delta_time_in_days))) == 0:
                boolean_constant_time_step = True
            else:
                boolean_constant_time_step = False

        if boolean_constant_time_step:
            print("-- Constant time step case ")
#            std_vector = np.sqrt( delta_time_in_days[0]) * model.time_sigma/nsigma_tau * model.sigma_fault * nf
            std_vector = model.time_sigma/nsigma_tau * model.sigma_fault * nf
            inv_cm_step = pyacs.lib.glinalg.corr_to_cov( corr_step,  std_vector )
            model.N[ :ncm , :ncm ] = model.N[ :ncm , :ncm ] + scipy.linalg.block_diag(*([inv_cm_step] * model.nstep))
        else:
            print("-- Variable time step case ")
            # THIS COULD BE IMPROVED: HERE THIS REQUIRES NSTEP INVERSION
            # ALTHOUGH MANY OF THEM ARE PROBABLY THE SAME

            bar = Bar('', max= model.nstep , suffix = '%(percent).1f%% - %(eta)ds')

            for i in np.arange( model.nstep ):

                bar.next()
#                std_vector = np.sqrt( delta_time_in_days[i] ) * model.fault_sigma/nsigma_tau * model.sigma_time/nsigma_tau * nf
                std_vector =  model.fault_sigma/nsigma_tau * model.sigma_time/nsigma_tau * nf
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
#            std_vector = np.sqrt( delta_time_in_days[i] ) * model.fault_sigma * model.sigma_time[i]/nsigma_tau
            std_vector = model.fault_sigma * model.sigma_time[i]/nsigma_tau
            cm_step = pyacs.lib.glinalg.corr_to_cov( corr_step,  std_vector )
            inv_cm_step = pyacs.lib.glinalg.syminv( cm_step )
            model.N[ i * model.nfaults:(i+1) * model.nfaults , i * model.nfaults:(i+1) * model.nfaults ] = \
                model.N[ i * model.nfaults:(i+1) * model.nfaults , i * model.nfaults:(i+1) * model.nfaults ] + inv_cm_step
        
        bar.finish()

###############################################################################
# CASE TEMPORAL SMOOTHING
###############################################################################

    if model.tau != 0:
        print("-- Building the model covariance matrix CM with temporal smoothing. Correlation time is tau = %.2lf days" % model.tau )
        
        # compute the time distance among model time step matrix
        print("-- Computing the model step time distance matrix" )
        delta_day_2D = pyeq.lib.regularization.time_to_2D_diff_time_in_days(  (model.np_model_date_s / 86400.)[:-1] + delta_time_in_days / 2 )
    
        # temporal correlation at all time step
        print("-- Computing the model step time correlation matrix" )
        CORR_TIME = np.exp( -delta_day_2D / model.tau )

        #######################################################################        
        #if model.sigma_type == 'time_constant/fault_constant':
        #######################################################################        
        print("-- Computing the time step covariance matrix - sigma_type is %s " % ( model.sigma_type ) )

        bar = Bar('', max= model.nstep , suffix = '%(percent).1f%% - %(eta)ds')

        DIAGONAL_REG = np.diag( model.N )[: model.nstep * model.nfaults ] - DIAGONAL_N_OBS

        # fill the element at their proper location - INV_COV_TIME is tri-diagonal
        for i in np.arange( model.nfaults ):
            bar.next()
            STD = 1./np.sqrt( DIAGONAL_REG[ i:model.nstep * model.nfaults:model.nfaults ] )
            COV_TIME = CORR_TIME * np.outer( STD, STD)
            INV_COV_TIME = pyacs.lib.glinalg.syminv( COV_TIME )
            
            for j in np.arange( model.nstep ):
                for k in np.arange( j , model.nstep ):
                    model.N[ j*model.nfaults + i , k*model.nfaults + i ] += INV_COV_TIME[ j , k ]
                    model.N[ k*model.nfaults + i , j*model.nfaults + i ] += INV_COV_TIME[ j , k ]
                
        bar.finish()

###############################################################################
# HANDLE ORIGIN TIME CONSTANTS
###############################################################################
    if model.offset > 0:
        print("-- Adding offset origin time constraint (mm): %.1f " % model.offset )
        # we need to augment model.N with rows and columns corresponding to the offset parameters
        # fill the block in model.N corresponding to the offset paramaters

        np.fill_diagonal( model.N[ -model.nconstant: , -model.nconstant: ] , np.diag(model.N[ -model.nconstant: , -model.nconstant: ])  + 1./model.offset**2 ) 
    
    return model


