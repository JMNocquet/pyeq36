def make_model_covariance_01( model ):
    """
    Create a model covariance matrix Cm
    
    cov_01 assumes:
    - no temporal smoothing: all correlation coefficients relating 2 parameters at different time are null
    - model.np_sigma is a scalar controlling the a priori average slip rate per day ( sigma**2 in mm**2 /day )
    - model.np_sigma is rescaled according to the model step duration
    
    - with this formulation, np_sigma can also be a vector of length model.nfaults
    - finally model.np_sigma can be a matrix whose columns are sigma vectors for each model time step
    """
    
    import pyeq.lib.regularization
    import numpy as np
    import pyacs.lib.glinalg
    import scipy.linalg
    
    # first get the correlation matrix
    corr_step = pyeq.lib.regularization.make_model_spatial_correlation_matrix( model.Dm , model.dc, cm_type =  model.cm_type )

    # computes normalization
    if model.cm_norm == 'd0':
        nf = pyeq.lib.regularization.normalization_factor( model.SGEOMETRY , model.dc )
    else:
        nf = 1.
    
    print("-- normalization factor for model step covariance matrix: %.3lf " % nf )
    
    # computes the renormlization factor according the model step duration
    
    delta_time_in_days = np.diff(model.np_model_date_s / 86400.)

    # compute the inversion correlation matrix
    inv_corr_step = pyacs.lib.glinalg.syminv( corr_step )
    
    # np_sigma is a scalar or a vector
    # normalization: Cm would be multiplied by delta_time_in_days * nf * sigma**2
    # inv_Cm is therefore mutliplied by its inverse
    # for each time step and store in a tuple of 2D-array
    # which is stored in a tuple to be then used by scipy.linalg.block_diag
    
    if model.np_sigma.ndim in [0,1]:
        
        normalization_vector = 1./ ( delta_time_in_days * nf * model.np_sigma**2 )
        bb = np.expand_dims(np.expand_dims( normalization_vector,axis=0), axis=0).T
        N = tuple( inv_corr_step * bb )
        
        model.INV_CM = scipy.linalg.block_diag( *N )

    # np_sigma is a matrix whose columns are sigma vectors for each model time step

    if model.np_sigma.ndim ==2:

        SIGMA = []
        for i in np.arange( model.np_sigma.shape[1] ):
            sigma = 1. / np.outer( model.np_sigma[:,i], model.np_sigma[:,i] ) * 1./( delta_time_in_days * nf )
            SIGMA.append( sigma )
        SIGMA = np.array(SIGMA)
        N = inv_corr_step * np.array(SIGMA)
        model.INV_CM = scipy.linalg.block_diag( *N )
        

    return model
