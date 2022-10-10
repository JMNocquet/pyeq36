def make_model_spatial_correlation_matrix( Dm , Dc, cm_type = 'exponential'):
    """
    From a matrix Dm_ij including the distances between pair of model parameters, 
    calculates the correlation matrix
    
    
    :param Dm: 2D numpy-array including the distances
    :param Dc: critical distance
    :param cm_type: choose among ['exponential','r_exp','cos_sin_exp','gaussian','inv_ch']
    
    :note: kernels taken from B. Valette, Inverse Problem Lectures.
    """
    
    # import
    import numpy as np
    
    r = Dm / Dc
    r2 = np.sqrt(2)
    
    if cm_type == 'exponential':
        Cm_corr = np.exp( -r )

    if cm_type == 'r_exp':
        Cm_corr = ( 1 + r ) * np.exp( -r )
    
    if cm_type == 'cos_sin_exp':
        Cm_corr = ( np.cos(r/r2) + np.sin(r/r2) ) * np.exp( -r/r2 )
        
    if cm_type == 'gaussian':
        Cm_corr = np.exp( -r**2 )
    
    if cm_type == 'inv_ch':
        Cm_corr = 1. / np.cosh( -r )
        
    return Cm_corr