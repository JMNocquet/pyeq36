###############################################################################
def dk_from_obs_tensor( k , T_OBS , component, normalized=True , verbose=False ):
###############################################################################
    """
    From an observation tensor, returns the k_th observation as a vector
    
    :param k: index of the requested date
    :param T_OBS: the observation tensor
    :param component: numpy array for the requested components. np.array([0,1]) for East and North, np.array([0,1,2]) for the 3 components
    :param normalized: boolean. If True, value are normalized by their uncertainty.
    :param verbose: verbose mode
    """
    
    # import
    import numpy as np
    
    # removes site with no data
    lindex = np.where( T_OBS[k,:,0] == np.nan )

    # extract the observation at date index k for site having data
    T_OBS_K = T_OBS[k,lindex,:]
    
    # normalize
    if normalized:
        T_OBS_K[:,0:3] = T_OBS_K[:,0:3] / T_OBS_K[:,3:7]
    
    # convert to 2D numpy array
    
    return T_OBS_K[:,0:3].reshape(-1,component.shape[0])
    
    
