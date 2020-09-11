###############################################################################
def set_zero_at_first_obs( T_OBS , verbose=False ):
###############################################################################
    """
    Makes observation to start with 0.
    
    :param T_OBS: the observation tensor
    :param verbose: verbose mode
    
    :return: new observation tensor
    """
    
    # import
    import numpy as np

    # loop on sites
    
    for i in np.arange( T_OBS.shape[1] ):
        # search for the first observation
        lindex = np.where( ~np.isnan(T_OBS[:,i,0]) )[0]
        index = lindex[0]

        # get the values
        NEU = T_OBS[index,i,0:3]

        # removes the values
        T_OBS[lindex,i,0:3] = T_OBS[lindex,i,0:3] - NEU
    
    return T_OBS 