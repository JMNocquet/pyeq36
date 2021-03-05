###############################################################################
def reorder_obs_tensor_from_names( T_OBS , np_names_t_obs, new_np_names, verbose=False ):
###############################################################################
    """
    reorder an observation tensor according to names
    
    :param T_OBS: the observation tensor
    :param np_names_t_obs: 1D array of names for the input observation tensor
    :param new_np_names: 1D array of names for the output observation tensor
    :param verbose: verbose mode
    
    :return: new observation tensor
    
    :note 1: names provided in new_np_names not present in np_names raise an Error
    :note 2: data for names present in np_names but not provided in new_np_names are deleted
    """
    
    # import
    import numpy as np

    
    np_index = np.array([], dtype = int )
    
    
    # loop on sites
    for i in np.arange( new_np_names.shape[0] ):
        
        # check all requested names are in np_names
        if new_np_names[i] not in np_names_t_obs:
            print("!!!ERROR: site %s requested but not present in observation. Exit." % new_np_names[i] )
            import sys
            sys.exit()
            
        else:
            # build the list of name index
            np_index = np.append(np_index , np.where(  np_names_t_obs == new_np_names[i] ) )
    
    # just to be sure
    np_index = np_index.flatten()
    
    return T_OBS[ : , np_index , : ]
    
