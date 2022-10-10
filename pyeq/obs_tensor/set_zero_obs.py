###############################################################################
def set_zero_obs(T_OBS, method='last_date'):
    ###############################################################################
    """
    Makes observation to start with 0 displacement at the first date.

    :param T_OBS: the observation tensor
    :param method: 'last_date' will shift displacement so that first date is zero. 'median,n','mean,n' will set the 0 as the n-epochs median, mean respectively
    :param verbose: verbose mode

    :return: new observation tensor
    """

    # import
    import numpy as np
    # VERBOSE
    import pyeq.message.message as MESSAGE
    import pyeq.message.verbose_message as VERBOSE
    import pyeq.message.error as ERROR
    import pyeq.message.warning as WARNING
    import pyeq.message.debug_message as DEBUG
    import pyacs.debug
    from icecream import ic



    # get the values
    if method == 'last_date':
        MESSAGE("Setting displacement reference from the first date")
        # loop on sites
        for i in np.arange(T_OBS.shape[1]):
            # search for the first observation
            lindex = np.where(~np.isnan(T_OBS[:, i, 0]))[0]
            index = lindex[0]
            NEU = T_OBS[index, i, 0:3]
            # shift T_OBS
            T_OBS[lindex, i, 0:3] = T_OBS[lindex, i, 0:3] - NEU
    if 'median' in method:
        ndays = 3
        if ',' in method:
            ndays = int(method.split(',')[-1])
        MESSAGE("Setting displacement reference from the median over %d epochs" % ndays)
        # loop on sites
        for i in np.arange(T_OBS.shape[1]):
            # search for the first observation
            lindex = np.where(~np.isnan(T_OBS[:, i, 0]))[0]
            index = lindex[0]
            NEU = np.nanmedian(T_OBS[index:index+ndays, i, 0:3],axis=0)
            # shift T_OBS
            T_OBS[lindex, i, 0:3] = T_OBS[lindex, i, 0:3] - NEU
        if pyacs.debug():
            ic(NEU)
    if 'mean' in method:
        ndays = 3
        if ',' in method:
            ndays = int(method.split(',')[-1])
        MESSAGE("Setting displacement reference from the men over %d epochs" % ndays)
        # loop on sites
        for i in np.arange(T_OBS.shape[1]):
            # search for the first observation
            lindex = np.where(~np.isnan(T_OBS[:, i, 0]))[0]
            index = lindex[0]
            NEU = np.nanmedian(T_OBS[index:index+ndays, i, 0:3],axis=0)
            # shift T_OBS
            T_OBS[lindex, i, 0:3] = T_OBS[lindex, i, 0:3] - NEU
            if pyacs.debug():
                ic(i,NEU)

    return T_OBS