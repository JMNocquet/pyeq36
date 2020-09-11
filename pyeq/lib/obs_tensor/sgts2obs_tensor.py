###############################################################################
def sgts2tensor( sgts, rounding='day' , verbose=False ):
###############################################################################
    """
    returns a numpy array including all gts (.data ENU) information as a 3D tensor and 1D array of site code
    
    :param sgts: Sgts instance
    :param rounding: 'day','hour','second'. all dates will be rounded to the nearest chosen day, hour or second. default is 'day'
    :param verbose: boolean for verbose mode
    
    :return: the numpy 3D array T_OBS, np_names, np_date_s
    
    :note: T_OBS[k,i,0] returns the East value at the k_th date for site i in mm
     
    """

    # import
    import numpy as np
    import pyacs.lib.astrotime as at
    
    # np print option for debug
    #np.set_printoptions(precision=2 , suppress=True)
    
    # np_names
    np_names = np.array( sorted(sgts.lcode() ) )
    
    
    # get all dates
    np_seconds_sol = np.array([],dtype=np.int64)
    for i in np.arange( np_names.shape[0] ):
        # code
        code = np_names[i]
        # get the gts for current site
        wts = sgts.__dict__[code]
        # update np of all dates        
        np_seconds_sol = np.unique( np.sort( np.append( np_seconds_sol, at.decyear2seconds( wts.data[:,0] , rounding=rounding )) ))
        

    # initialize T
    T = np.full( (np_seconds_sol.shape[0] , np_names.shape[0], 6  ) , np.nan )
    
    if verbose:
        print("-- shape of T " , T.shape)
    
    # loop on gts in sgts
    
    for i in np.arange( np_names.shape[0] ):

        # code
        code = np_names[i]
        
        if verbose:
            print("-- %04d / %s " % (i,code) )
        
        # get the gts for current site
        wts = sgts.__dict__[code]

        # date of the current wts in seconds
        np_seconds_site = at.decyear2seconds( wts.data[:,0] , rounding=rounding )
        
        # check for non duplicated dates
        if np.min( np.diff( np_seconds_site) ) <= 0:
            print("!!!ERROR: there is a date problem in gts: %s" % code)
            print("!!!ERROR: most probably, your round option (%s) leads to two successive equal dates" % (rounding))
            print(np.min( np.diff( np_seconds_site) ) )
            print(np.diff( np_seconds_site) )
            print( np_seconds_site )
            
            import sys
            sys.exit()
        
        # current gts date index in T
        
        lindex= np.intersect1d(np_seconds_sol, np_seconds_site, assume_unique=True, return_indices=True)[1]

        if verbose:        
            print("-- number of observation in ts: %s:%d" % (code, wts.data.shape[0]) )
            print("-- number of observation in T : %d" % (T.shape[0]) )
            print("-- number of observation to be filled in T from lindex : %d" % (lindex.shape[0]) )
            print("-- size of T slice to be filled: " , T[lindex,i,:].shape )
            
        # fill T - ENU - SE SN SU
        T[lindex,i,0] = wts.data[:,2]*1000.
        T[lindex,i,1] = wts.data[:,1]*1000.
        T[lindex,i,2] = wts.data[:,3]*1000.

        T[lindex,i,3] = wts.data[:,5]*1000.
        T[lindex,i,4] = wts.data[:,4]*1000.
        T[lindex,i,5] = wts.data[:,6]*1000.
        
        # return

        
        
    return T , np_names , np_seconds_sol

