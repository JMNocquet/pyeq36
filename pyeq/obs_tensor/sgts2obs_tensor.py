###############################################################################
def sgts2tensor( sgts, rounding='day' , verbose=False , ncpu=1 ):
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

    import pyeq.message.message as MESSAGE
    import pyeq.message.verbose_message as VERBOSE
    import pyeq.message.error as ERROR
    import pyeq.message.debug_message as DEBUG

    from tqdm import tqdm

    # Check data
    for code in sorted(sgts.lcode()):
        try:
            if sgts.__dict__[code].data.shape[0] <2 :
                ERROR(("%s has less than two epochs." % code),exit=True)
        except:
            ERROR(("%s has problem with .data." % code), exit=True)


    MESSAGE("Converting Sgts with %d time series into Observation tensor" % sgts.n())

    # np print option for debug
    #np.set_printoptions(precision=2 , suppress=True)
    
    # np_names
    np_names = np.array( sorted(sgts.lcode() ) )
    
    ###########################################################################
    # GET ALL DATES
    ###########################################################################

    # if rounding is day then we just get the min-max date and iterate a 1 days sampling

    if rounding == 'day':
        lstart_decyear = []
        lend_decyear = []
        for code in sorted(sgts.lcode()):
                lstart_decyear.append( sgts.__dict__[code].data[0, 0] )
                lend_decyear.append(sgts.__dict__[code].data[-1, 0])
        start_decyear = np.min( np.array( lstart_decyear ) )
        end_decyear   = np.max( np.array( lend_decyear ) )
        start_seconds = at.decyear2seconds( start_decyear, rounding='day' )
        end_seconds = at.decyear2seconds( end_decyear, rounding='day' )
        np_seconds_sol = np.arange( start_seconds, end_seconds+1, 86400, dtype=np.int64 )

    else:

        # get all present dates from the time series
        np_seconds_sol = np.array([],dtype=np.int64)
        MESSAGE("Collecting dates in all time series")

        for i in tqdm(np.arange( np_names.shape[0] )):
            # code
            code = np_names[i]
            # get the gts for current site
            wts = sgts.__dict__[code]
            # update np of all dates
            np_seconds_sol = np.unique( np.sort( np.append( np_seconds_sol, at.decyear2seconds( wts.data[:,0] , rounding=rounding )) ))

    # initialize T
    T = np.full( (np_seconds_sol.shape[0] , np_names.shape[0], 6  ) , np.nan )
    
    DEBUG(("shape of T: %s " % (T.shape,)))
    
    # loop on gts in sgts
    if ncpu == 1:
        for i in tqdm( np.arange( np_names.shape[0] ) ):

            # code
            code = np_names[i]

            DEBUG("%04d / %s " % (i,code) )

            # get the gts for current site
            wts = sgts.__dict__[code]

            # date of the current wts in seconds
            np_seconds_site = at.decyear2seconds( wts.data[:,0] , rounding=rounding )

            # check for non duplicated dates
            try:
                if np.min( np.diff( np_seconds_site) ) <= 0:
                    ERROR("there is a date problem in gts: %s" % code)
                    ERROR("most probably, your round option (%s) leads to two successive equal dates" % (rounding))
                    ERROR(("%d" % np.min( np.diff( np_seconds_site) ) ))
                    print(np.diff( np_seconds_site) )
                    print( np_seconds_site )
                    ERROR("",exit=True)
            except:
                ERROR(("Problem with %s" % code ),exit=True)

            # current gts date index in T

            lindex= np.intersect1d(np_seconds_sol, np_seconds_site, assume_unique=True, return_indices=True)[1]

            DEBUG("number of observation in ts: %s:%d" % (code, wts.data.shape[0]) )
            DEBUG("number of observation in T : %d" % (T.shape[0]) )
            DEBUG("number of observation to be filled in T from lindex : %d" % (lindex.shape[0]) )
            DEBUG("size of T slice to be filled: %s" % (T[lindex,i,:].shape,) )

            # fill T - ENU - SE SN SU
            T[lindex,i,0] = wts.data[:,2]*1000.
            T[lindex,i,1] = wts.data[:,1]*1000.
            T[lindex,i,2] = wts.data[:,3]*1000.

            T[lindex,i,3] = wts.data[:,5]*1000.
            T[lindex,i,4] = wts.data[:,4]*1000.
            T[lindex,i,5] = wts.data[:,6]*1000.

    ## PARALLEL PROCESSING ##

    if ncpu > 1:
        MESSAGE("Parallel processing using %d CPU for %d time series" % (ncpu,np_names.shape[0]))

        import ipyparallel as ipp

        lgts = []
        for i in np.arange(np_names.shape[0]):
            code = np_names[i]
            lgts.append(sgts.__dict__[code].data)

        def task( data, np_s):
            import pyacs.lib.astrotime as at
            import numpy as np

            # date of the current wts in seconds
            np_seconds_site = at.decyear2seconds(data[:, 0], rounding=rounding)

            # check for non duplicated dates
            if np.min(np.diff(np_seconds_site)) <= 0:
                return None
            # current gts date index in T
            lidx = np.intersect1d(np_s, np_seconds_site, assume_unique=True, return_indices=True)[1]

            # fill T - ENU - SE SN SU
            TT = np.full((np_s.shape[0], 6), np.nan)
            TT[lidx, 0] = data[:, 2] * 1000.
            TT[lidx, 1] = data[:, 1] * 1000.
            TT[lidx, 2] = data[:, 3] * 1000.

            TT[lidx, 3] = data[:, 5] * 1000.
            TT[lidx, 4] = data[:, 4] * 1000.
            TT[lidx, 5] = data[:, 6] * 1000.

            return TT

        # request a cluster
        with ipp.Cluster(n=ncpu) as rc:
            # get a view on the cluster
            view = rc.load_balanced_view()
            # submit the tasks
            asyncresult = view.map_async(task,  lgts, [np_seconds_sol]*len(lgts) )
            # wait interactively for results
            asyncresult.wait_interactive()
            # retrieve actual results
            result = asyncresult.get()
        # at this point, the cluster processes have been shutdown

        # populates with results
        for i in np.arange(np_names.shape[0]):
            T[:, i, :] = result[i]

    # return
    return T, np_names, np_seconds_sol

