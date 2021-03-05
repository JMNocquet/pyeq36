def decipher_sigma_arg( model ):
    """
    decipher and check sigma argument
    """

    ###########################################################################
    # IMPORT
    ###########################################################################

    import numpy as np
    from scipy import interpolate
    import sys , os
    import pandas
    from datetime import datetime

    import pyeq.message.message as MESSAGE
    import pyeq.message.verbose_message as VERBOSE
    import pyeq.message.error as ERROR
    import pyeq.message.warning as WARNING
    import pyeq.message.debug_message as DEBUG

    ###########################################################################
    # CHECK ALL REGULARIZATION PARAMETERS DO EXIST AND ARE OF CORRECT TYPE
    ###########################################################################

    # reference date is 1980.0
    ref_date_time = datetime(1980, 1, 1, 0, 0, 0)

    # regularization
    if 'regularization' not in model.__dict__:
        ERROR("regularization argument not provided.",exit=True)

    lregularization = ['laplacian','covariance']
    if model.regularization not in lregularization:
        ERROR(("regularization argument not understood. regularieation should be in %s" % lregularization),exit=True)

    # covariance type
    if 'covariance' in model.regularization:
        larg_common = ['sigma','dc','tau','cm_type']

        for arg in larg_common:
            if arg not in model.__dict__:
                ERROR(("regularization argument %s missing. Exit." % ( arg ) ),exit=True)

        try:
            model.dc = float( model.dc )
        except:
            ERROR(("dc arg is not a float number. dc=%s. Exit." % ( model.dc ) ),exit=True)

        try:
            model.tau = float( model.tau )
        except:
            ERROR(("ERROR: tau arg is not a float number. tau=%s. Exit." % ( model.tau ) ), exit=True)
            sys.exit()

        if model.cm_type not in ['exponential', 'r_exp', 'cos_sin_exp', 'gaussian', 'inv_ch']:
            ERROR(("regularization covariance argument cm_type not understood. cmtype=%s missing. Exit." % (model.cm_type)),exit=True)
            sys.exit()

    # laplacian
    if model.regularization == 'laplacian':
            
        larg_laplacian = ['lambda_temporal_smoothing','lambda_spatial_smoothing']
        
        for arg in larg_laplacian:
            if arg not in model.__dict__:
                ERROR(("regularization argument %s missing. Exit." % ( arg ) ),exit=True)

        try:
            model.sigma_spatial_smoothing = float( model.lambda_spatial_smoothing )
        except:
            ERROR(("lambda_spatial_smoothing arg is not a float number: %s" % ( model.lambda_spatial_smoothing ) ), exit=True)

        try:
            model.sigma_temporal_smoothing = float( model.lambda_temporal_smoothing )
        except:
            ERROR(("lambda_temporal_smoothing arg is not a float number: %s" % ( model.lambda_temporal_smoothing ) ), exit=True)

    ###########################################################################
    # DECIPHER SIGMA AND GENERATE NP_SIGMA
    ###########################################################################
    try:
        float( model.sigma )
        sigma_type = 'constant'
    except:
        sigma_type = 'variable'

    DEBUG("sigma type is: %s" % sigma_type )
    
    # sigma is float = constant in space and time
    ###########################################################################
    
    if sigma_type == 'constant':
        model.sigma_time = float( model.sigma )
        model.sigma_fault = float( model.sigma )
        model.sigma_type = 'time_constant/fault_constant'
        model.np_sigma = np.array([float( model.sigma )] * (model.nfaults*model.nstep) )
        return model
    
    # sigma is variable
    ###########################################################################
    if sigma_type == 'variable':
    
        [sigma_time,sigma_fault] = model.sigma.replace('[','').replace(']','').split(',')
        
        try:
            model.sigma_time = float( sigma_time )
            model.sigma_type = 'time_constant'
            VERBOSE("sigma is time constant = %.2f" % ( model.sigma_time ) )
        except:
            VERBOSE("sigma is time variable. Getting sigma_time from file: %s" % ( sigma_time ) )
            model.sigma_type = 'time_variable'

            ###################################################################
            # WE EXPECT SIGMA TO CHANGE THROUGH TIME
            # EXPECTED FORMAT DATE - SIGMA
            #2016-04-17 12:00:00 100.
            #2016-05-18 12:00:00  90.
            #2016-06-13 12:00:00   8.
            # ....
            ###################################################################

            if os.path.isfile( sigma_time ) :
            
                # load file for dates
                d = np.genfromtxt( sigma_time , dtype=str , usecols=(0,1))
                # join the two columns
                l = []
                for i in np.arange( d.shape[0] ):
                     l.append( ("%s %s" % ( d[i,0],d[i,1] ) ) )
                d = np.array( l )
                # parser
                # the line below raises a python fure warning
                #parse = lambda x: pandas.datetime.strptime(x, '%Y-%m-%d %H:%M:%S')
                # new line not tested 11/01/2021
                import datetime
                parse = lambda x: datetime.datetime.strptime(x, '%Y-%m-%d %H:%M:%S')
                # get np_datetime
                np_datetime = np.array(list(map(parse, d)))
                # convert to seconds
                np_date_s =  np.array(list(map(int,[x.total_seconds() for x in ( np_datetime - ref_date_time) ])), dtype=np.int64)
                
                # load file for sigma
                np_sigma_s = np.genfromtxt( sigma_time , usecols=(2)) 

            
            else:
                # not a file
                ERROR(("sigma (sigma_time) option error. %s is not a file" % sigma_time ),exit=True)
                sys.exit()
        
            # make interpolation for each model time step
            np_mid_model_date_s = ( model.np_model_date_s[:-1] +  np.diff( model.np_model_date_s)/2. ).astype( int )
            f = interpolate.interp1d( np_date_s , np_sigma_s )
            model.sigma_time = f( np_mid_model_date_s )   # use interpolation function returned by `interp1d`

            VERBOSE("sigma_time range: %.2f - %.2f" % ( np.min(model.sigma_time) , np.max(model.sigma_time) ) )
            model.np_sigma = (np.ones(model.nfaults) * model.sigma_time.reshape(model.nstep,-1)).flatten()
            model.sigma_fault = float( sigma_fault )
            model.sigma_type = model.sigma_type +'/fault_constant'
            VERBOSE("sigma is fault constant = %.2f" % ( model.sigma_fault ) )
            return model

        # sigma is fault variable
        VERBOSE("sigma is fault variable. Getting sigma_fault from file: %s" % ( sigma_fault ) )
        try:
            model.sigma_fault = np.genfromtxt( sigma_fault )
            model.sigma_type = model.sigma_type +'/fault_variable'
            model.np_sigma = np.repeat([model.sigma_fault],model.nstep,axis=0).flatten()
        except:
            ERROR(("sigma (sigma_fault) option error. %s is not a file" % sigma_fault ),exit=True)

    if model.np_sigma.ndim != 1:
        ERROR("in decipher_sigma_arg. np_sigma is not a 1D numpy array.")

    if model.np_sigma.shape[0] != model.nfaults*model.nstep:
        ERROR("in decipher_sigma_arg. np_sigma has a wrong dimension:%d. nfaults * nstep: %d * %d = %d" %
              (model.np_sigma.shape[0],model.nfaults,model.nstep,model.nfaults*model.nstep))

    return model    
            
        
        
    
    
    
    
    