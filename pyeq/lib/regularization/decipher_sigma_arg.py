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

    ###########################################################################
    # CHECK ALL REGULARIZATION PARAMETERS DO EXIST AND ARE OF CORRECT TYPE
    ###########################################################################

    # reference date is 1980.0
    ref_date_time = datetime(1980, 1, 1, 0, 0, 0)

    # common
    larg_common = ['regularization','sigma','dc','tau']
    
    for arg in larg_common:
        if arg not in model.__dict__:
            print("!!!ERROR: regularization argument %s missing. Exit." % ( arg ) )
            sys.exit()

    try:
        model.dc = float( model.dc )
    except:
        print("!!!ERROR: dc arg is not a float number. dc=%s. Exit." % ( model.dc ) )
        sys.exit()
        
    try:
        model.tau = float( model.tau )
    except:
        print("!!!ERROR: tau arg is not a float number. tau=%s. Exit." % ( model.tau ) )
        sys.exit()
    
    # check regularization argument
    if not ( model.regularization in ['laplacian','covariance']):
        print("!!!ERROR: regularization must be either 'laplacian' or 'covariance'. regularization argument is %s. Exit." % ( model.regularization ) )
        sys.exit()
    
    # check laplacian argument
    if model.regularization == 'laplacian':
            
        larg_laplacian = ['sigma_smoothing']
        
        for arg in larg_laplacian:
            if arg not in model.__dict__:
                print("!!!ERROR: regularization laplacian argument %s missing. Exit." % ( arg ) )
                sys.exit()

        try:
            model.sigma_smoothing = float( model.sigma_smoothing )
        except:
            print("!!!ERROR: sigma_smoothing arg is not a float number. sigma_smoothing=%s. Exit." % ( model.sigma_smoothing ) )
            sys.exit()


    # check covariance arguments
    if model.regularization == 'covariance':

        larg_covariance = ['cm_type','cm_norm']
        
        for arg in larg_covariance:
            if arg not in model.__dict__:
                print("!!!ERROR: regularization covariance argument %s missing. Exit." % ( arg ) )
                sys.exit()
    
        if model.cm_type not in ['exponential','r_exp','cos_sin_exp','gaussian','inv_ch']:
            print("!!!ERROR: regularization covariance argument cm_type not understood. cmtype=%s missing. Exit." % ( model.cm_type ) )
            sys.exit()
        
    
    ###########################################################################
    # DECIPHER SIGMA AND GENERATE NP_SIGMA
    ###########################################################################
    try:
        float( model.sigma )
        sigma_type = 'constant'
    except:
        sigma_type = 'variable'

    print("-- sigma type is: %s" % sigma_type )
    
    # sigma is constant
    ###########################################################################
    
    if sigma_type == 'constant':
        model.sigma_time = float( model.sigma )
        model.sigma_fault = float( model.sigma )
        model.sigma_type = 'time_constant/fault_constant'
        return model
    
    # sigma is variable
    ###########################################################################
    if sigma_type == 'variable':
    
        [sigma_time,sigma_fault] = model.sigma.replace('[','').replace(']','').split(',')
        
        try:
            model.sigma_time = float( sigma_time )
            model.sigma_type = 'time_constant'
            print("-- sigma is time constant = %.2f" % ( model.sigma_time ) )
        except:
            print("-- sigma is time variable. Getting sigma_time from file: %s" % ( sigma_time ) )
            model.sigma_type = 'time_variable'
            ###################################################################
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
                parse = lambda x: pandas.datetime.strptime(x, '%Y-%m-%d %H:%M:%S') 
                # get np_datetime
                np_datetime = np.array(list(map(parse, d)))
                # convert to seconds
                np_date_s =  np.array(list(map(int,[x.total_seconds() for x in ( np_datetime - ref_date_time) ])), dtype=np.int64)
                
                # load file for sigma
                np_sigma_s = np.genfromtxt( sigma_time , usecols=(2)) 

            
            else:
                # not a file
                print("!!!ERROR sigma (sigma_time) option error. %s is not a file" % sigma_time ) 
                sys.exit()
        
            # make interpolation for each model time step
            np_mid_model_date_s = ( model.np_model_date_s[:-1] +  np.diff( model.np_model_date_s)/2. ).astype( int )
            f = interpolate.interp1d( np_date_s , np_sigma_s )
            model.sigma_time = f( np_mid_model_date_s )   # use interpolation function returned by `interp1d`

            print("-- sigma_time range: %.2f - %.2f" % ( np.min(model.sigma_time) , np.max(model.sigma_time) ) )

        try:
            model.sigma_fault = float( sigma_fault )
            model.sigma_type = model.sigma_type +'/fault_constant'
            print("-- sigma is fault constant = %.2f" % ( model.sigma_fault ) )
        except:
            print("-- sigma is fault variable. Getting sigma_fault from file: %s" % ( sigma_fault ) )
            try:
                model.sigma_fault = np.genfromtxt( sigma_fault )
                model.sigma_type = model.sigma_type +'/fault_varaiable'
            except:
                print("!!!ERROR sigma (sigma_fault) option error. %s is not a file" % sigma_fault ) 
                sys.exit()
                
    return model    
            
        
        
    
    
    
    
    