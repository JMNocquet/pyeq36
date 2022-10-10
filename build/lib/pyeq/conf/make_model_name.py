def make_model_name( model ):
    """
    produces output name according to input parameters
    current convention is

    """
    # import
    
    import pyacs.lib.astrotime as at
    import numpy as np


    # start_date
    model_datetime = at.seconds2datetime( model.np_model_date_s[0] )
    str_yr = str(model_datetime.year)[-2:]
    doy = at.datetime2dayno( model_datetime)[1]
    duration = int( (model.np_model_date_s[-1] - model.np_model_date_s[-0] ) /86400. )

    str_date = ("%s%03d_%d" % ( str_yr , doy , int( duration ) ) )


    # regularization covariance case
    if model.regularization == 'covariance':
        str_reg = 'cov'

        if model.sigma_type == 'time_constant/fault_constant':
            if model.sigma_time < 1:
                tmp_str='{0:1.0E}'.format(model.sigma_time)
                number = int(tmp_str.split('E')[0])
                exponent = int(tmp_str.split('-')[1])
                str_sigma = ("%dE%d" % ( number , exponent ) )
            else:
                str_sigma = ("%d" % ( int( model.sigma_time ) ) )
        else:
            if 'time_constant' in model.sigma_type:
                if model.sigma_time < 1:
                    str_sigma = ("%dE1v" % ( int( model.sigma_time * 10 ) ) )
                else:
                    str_sigma = ("%dv" % ( int( model.sigma_time ) ) )
            else:
                if 'fault_constant' in model.sigma_type:
                    str_sigma = ("vc" )
                else:
                    str_sigma = ("vv" )

        # dc
        str_dc = ("%d" % ( int(model.dc) ) )

        # tau
        str_tau = ("%d" % (int(model.tau)))

        # final name
        name = ("%s_%s_%s_%s_%s_%d_%d" % ( str_reg , str_sigma , str_dc , str_tau , str_date , model.nstep, model.geometry.shape[0] ) )

    # regularization laplacian_like case

    if model.regularization == 'laplacian_like':
        str_reg = 'll'
    if model.regularization == 'cvxopt':
        str_reg = 'cvx'
    if model.regularization == 'laplacian':
        str_reg = 'lp'

    if model.regularization in[ 'laplacian_like' ,'laplacian', 'laplace','cvxopt']:

        # lambda_spatial_smoothing
        if isinstance(model.lambda_spatial_smoothing,np.ndarray):
            str_lss = ("vs")
        else:
            if (model.lambda_spatial_smoothing < 1) and (model.lambda_spatial_smoothing !=0) :
                tmp_str='{0:1.0E}'.format(model.lambda_spatial_smoothing)
                number = int(tmp_str.split('E')[0])
                exponent = int(tmp_str.split('-')[1])
                str_lss = ("%dE%d" % ( number , exponent ) )
            else:
                str_lss = ("%d" % ( int( model.lambda_spatial_smoothing ) ) )

        # lambda_temporal_smoothing
        if (model.lambda_temporal_smoothing) < 1 and (model.lambda_temporal_smoothing !=0):
            tmp_str='{0:1.0E}'.format(model.lambda_temporal_smoothing)
            number = int(tmp_str.split('E')[0])
            exponent = int(tmp_str.split('-')[1])
            str_lts = ("%dE%d" % ( number , exponent ) )
        else:
            str_lts = ("%d" % ( int( model.lambda_temporal_smoothing ) ) )

        if model.sigma_type == 'time_constant/fault_constant':
            if model.sigma_time == 0:
                str_sigma ='0'
            else:
                if model.sigma_time < 1:
                    tmp_str='{0:1.0E}'.format(model.sigma_time)
                    number = int(tmp_str.split('E')[0])
                    exponent = int(tmp_str.split('-')[1])
                    str_sigma = ("%dE%d" % ( number , exponent ) )
                else:
                    str_sigma = ("%d" % ( int( model.sigma_time ) ) )
        else:
            if 'time_constant' in model.sigma_type:
                if model.sigma_time < 1:
                    str_sigma = ("%dE1v" % ( int( model.sigma_time * 10 ) ) )
                else:
                    str_sigma = ("%dv" % ( int( model.sigma_time ) ) )
            else:
                if 'fault_constant' in model.sigma_type:
                    str_sigma = ("vc" )
                else:
                    str_sigma = ("vv" )


        # final name
        name = ("%s_%s_%s_%s_%s_%d_%d" % ( str_reg , str_lss , str_lts , str_sigma , str_date , model.nstep, model.geometry.shape[0] ) )

    return( name  )
    
