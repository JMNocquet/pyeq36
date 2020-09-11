def make_model_name( model ):
    """
    produces output name according to input parameters
    """
    # import
    
    import pyacs.lib.astrotime as at
    import numpy as np
    
    
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

    # start_date
    model_datetime = at.seconds2datetime( model.np_model_date_s[0] ) 
    str_yr = str(model_datetime.year)[-2:]
    doy = at.datetime2dayno( model_datetime)[1]
    duration = int( (model.np_model_date_s[-1] - model.np_model_date_s[-0] ) /86400. ) 

    str_date = ("%s%03d_%d" % ( str_yr , doy , int( duration ) ) )
 
    name = ("%s_%s_%s_%s_%d_%d" % ( str_sigma , str_dc , str_tau , str_date , model.nstep, model.geometry.shape[0] ) )
    
    return( name  )
    
