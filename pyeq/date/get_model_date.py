def get_model_date( model ):

    import pyeq.date
    import pyacs.lib.astrotime as at
    import copy
    import numpy as np

    # old build
    if model.build == 1:
        print("-- Interpretating user requested model date using pyeq.lib.date_args.get_np_dates_from_arg: %s" % ( model.dates ) )
        
        np_dates = pyeq.date.get_np_dates_from_arg_old(model.dates, model.ts, model.verbose)
        date_ref=np_dates[0]
    
        if model.verbose:
            print('--- requested dates for model')
            
            for i in np.arange(np_dates.shape[0]):
        
                date=np_dates[i]
                (mday,month,ut) = at.decyear2cal( date )
                (noday,ut)      = at.decyear2dayno( date )    
                print(("%04d %10.5lf +%4.1lf days %02d-%02d-%04d-%.1lf %03d" % (i,date,(date-date_ref)*365.25,mday,month,int(date),ut,noday)))
            
        np_dates_save=copy.deepcopy(np_dates)

    # for the new pyacs >= 0.50.3
    if model.build in [ 2 , 3 , 4 ]:
        print("-- Interpretating user requested model date using pyeq.lib.date_args2.get_np_dates_from_arg2: %s" % ( model.dates ) )
        
        np_dates_save = pyeq.date.get_np_dates_from_arg(model.dates, model.np_obs_date_s, rounding=model.rounding, verbose=False)

    print("-- model dates requested by user: %d" % (np_dates_save.shape[0]) )
    for i in np.arange( np_dates_save.shape[0] ):
        print("%s" % ( at.seconds2datetime( np_dates_save[i] ).isoformat(' ')   ) )

    return np_dates_save