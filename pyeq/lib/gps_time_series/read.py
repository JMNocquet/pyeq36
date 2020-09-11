def read( model ):

    # import
    import os
    import textwrap
    import numpy as np

    # new handling of dates in pyeq since October 2019 >= pyeq.0.50.3

    if os.path.isdir( model.dir_ts ):
        print("-- Reading GPS time series from directory: %s" % ( model.dir_ts ) )
        
        from pyacs.gts.Sgts import Sgts
        tts=Sgts( ts_dir=model.dir_ts, xyz=False , lexclude = model.lexclude_gps , verbose=model.verbose )
        print("--- Number of time series read in %s : %d " % ( model.dir_ts , len(tts.lcode())) )
        print("%s" % ( textwrap.fill(" ".join( tts.lcode() ) , width=80 ) ))

        # kept for backward compatibility
        if model.build == 1:
            if model.rounding=='day':
                print("--- Forcing time series to be daily")
                tts=tts.gts('force_daily')
        
        # keeps only sites present in GREEN
        
        np_code_in_ts_green = np.intersect1d(np.array( sorted( tts.lcode() ) , dtype=str ), model.name_obs, assume_unique=True, return_indices=True )

        if not( np.array_equal( np_code_in_ts_green[0], np.array(sorted( tts.lcode() ) ) )):
            print("!!!WARNING: some sites have time series but are not present in GREEN and therefore have no Green function values")
            print("!!! list of time series codes: %d" % ( len( tts.lcode() )))
            print("%s" % ( textwrap.fill(" ".join( tts.lcode() ) , width=80 ) ))
            print("!!! list of GREEN codes: %d" % ( model.name_obs.shape[0] ))
            print("%s" % ( textwrap.fill(" ".join( sorted( model.name_obs.tolist() ) ), width=80 ) ))
            print("!!! Keeping sites only present in GREEN and time series")
            tts =  tts.sub( linclude = np_code_in_ts_green[0].tolist() )
            print("--- Number of time series kept for inversion: %d " % ( len(tts.lcode())) )
            print("%s" % ( textwrap.fill(" ".join( tts.lcode() ) , width=80 ) ))

        ts = tts
        print("-- converting sgts time series to observation tensor using pyeq.lib.obs_tensor.sgts2tensor")

        import pyeq.lib.obs_tensor.sgts2obs_tensor
        T_OBS_RAW , np_names_t_obs, np_obs_date_s = pyeq.lib.obs_tensor.sgts2obs_tensor.sgts2tensor( ts, rounding=model.rounding , verbose=model.debug )

        model.t_obs_raw = T_OBS_RAW
        model.np_names_t_obs = np_names_t_obs 
        model.np_obs_date_s = np_obs_date_s

        if model.verbose:
            print("-- raw_observation_tensor shape: " , model.t_obs_raw.shape)


    model.n_site_input_ts = np_names_t_obs.shape[0]
    model.n_date_input_ts = np_obs_date_s.shape[0]

    return model
