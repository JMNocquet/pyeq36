def read( model ):

    ###################################################################
    # import
    ###################################################################

    import os
    import textwrap
    import numpy as np

    import pyeq.message.warning as WARNING
    import pyeq.message.message as MESSAGE
    import pyeq.message.verbose_message as VERBOSE
    import pyeq.message.debug_message as DEBUG
    import pyeq.message.error as ERROR

    # ensure model.lexclude_gps is unique
    model.lexclude_gps = list(set(sorted(model.lexclude_gps)))

    model.warning +="#GPS sites excluded from read in gps_time_series.read\n"
    model.warning +=("%s" % ( textwrap.fill(" ".join( model.lexclude_gps ) , width=80 ) ))
    model.warning += "\n"

    # new handling of dates in pyeq since October 2019 >= pyeq.0.50.3

    if os.path.isdir( model.dir_ts ):

        from pyacs.gts.Sgts import Sgts
        VERBOSE("user request excluded sites: %d" % len(model.lexclude_gps))
        DEBUG("%s" % (textwrap.fill(" ".join(model.lexclude_gps), width=80)))
        VERBOSE("Reading GPS time series from: %s" % ( model.dir_ts ) )
        try:
            tts=Sgts( ts_dir=model.dir_ts, xyz=False , lexclude = model.lexclude_gps , verbose=False )
        except:
            ERROR(("Could not read time series from directory: %s" % (model.dir_ts)),exit=True)

        VERBOSE("Number of time series read: %d " % ( len(tts.lcode())) )
        DEBUG("%s" % ( textwrap.fill(" ".join( tts.lcode() ) , width=80 ) ))

        # kept for backward compatibility
        if model.build == 1:
            if model.rounding=='day':
                MESSAGE("Old PYEQ feature for build 1. Forcing time series to be daily")
                tts=tts.gts('force_daily')
        
        # keeps only sites present in GREEN

        VERBOSE("Checking that GREEN and time series have the same sites")
        np_code_in_ts_green = np.intersect1d(np.array( sorted( tts.lcode() ) , dtype=str ),\
                                             model.name_obs, assume_unique=True, return_indices=True )

        if not( np.array_equal( np_code_in_ts_green[0], np.array(sorted( tts.lcode() ) ) )):
            WARNING("Check your green functions.")
            WARNING("some sites have time series but are not present in GREEN.")
            WARNING("list of time series codes: %d" % ( len( tts.lcode() )))
            WARNING("%s" % ( textwrap.fill(" ".join( tts.lcode() ) , width=80 ) ))
            WARNING("list of GREEN codes: %d" % ( model.name_obs.shape[0] ))
            WARNING("%s" % ( textwrap.fill(" ".join( sorted( model.name_obs.tolist() ) ), width=80 ) ))
            WARNING("Keeping sites only present in GREEN and time series")
            tts =  tts.sub( linclude = np_code_in_ts_green[0].tolist() )
            VERBOSE("Number of time series kept for inversion: %d " % ( len(tts.lcode())) )
            DEBUG("%s" % ( textwrap.fill(" ".join( tts.lcode() ) , width=80 ) ))

        ts = tts
        VERBOSE("converting sgts time series to observation tensor")

        import pyeq.obs_tensor.sgts2obs_tensor
        T_OBS_RAW , np_names_t_obs, np_obs_date_s = pyeq.obs_tensor.sgts2obs_tensor.sgts2tensor(ts, rounding=model.rounding, verbose=model.debug)

        model.t_obs_raw = T_OBS_RAW
        model.np_names_t_obs = np_names_t_obs 
        model.np_obs_date_s = np_obs_date_s

        VERBOSE(("observation_tensor shape (n_obs_dates, nsites, ncomponents+nsigma): %s" % (model.t_obs_raw.shape,)))


    else:
        ERROR(("%s is not a directory. Could not read time series" % (model.dir_ts)),exit=True)


    model.n_site_input_ts = np_names_t_obs.shape[0]
    model.n_date_input_ts = np_obs_date_s.shape[0]

    return model
