def print_slip_time_series( model ):
    """
    Print rate, incremental and cumulative slip time series as text files
    """
    
    # import
    import numpy as np
    import pyacs.lib.utils

    # START LOOP ON FAULT ELEMENTS
    for isubfault in np.arange(  model.nfaults ):
        # SLIP RATE
        model.DELTA_D_AND_STS_RATE = np.vstack( ( model.np_mid_model_delta_d , model.RATE_SLIP_PER_TIME_STEP[:,isubfault,:].T ) ).T
        fname_rate = ("%s/slip_time_series/rate/%04d_rate.dat" % (model.odir,isubfault))
        # 1 or 2 rakes
        if model.RATE_SLIP_PER_TIME_STEP.shape[2] == 1:
            format = "%5.2lf  %10.3lf  %s"
        else:
            format = "%5.2lf  %10.3lf  %10.3lf  %s"
            
        comment_rate = ("subfault: #%04d . units are mm/day. dates are at the middle of model time steps. Decimal days since model date start: %s" % \
                        ( isubfault, model.np_model_date_isoformat[0] ) )

        pyacs.lib.utils.save_np_array_with_string( model.DELTA_D_AND_STS_RATE , \
                                                  model.np_mid_model_isoformat ,\
                                                  format, \
                                                  fname_rate, \
                                                  comment_rate )


        # INCREMENTAL SLIP
        model.DELTA_D_AND_STS_INC = np.vstack( ( model.np_mid_model_delta_d , model.INC_SLIP_PER_TIME_STEP[:,isubfault,:].T ) ).T
        fname_inc  = ("%s/slip_time_series/incremental/%04d_incremental.dat" % (model.odir,isubfault))
        comment_inc  = ("subfault: #%04d . units are mm. dates are at the middle of model time steps. Decimal days since model date start: %s" % \
                        ( isubfault, model.np_model_date_isoformat[0] ) )
        pyacs.lib.utils.save_np_array_with_string( model.DELTA_D_AND_STS_INC , \
                                                  model.np_mid_model_isoformat ,\
                                                  format, \
                                                  fname_inc, \
                                                  comment_inc )

        # CUMULATIVE SLIP
        model.DELTA_D_AND_STS_CUM = np.vstack( ( model.np_model_delta_d , model.CUMULATIVE_SLIP_PER_TIME_STEP[:,isubfault,:].T ) ).T
        fname_cum  = ("%s/slip_time_series/cumulative/%04d_cumulative.dat" % (model.odir,isubfault))
        comment_cum  = ("subfault: #%04d . units are mm. dates are model dates. Decimal days since model date start: %s" % \
                        ( isubfault, model.np_model_date_isoformat[0] ) )
        pyacs.lib.utils.save_np_array_with_string( model.DELTA_D_AND_STS_CUM , \
                                                  model.np_model_date_isoformat ,\
                                                  format, \
                                                  fname_cum, \
                                                  comment_cum )


    return model

