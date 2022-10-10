def plot_stf( model ):
    """
    Make plot for stf
    """
    
    # import
    import matplotlib.pyplot as plt
    import numpy as np
    from pandas.plotting import register_matplotlib_converters
    import os

    register_matplotlib_converters()

    # create output directpory
    os.makedirs( os.path.dirname(model.odir+'/plots/stf/'), exist_ok=True )

    # MAKE PLOTS
    # STF
    fig, (ax1 , ax2) = plt.subplots(2,1, figsize=(6,8))

    ax1.plot( model.np_mid_model_delta_d , model.STF , '-bo' , markersize=2)
    ax1.set_title( ("STF - %s " % model.name ) )
    ax1.set_xlabel('days')
    ax1.set_ylabel('Moment rate (N.m)')
    # avoid infinity in the logarithm
    mag = np.where( model.STF == 0, np.nan, model.STF ) 
    mag = 2./3.*(np.log10( mag ) -9.1 )  
    ax2.set_xlim( model.np_mid_model_datetime[0] , model.np_mid_model_datetime[-1] )
    ax2.plot( model.np_mid_model_datetime , mag , '-bo' , markersize=2)

    for label in ax2.get_xticklabels():
        label.set_rotation(70)
        label.set_horizontalalignment('right')
    ax2.set_ylabel('Equivalent Magnitude')
    plt.savefig( model.odir+'/plots/stf/stf.pdf')
    plt.close( fig )

    # FAULT_AVERAGE_SLIP_RATE
    fig, (ax1 , ax2) = plt.subplots(2,1, figsize=(6,8))
    ax1.plot( model.np_mid_model_delta_d , model.FAULT_AVERAGE_SLIP_RATE , '-bo' , markersize=2)
    ax1.set_title( ("FAULT_AVERAGE_SLIP_RATE - %s " % model.name ))
    ax1.set_xlabel('days')
    ax1.set_ylabel('mm')
    ax2.plot( model.np_mid_model_datetime , model.FAULT_AVERAGE_SLIP_RATE , '-bo' , markersize=2 )

    for label in ax2.get_xticklabels():
        label.set_rotation(70)
        label.set_horizontalalignment('right')
    ax2.set_ylabel('mm')
    plt.savefig( model.odir+'/plots/stf/fault_average_slip_rate.pdf')
    plt.close( fig )

    # CSTF
    fig, (ax1 , ax2) = plt.subplots(2,1, figsize=(6,8))
    ax1.plot( model.np_model_delta_d , model.CSTF , '-bo' , markersize=2)
    ax1.set_title( ("CSTF - %s " % model.name ) )
    ax1.set_xlabel('days')
    ax1.set_ylabel('Cumulated Moment (N.m)')
    # avoid infinity in the logarithm
    mag = np.where( model.CSTF == 0, np.nan, model.CSTF ) 
    mag = 2./3.*(np.log10( mag ) -9.1 )  
    ax2.set_xlim( model.np_model_datetime[0] , model.np_model_datetime[-1] )
    ax2.plot( model.np_model_datetime , mag , '-bo' , markersize=2)
    ax2.set_xlim( model.np_model_datetime[0] , model.np_model_datetime[-1] )
    for label in ax2.get_xticklabels():
        label.set_rotation(70)
        label.set_horizontalalignment('right')
    ax2.set_ylabel('Equivalent Magnitude')
    plt.savefig( model.odir+'/plots/stf/cstf.pdf')
    plt.close( fig )

    # RMS
    fig, (ax1 , ax2) = plt.subplots(2,1, figsize=(6,8))
    ax1.plot( ( model.np_obs_date_s -  model.np_obs_date_s[0]) / 86400. , model.ts_rms_2D , '-bo' , markersize=2)
    ax1.set_title( ("RMS 2D - %s " % model.name ) )
    ax1.set_xlabel('days')
    ax1.set_ylabel('mm')
    ax2.plot( model.np_obs_datetime , model.ts_rms_2D , '-bo' , markersize=2)
    for label in ax2.get_xticklabels():
        label.set_rotation(70)
        label.set_horizontalalignment('right')
    plt.savefig( model.odir+'/plots/stf/rms_2D.pdf')
    plt.close( fig )
    
    # MAX SLIP (CUMULATIVE AND RATE)
    
    ts_max_slip_rate = model.RATE_SLIP_PER_TIME_STEP[:,model.idx_cumulative_slip_max,:]
    ts_max_slip_cum  = model.CUMULATIVE_SLIP_PER_TIME_STEP[:,model.idx_cumulative_slip_max,:]

    fig, (ax1 , ax2) = plt.subplots(2,1, figsize=(6,8.5))
    ax1.plot( model.np_model_delta_d , ts_max_slip_cum , '-bo' , markersize=2)
    title = ( ("max cumulative slip - subfault #%d - %s " % ( model.idx_cumulative_slip_max , model.name ) ) )
    ax1.set_title(title)
    ax1.set_xlabel('days')
    ax1.set_ylabel('mm')
    
    ax2.plot( model.np_mid_model_datetime , ts_max_slip_rate , '-bo' , markersize=2 )
    title = ("max slip rate - subfault #%d" % (model.idx_cumulative_slip_max) )
    ax2.set_title(title)
    ax2.set_ylabel('mm/day')

    for label in ax2.get_xticklabels():
        label.set_rotation(70)
        label.set_horizontalalignment('right')
    plt.savefig( model.odir+'/plots/stf/max_slip.pdf')
    plt.close( fig )
    
    