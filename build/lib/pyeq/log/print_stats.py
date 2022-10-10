def print_stats( model ):
    """
    Print inversion results statistics
    """
    
    ###########################################################################
    # IMPORT
    ###########################################################################

    import pyeq.obs_tensor.sgts2obs_tensor
    import pyeq.obs_tensor.obs_tensor2sgts
    import numpy as np
    import pyacs.lib.astrotime as at
    import pyacs.lib.utils
    
    ###########################################################################
    # INITIALIZE
    ###########################################################################
    
    T_MOD , np_names_t_obs, np_obs_date_s = pyeq.obs_tensor.sgts2obs_tensor.sgts2tensor(model.Mgts)
    T_RES , np_names_t_obs, np_obs_date_s = pyeq.obs_tensor.sgts2obs_tensor.sgts2tensor(model.Rgts)
    
    ###########################################################################
    # RMS PER SITE
    ###########################################################################
    
    format =("%10.2lf %10.2lf %10.2lf     %s")
    comment = '    East      North         Up  site (rms in mm)'
    
    rms   = np.sqrt( np.nansum (  T_RES[:,:,:3]**2 , axis=0 ) / np.nansum( T_RES[:,:,:3]*0.0 +1 , axis=0 ) )
    
    
    pyacs.lib.utils.save_np_array_with_string( rms , \
                                               model.np_gps_site ,\
                                              format, \
                                              model.odir+'/stats/rms.dat', \
                                              comment )

    ###########################################################################
    # STD PER SITE
    ###########################################################################
    
    format =("%10.2lf %10.2lf %10.2lf     %s")
    comment = '    East      North         Up  site (std in mm)'
    
    std = np.nanstd ( T_RES[:,:,:3] , axis=0 )
    
    
    pyacs.lib.utils.save_np_array_with_string( std , \
                                               model.np_gps_site ,\
                                              format, \
                                              model.odir+'/stats/std.dat', \
                                              comment )

    ###########################################################################
    # WRMS PER SITE
    ###########################################################################
    
    denom = np.nansum (  1./T_RES[:,:,3:]**2 , axis=0 )
    num   = np.nansum (  (T_RES[:,:,:3]/T_RES[:,:,3:])**2 , axis=0 )
    
    wrms = np.sqrt( num / denom )
    
    comment = '    East     North     Up  site (wrms in mm)'
    pyacs.lib.utils.save_np_array_with_string( wrms , \
                                               model.np_gps_site ,\
                                              format, \
                                              model.odir+'/stats/wrms.dat', \
                                              comment )
    
    ###########################################################################
    # BIAS PER SITE
    ###########################################################################

    bias = np.nanmean ( T_RES[:,:,:3] , axis=0 )

    
    comment = '    East     North     Up  site (bias in mm)'
    pyacs.lib.utils.save_np_array_with_string( bias , \
                                               model.np_gps_site ,\
                                              format, \
                                              model.odir+'/stats/bias.dat', \
                                              comment )
    
    ###########################################################################
    # CHI2
    ###########################################################################

    model.chi2_h   = np.nansum( ( T_RES[:,:,:2] / T_RES[:,:,3:5] )**2 ) 
    model.chi2_all = np.nansum( ( T_RES[:,:,:3] / T_RES[:,:,3:6] )**2 ) 
    
    model.reduced_chi2_h   = np.sqrt( np.nansum( ( T_RES[:,:,:2] / T_RES[:,:,3:5] )**2 ) / np.count_nonzero(~np.isnan( T_RES[:,:,:2] )) )
    model.reduced_chi2_all = np.sqrt( np.nansum( ( T_RES[:,:,:3] / T_RES[:,:,3:6] )**2 ) / np.count_nonzero(~np.isnan( T_RES[:,:,:2] )) )

    model.rms_h   = np.sqrt( np.nansum (  T_RES[:,:,:2]**2 ) / np.nansum( T_RES[:,:,:2]*0.0 +1  ) )
    model.rms_all = np.sqrt( np.nansum (  T_RES[:,:,:3]**2 ) / np.nansum( T_RES[:,:,:3]*0.0 +1  ) )

    # wrms    
    
    mask = np.any(np.isnan(T_RES[:,:,3:6]) | np.equal(T_RES[:,:,3:6], 0), axis=2)
    
    denom = np.nansum (  1./T_RES[:,:,3:6][~mask]**2 )
    num   = np.nansum (  (T_RES[:,:,:3][~mask]/T_RES[:,:,3:6][~mask])**2 )
    
    model.wrms_all = np.sqrt( num / denom )

    
    mask = np.any(np.isnan(T_RES[:,:,3:5]) | np.equal(T_RES[:,:,3:5], 0), axis=2)
    denom = np.nansum (  1./T_RES[:,:,3:5][~mask]**2 )
    num   = np.nansum (  (T_RES[:,:,:2][~mask]/T_RES[:,:,3:5][~mask])**2 )
    
    model.wrms_h = np.sqrt( num / denom )
    
    # bias (keep the component separated)
    model.bias = np.nanmean( np.nanmean ( T_RES[:,:,:3] , axis=0 ) , axis=0 )

    # WORST RMS SITE
    idx_worst = np.argmax( rms[:,0]**2 + rms[:,1]**2 )
    model.worst_site = model.np_gps_site[ idx_worst ]
    model.rms_worst = rms[ idx_worst, :]
    model.n_obs_inversion = np.count_nonzero(~np.isnan( T_RES[:,:,:2] ))

    ###########################################################################
    # PRINT STAT SUMMARY
    ###########################################################################
    fsumn = model.odir+'/stats/stats_sum.dat'
    fsum=open(fsumn,'w')

    fsum.write('#--------------------------------------------------------------------------\n')
    fsum.write('#   Inversion statistics: \n')
    fsum.write('#--------------------------------------------------------------------------\n')

    fsum.write("    rms  horizontal  (mm)            : %8.1f\n" % ( model.rms_h ))
    fsum.write("    rms  3D          (mm)            : %8.1f\n" % ( model.rms_all ))

    fsum.write("    wrms  horizontal (mm)            : %8.1f\n" % ( model.wrms_h ))
    fsum.write("    wrms  3D         (mm)            : %8.1f\n" % ( model.wrms_all ))


    fsum.write("    bias overall  east (mm)          : %8.1f\n" % model.bias[1])
    fsum.write("    bias overall north (mm)          : %8.1f\n" % model.bias[0])
    fsum.write("    bias overall    up (mm)          : %8.1f\n" % model.bias[2])

    fsum.write("    chi2 obs horizontal              : %15.1f\n" %  model.chi2_h )
    fsum.write("    chi2 obs 3D                      : %15.1f\n" %  model.chi2_all )

    fsum.write("    reduced chi2 obs horizontal      : %8.1f\n" %  model.reduced_chi2_h )
    fsum.write("    reduced chi2 obs 3D              : %8.1f\n" %  model.reduced_chi2_all )
    
    fsum.write("    number of observations           : %8d\n" %  model.n_obs_inversion )
    fsum.write("    worst site rms ENU (mm)          : %s %8.2lf %8.2lf %8.2lf\n" % ( model.worst_site , model.rms_worst[1] , model.rms_worst[0] , model.rms_worst[2] ) )

    fsum.close()

    ###########################################################################
    # RMS TIME SERIES
    ###########################################################################

    comment = ("Decimal days since model date start: %s" % ( at.seconds2datetime( model.np_model_date_s[0] ).isoformat(' ')))

    # HORIZONTAL
    rms_ts_h = np.sqrt( np.nansum( ( T_RES[:,:,0]**2 + T_RES[:,:,1]**2 ) , axis=1  )  / (2 * (~np.isnan( T_RES[:,:,0] )).sum(axis=1) ) )

    format = "%5.2lf  8%.2f  %s"

    to_print= np.vstack( ( ( model.np_obs_date_s -  model.np_obs_date_s[0]) / 86400. , rms_ts_h ) ).T

    model.ts_rms_2D = rms_ts_h

    pyacs.lib.utils.save_np_array_with_string( to_print , \
                                               model.np_obs_date_isoformat ,\
                                              format, \
                                              model.odir+'/stats/rms_ts_h.dat', \
                                              comment )

    # 3D
    rms_ts_all = np.sqrt( np.nansum( ( T_RES[:,:,0]**2 + T_RES[:,:,1]**2 + T_RES[:,:,2]**2 ) , axis=1  )  / (3 * (~np.isnan( T_RES[:,:,0] )).sum(axis=1) ) )

    format = "%5.2lf  8%.2f  %s"

    to_print= np.vstack( ( ( model.np_obs_date_s -  model.np_obs_date_s[0]) / 86400. , rms_ts_all ) ).T

    pyacs.lib.utils.save_np_array_with_string( to_print , \
                                               model.np_obs_date_isoformat ,\
                                              format, \
                                              model.odir+'/stats/rms_ts_3D.dat', \
                                              comment )

    
    
    return model 

