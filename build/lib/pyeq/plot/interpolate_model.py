def interpolate_model( model , geometry ):
    
    
    ###########################################################################
    # IMPORT
    ###########################################################################
    from os import mkdir, path
    import numpy as np
    import scipy.interpolate
    import pyacs.lib.astrotime as at
    import pyeq.plot

    ###########################################################################
    # DIRECTORY
    ###########################################################################

    if not path.exists( model.odir+'/slip/i_incremental'):mkdir( model.odir+'/slip/i_incremental')
    if not path.exists( model.odir+'/slip/i_cumulative'):mkdir( model.odir+'/slip/i_cumulative')
    if not path.exists( model.odir+'/slip/i_rate'):mkdir( model.odir+'/slip/i_rate')

    if not path.exists( model.odir+'/shapefile/i_slip_incremental'):mkdir( model.odir+'/shapefile/i_slip_incremental')
    if not path.exists( model.odir+'/shapefile/i_slip_cumulative'):mkdir( model.odir+'/shapefile/i_slip_cumulative')
    if not path.exists( model.odir+'/shapefile/i_slip_rate'):mkdir( model.odir+'/shapefile/i_slip_rate')

    if not path.exists( model.odir+'/gmt/i_slip_incremental'):mkdir( model.odir+'/gmt/i_slip_incremental')
    if not path.exists( model.odir+'/gmt/i_slip_cumulative'):mkdir( model.odir+'/gmt/i_slip_cumulative')
    if not path.exists( model.odir+'/gmt/i_slip_rate'):mkdir( model.odir+'/gmt/i_slip_rate')


    ###########################################################################
    # CUMULATIVE SLIP
    ###########################################################################

    CUMULATIVE_SLIP_TIME_STEP=np.zeros(( model.nfaults, 3 ))
    CUMULATIVE_SLIP_TIME_STEP[:,0] = model.sgeometry.centroid_long
    CUMULATIVE_SLIP_TIME_STEP[:,1] = model.sgeometry.centroid_lat


    for i in np.arange( model.np_model_date_s.shape[0] ):
        isodate = model.np_model_date_isoformat[i]

        if model.verbose:
            print("--- writing cumulative slip for model time : %04d %s" % (i ,  isodate ) )

        date_info=("# model date  %04d %s" % (i ,  isodate ) )

        # CUMULATIVE SLIP
        CUMULATIVE_SLIP_TIME_STEP[:,2] = model.CUMULATIVE_SLIP_PER_TIME_STEP[i].reshape( model.nfaults, -1 ).flatten()
        fname=model.odir+"/slip/i_cumulative/"+("%04d_cumulative_slip.dat" % (i))

        # INTERPOLATION
        ip = scipy.interpolate.LinearNDInterpolator( model.geometry[:,9:11], CUMULATIVE_SLIP_TIME_STEP[:,2], fill_value=0., rescale=False)
        ISLIP = ip( geometry[:,9:11] )

        # SAVE RESULTS
        ICUMULATIVE_SLIP_TIME_STEP = np.hstack( ( geometry[:,9:11] , ISLIP.reshape(-1,1) )   )
        np.savetxt(fname, ICUMULATIVE_SLIP_TIME_STEP, fmt="%10.5lf %10.5lf %10.3lf", header=date_info)
    
    
    
    ###########################################################################
    # LOOP ON MODEL TIME STEP FOR SLIP RATE & INCREMENTAL SLIP
    ###########################################################################

    RATE_SLIP_TIME_STEP = np.zeros(( model.nfaults, 3 ))
    RATE_SLIP_TIME_STEP[:,0] = model.sgeometry.centroid_long
    RATE_SLIP_TIME_STEP[:,1] = model.sgeometry.centroid_lat

    INC_SLIP_TIME_STEP = np.zeros(( model.nfaults, 3 ))
    INC_SLIP_TIME_STEP[:,0] = model.sgeometry.centroid_long
    INC_SLIP_TIME_STEP[:,1] = model.sgeometry.centroid_lat
    
    
    for i in np.arange( model.np_model_step_duration_s.shape[0] ):
        if model.verbose:
            print("--- calculating incremental slip, slip rate for model time step: %d" % (i) )
        sdatetime = model.np_model_datetime[i]
        edatetime = model.np_model_datetime[i+1]
        iso_sdate = model.np_model_date_isoformat[i]
        iso_edate = model.np_model_date_isoformat[i+1]
        delta_d   = model.np_model_step_duration_days[i]
        delta_t0  = model.np_model_delta_d[i]
        syear, sdoy, _ut  = at.datetime2dayno( sdatetime )
        eyear, edoy, _ut  = at.datetime2dayno( edatetime )

        sdecyear = at.datetime2decyear( sdatetime )
        edecyear = at.datetime2decyear( edatetime )

        date_info=("step #%04d  %s -> %s %8.3lf %8.3lf %04d %03d -> %04d %03d %15.10lf -> %15.10lf \n" % \
                   ( i, iso_sdate, iso_edate, delta_d, delta_t0, syear, sdoy, eyear, edoy, sdecyear, edecyear ))
        
        # SLIP RATE
        RATE_SLIP_TIME_STEP[ : , 2 ]  =  model.RATE_SLIP_PER_TIME_STEP[i].reshape( model.nfaults, -1 ).flatten()
        fname=model.odir+"/slip/i_rate/"+("%04d_slip_rate.dat" % ( i ))

        # INTERPOLATION
        ip = scipy.interpolate.LinearNDInterpolator( model.geometry[:,9:11], RATE_SLIP_TIME_STEP[:,2], fill_value=0., rescale=False)
        ISLIP = ip( geometry[:,9:11] )

        # SAVE RESULTS
        IRATE_SLIP_TIME_STEP = np.hstack( ( geometry[:,9:11] , ISLIP.reshape(-1,1) )   )
        np.savetxt(fname, IRATE_SLIP_TIME_STEP, fmt="%10.5lf %10.5lf %10.3lf", header=date_info)
            
        # INCREMENTAL SLIP
        INC_SLIP_TIME_STEP[:,2] = model.INC_SLIP_PER_TIME_STEP[i].reshape( model.nfaults, -1 ).flatten()
        fname=model.odir+"/slip/i_incremental/"+("%04d_delta_slip.dat" % ( i ))

        # INTERPOLATION
        ip = scipy.interpolate.LinearNDInterpolator( model.geometry[:,9:11], INC_SLIP_TIME_STEP[:,2], fill_value=0., rescale=False)
        ISLIP = ip( geometry[:,9:11] )

        # SAVE RESULTS
        IINC_SLIP_TIME_STEP = np.hstack( ( geometry[:,9:11] , ISLIP.reshape(-1,1) )   )
        np.savetxt(fname, IINC_SLIP_TIME_STEP, fmt="%10.5lf %10.5lf %10.3lf", header=date_info)
    
    ###########################################################################
    # SHAPEFILE
    ###########################################################################
    
    print('-- writes GMT a shapefiles for visualization in QGIS')
    one_degree=111.1
    TRIANGLE=True    
    
    import glob
    
    lslip_dat=glob.glob(model.odir+"/slip/i_cumulative/*_cumulative_slip.dat")
    pyeq.plot.model2shp_gmt(geometry, 'tde', lslip_dat, out_dir_shp=model.odir + '/shapefile/i_slip_cumulative', out_dir_gmt=model.odir + '/gmt/i_slip_cumulative', verbose=model.verbose)
    
    lslip_dat=glob.glob(model.odir+"/slip/i_rate/*_slip_rate.dat")
    pyeq.plot.model2shp_gmt(geometry, 'tde', lslip_dat, out_dir_shp=model.odir + '/shapefile/i_slip_rate', out_dir_gmt=model.odir + '/gmt/i_slip_rate', verbose=model.verbose)
    
    ###########################################################################
    # PLOT
    ###########################################################################

    pyeq.plot.plot_model(model, interpolation=True)
    
    return
    