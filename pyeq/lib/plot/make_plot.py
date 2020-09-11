def make_plot( model , interpolation=False ):
    
    ###########################################################################
    # IMPORT
    ###########################################################################

    import pyeq.lib.plot
    import glob
    
    ###########################################################################
    # MAKE PLOT FOR TIME SERIES
    ###########################################################################
    
    #if not interpolation:
    #    print("-- Plotting time series")
    #    pyeq.lib.plot.plot_time_series( model )

    ###########################################################################
    # GENERATE SHAPEFILES OF CUMULATIVE SLIP FOR MODEL DISPLAY IN QGIS
    ###########################################################################
    
    # print('-- writes GMT and shapefiles for visualization in QGIS')
    # one_degree=111.1
    # TRIANGLE=True
    #
    # if interpolation:
    #
    #     lslip_dat=glob.glob(model.odir+"/slip/i_cumulative/*_cumulative_slip.dat")
    #     pyeq.lib.plot.model2shp_gmt(model.geometry, 'tde', lslip_dat, out_dir_shp=model.odir+'/shapefile/i_slip_cumulative', out_dir_gmt=model.odir+'/gmt/i_slip_cumulative' , verbose=model.verbose )
    #
    #     lslip_dat=glob.glob(model.odir+"/slip/i_rate/*_slip_rate.dat")
    #     pyeq.lib.plot.model2shp_gmt(model.geometry, 'tde', lslip_dat, out_dir_shp=model.odir+'/shapefile/i_slip_rate', out_dir_gmt=model.odir+'/gmt/i_slip_rate' , verbose=model.verbose )
    #
    # if not interpolation:
    #
    #     lslip_dat=glob.glob(model.odir+"/slip/cumulative/*_cumulative_slip.dat")
    #     pyeq.lib.plot.model2shp_gmt(model.geometry, 'tde', lslip_dat, out_dir_shp=model.odir+'/shapefile/slip_cumulative', out_dir_gmt=model.odir+'/gmt/slip_cumulative' , verbose=model.verbose )
    #
    #     lslip_dat=glob.glob(model.odir+"/slip/rate/*_slip_rate.dat")
    #     pyeq.lib.plot.model2shp_gmt(model.geometry, 'tde', lslip_dat, out_dir_shp=model.odir+'/shapefile/slip_rate', out_dir_gmt=model.odir+'/gmt/slip_rate' , verbose=model.verbose )

    ###########################################################################
    # MAKE PLOT FOR MODELS
    ###########################################################################
    
    print("-- plotting model")
    #if model.interseismic !=0:
    #    pyeq.lib.plot.plot_model_interseismic( model )
    #else:
    pyeq.lib.plot.plot_model( model , interpolation = interpolation )
    
    ###########################################################################
    # MAKE PLOT FOR STF
    ###########################################################################

    print("-- making plots for stf results in %s" % ( model.odir+'/plots/stf' ) )    
    pyeq.lib.plot.plot_stf( model )
    