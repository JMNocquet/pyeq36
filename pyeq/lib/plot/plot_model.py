def plot_model( model , interpolation=False ):
    """
    plot models
    """
    
    ###################################################################
    # IMPORT
    ###################################################################

    import geopandas
    from glob import glob
    import matplotlib as mpl
    import matplotlib.pylab as plt
    from geopandas.plotting import plot_polygon_collection    
    import matplotlib.colors as colors
    from matplotlib.tri import Triangulation, TriAnalyzer, UniformTriRefiner    
    import numpy as np
    import multiprocessing
    from joblib import Parallel, delayed
    from tqdm import tqdm
    from os import mkdir, path
    from pyeq.lib.plot import plot_model_shp

    ###################################################################
    # NUMBER OF CPU FOR PARALLEL JOBS
    ###################################################################
    num_cores = np.max( [ multiprocessing.cpu_count()-10 , 2 ] )

    ###################################################################
    # DEFINE THE UPPER BOUNDS VALUE FOR COLOR SCALE
    ###################################################################

    def get_max_scale( mx ):
        import numpy as np
        scale =  np.arange( 10**int(  np.log10( mx )  ) ,  10**(int(  np.log10( mx )  )+1)+ 10** int(  np.log10( mx )) , 10** int(  np.log10( mx )  ))
        return scale[np.where( scale > mx )][0]

    ###################################################################
    # CASE INTERPOLATION
    ###################################################################

    if interpolation:
        if not path.exists(  model.odir +'/plots/i_model_cumulative' ):mkdir(  model.odir +'/plots/i_model_cumulative' )
        if not path.exists(  model.odir +'/plots/i_model_rate' ):mkdir(  model.odir +'/plots/i_model_rate' )
    else:
        if not path.exists(  model.odir +'/plots/model_cumulative' ):mkdir(  model.odir +'/plots/model_cumulative' )
        if not path.exists(  model.odir +'/plots/model_rate' ):mkdir(  model.odir +'/plots/model_rate' )


    ###################################################################
    # VARIABLES FOR LOOP OF CUMULATIVE MODELS
    ###################################################################
        
    mxslip = model.cumulative_slip_max * 1.1
    bounds = np.linspace(0, get_max_scale( mxslip ), 21)

    # title
    # load cstf
    H_title = {}
    H_shp   = {}
    H_ishp   = {}
    H_disp   = {}

    cstf = np.genfromtxt( model.odir+'/stf/cstf.dat' , dtype=str)

    for i in np.arange( cstf.shape[0] ):
        H_title[i] = ("%s -- time step %04d \n\n day %04d - %s %s - Mo=%s N.m" % ( model.name , i, int(float((cstf[i,0]))), cstf[i,2] , cstf[i,3], cstf[i,1] ))
        H_shp[i] = ("%s/shapefile/slip_cumulative/%04d_cumulative_slip.shp" % ( model.odir , i ) )
        H_ishp[i] = ("%s/shapefile/i_slip_cumulative/%04d_cumulative_slip.shp" % ( model.odir , i ) )
        H_disp[i] = ("%s/shapefile/disp_cumulative/%04d_model_cum_disp.shp" % ( model.odir , i ) )

    if interpolation:
        outdir = model.odir+'/plots/i_model_cumulative'
        print('-- making plot for interpolated cumulative slip in ' , outdir )

    if not interpolation:
        outdir = model.odir+'/plots/model_cumulative'
        if not path.exists(  outdir ):mkdir(  outdir )
        print('-- making plot for cumulative slip in ' , outdir )

    # cumulative slip plot
    processed_list = Parallel(n_jobs=num_cores)\
        (delayed( plot_model_shp ) \
                                    (H_shp[i], \
                                     cmap='magma_r', \
                                     ncolor=20, \
                                     contour=None, \
                                     crs=None, \
                                     log=False, \
                                     title=H_title[i], \
                                     bounds=bounds, \
                                     outdir=outdir, \
                                     outfile=None, \
                                     shp_poly=model.external_shapefile_poly, \
                                     shp_line=model.external_shapefile_line, \
                                     shp_point=H_disp[i])
        for i in tqdm(np.arange(cstf.shape[0])))

    outdir = model.odir+'/plots/model_cumulative_contour'
    if not path.exists( outdir ):mkdir( outdir )
    print('-- making contour plot for cumulative slip in ' , outdir )

    # cumulative slip contour plot
    processed_list = Parallel(n_jobs=num_cores)\
        (delayed( plot_model_shp ) \
                                    (H_shp[i], \
                                     cmap='magma_r', \
                                     ncolor=20, \
                                     contour=2, \
                                     crs=None, \
                                     log=False, \
                                     title=H_title[i], \
                                     bounds=[0.,mxslip], \
                                     outdir=outdir, \
                                     outfile=None, \
                                     shp_poly=model.external_shapefile_poly, \
                                     shp_line=model.external_shapefile_line, \
                                     shp_point=H_disp[i])
        for i in tqdm(np.arange(cstf.shape[0])))

    ###################################################################
    # PLOT MODEL RATE
    ###################################################################

    ###################################################################
    # VARIABLES FOR LOOP OF RATE MODELS
    ###################################################################
        
    mxslip = model.rate_slip_max * 1.1
    bounds = np.linspace(0, get_max_scale( mxslip ), 21)
    H_title = {}
    H_shp   = {}
    H_ishp   = {}
    stf = np.genfromtxt( model.odir+'/stf/stf.dat' , dtype=str)

    for i in np.arange( stf.shape[0] ):
        H_title[i] = ("%s -- time step %04d \n\n day %04d - %s %s - Mo=%s N.m" % ( model.name , i, int(float((stf[i,0]))), stf[i,2] , stf[i,3], stf[i,1] ))
        H_shp[i] = ("%s/shapefile/slip_rate/%04d_slip_rate.shp" % ( model.odir , i ) )
        H_ishp[i] = ("%s/shapefile/i_slip_rate/%04d_slip_rate.shp" % ( model.odir , i ) )
    # interpolation/no interpolation case
    if interpolation:
        outdir = model.odir+'/plots/i_model_rate'
        print('-- making plot for interpolated rate slip in ' , outdir )

    if not interpolation:
        outdir = model.odir+'/plots/model_rate'
        print('-- making plot for rate slip in ' , outdir )

    # model rate plot (linear scale)
    processed_list = Parallel(n_jobs=num_cores)\
        (delayed( plot_model_shp ) \
                                    (H_shp[i], \
                                     cmap='magma_r', \
                                     ncolor=None, \
                                     contour=None, \
                                     crs=None, \
                                     log=False, \
                                     title=H_title[i], \
                                     bounds=bounds, \
                                     outdir=outdir, \
                                     outfile=None, \
                                     shp_poly=model.external_shapefile_poly, \
                                     shp_line=model.external_shapefile_line, \
                                     shp_point=H_disp[i+1])
         for i in tqdm(np.arange(stf.shape[0])))

    # model rate contour plot (log scale)
    bounds = np.linspace(np.log10(0.1), np.log10(get_max_scale( mxslip )), 1001)
    outdir = model.odir+'/plots/model_rate_contour'
    if not path.exists(  outdir ):mkdir(  outdir )
    print('-- making contour plot for rate slip in ' , outdir )

    processed_list = Parallel(n_jobs=num_cores)\
        (delayed( plot_model_shp ) \
                                    (H_shp[i], \
                                     cmap='jet', \
                                     ncolor=None, \
                                     contour=1, \
                                     crs=None, \
                                     log=True, \
                                     title=H_title[i], \
                                     bounds=bounds, \
                                     outdir=outdir, \
                                     outfile=None, \
                                     shp_poly=model.external_shapefile_poly, \
                                     shp_line=model.external_shapefile_line, \
                                     shp_point=H_disp[i+1])
         for i in tqdm(np.arange(stf.shape[0])))

