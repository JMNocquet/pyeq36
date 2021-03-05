def plot_model_interseismic( model ):
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
    import numpy as np
    import multiprocessing
    from joblib import Parallel, delayed
    from tqdm import tqdm


    ###################################################################
    # DEFINE THE UPPER BOUNDS VALUE FOR COLOR SCALE
    ###################################################################

    def get_max_scale( mx ):
        import numpy as np
        
        scale =  np.arange( 10**int(  np.log10( mx )  ) ,  10**(int(  np.log10( mx )  )+1)+ 10** int(  np.log10( mx )) , 10** int(  np.log10( mx )  ))
        
        return scale[np.where( scale > mx )][0] 
    
    ###################################################################
    # VALUES FOR PLOTS
    ###################################################################


    # custum color map
    cmap = plt.cm.magma_r  # define the colormap
    # extract all colors from the .jet map
    cmaplist = [cmap(i) for i in range(cmap.N)]
    # force the first color entry to be grey
    cmaplist[0] = (1., 1., 1., 1. )
    
    # create the new map
    cmap = mpl.colors.LinearSegmentedColormap.from_list(
        'my_cmap', cmaplist, cmap.N)


    lcumulative = glob( model.odir + '/shapefile/slip_cumulative/*_cumulative_slip.shp' )
    lrate = glob( model.odir + '/shapefile/slip_rate/*_slip_rate.shp' )

    exp = model.name
    
    num_cores = np.max( [ multiprocessing.cpu_count()-10 , 2 ] )
    
    ###################################################################
    # PLOT MODEL CUMULATIVE
    ###################################################################


    def plot_model_cumulative( shp , cmap , norm , exp , mxslip , outdir ):

        slip = geopandas.read_file( shp )
        
        fig, ax = plt.subplots()
        ax.set_aspect('equal')
        plot_polygon_collection(ax, slip['geometry'],alpha=230)
        plot_polygon_collection(ax, slip['geometry'][slip.slip>0], slip['slip'][slip.slip>0],cmap=cmap, norm=norm )

        sm = plt.cm.ScalarMappable(cmap=cmap , norm=norm )
        sm.set_array([0., get_max_scale( mxslip ) ])
        plt.colorbar(sm)
        title = ( "step #%s - %s" % ( shp.split('/')[-1][:4] , exp )  )
        plt.title( title )
        plot_name = ("%s/plots/model_cumulative/%s_cumulative.pdf" % ( outdir , shp.split('/')[-1][:4] ) )
        plt.savefig( plot_name )
        plt.close( fig )
    ###################################################################


        
    inputs = tqdm( lcumulative )
    mxslip = model.cumulative_slip_max
    bounds = np.linspace(0, get_max_scale( mxslip ), 21)
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

    processed_list = Parallel(n_jobs=num_cores)(delayed( plot_model_cumulative )( i, cmap , norm , exp , mxslip , model.odir ) for i in inputs )


    ###################################################################
    # PLOT MODEL RATE
    ###################################################################

    def plot_model_rate( shp , cmap , norm , exp , mxslip , outdir ):

        slip = geopandas.read_file( shp )
        fig, ax = plt.subplots()
        ax.set_aspect('equal')
        plot_polygon_collection(ax, slip['geometry'],alpha=230)
        
        plot_polygon_collection(ax, slip['geometry'], slip['slip'],cmap="magma" )

        sm = plt.cm.ScalarMappable(cmap="magma")
        sm.set_array([-47.5,mxslip])
        plt.colorbar(sm)
        title = ( "step #%s - %s" % ( shp.split('/')[-1][:4] , exp ) )
        plt.title( title )
        plot_name = ("%s/plots/model_rate/%s_rate.pdf" % ( model.odir , shp.split('/')[-1][:4] ) )
        plt.savefig( plot_name )
        plt.close( fig )
    ###################################################################

    
    mxslip = model.rate_slip_max
    inputs = tqdm( lrate )
    bounds = np.linspace(0, get_max_scale( mxslip ), 21)
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

    processed_list = Parallel(n_jobs=num_cores)(delayed( plot_model_rate )( i, cmap , norm , exp , mxslip , model.odir ) for i in inputs )
    
