def plot_time_series( model ):
    """
    Plot time series
    """

    ###########################################################################
    # IMPORT
    ###########################################################################

    from progress.bar import Bar
    from time import time
    import multiprocessing
    from joblib import Parallel, delayed
    from tqdm import tqdm
    import numpy as np
    from os import mkdir, path

    import pyeq.message.message as MESSAGE
    import pyeq.message.verbose_message as VERBOSE
    import pyeq.message.error as ERROR
    import pyeq.message.warning as WARNING
    import pyeq.message.debug_message as DEBUG

    ###########################################################################
    # CHECK OUTPUT DIR
    ###########################################################################

    if not path.exists(model.odir + '/plots'): mkdir(model.odir + '/plots')

    ###########################################################################
    # MAP OF GPS SITES
    ###########################################################################

    if model.plot_settings.map:
        MESSAGE("Making map of GPS sites location: %s"% (model.odir + '/plots/map/map.png'))
        if not path.exists(  model.odir +'/plots/map' ):mkdir(  model.odir +'/plots/map' )
        model.Ogts.show_map( show=False, save = model.odir + '/plots/map/map.png' )

    num_cores = multiprocessing.cpu_count()
    inputs = tqdm( model.Ogts.lcode() )

    # plot
    def plot_ts_pyeq( site , model ):
        _dev_nul = model.Ogts.__dict__[ site ].plot(
            min_yaxis= float(model.plot_settings.ts_min_yaxis),
            superimposed=[ model.Mgts.__dict__[ site ] ],
            date_unit='cal',
            save_dir_plots = model.odir + '/plots/ts',
            save = True,
            show = False,
            center = False,
            verbose=  model.debug,
        )
        
    processed_list = Parallel(n_jobs=num_cores)(delayed( plot_ts_pyeq )(i, model ) for i in inputs)

