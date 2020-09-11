def plot_time_series( model ):
    """
    Plot time series
    """
    from progress.bar import Bar
    from time import time

    import multiprocessing
    from joblib import Parallel, delayed
    from tqdm import tqdm
    
    num_cores = multiprocessing.cpu_count()
    inputs = tqdm( model.Ogts.lcode() )
    
    def plot_ts_pyeq( site , model ):

        _dev_nul = model.Ogts.__dict__[ site ].plot( \
                                          superimposed=[ model.Mgts.__dict__[ site ] ], \
                                          date_unit='cal',\
                                          save_dir_plots = model.odir + '/plots/ts', \
                                          save = True, \
                                          show = False, \
                                          center = False, \
                                          verbose=  model.verbose,\
                                           )
        
    processed_list = Parallel(n_jobs=num_cores)(delayed( plot_ts_pyeq )(i, model ) for i in inputs)
