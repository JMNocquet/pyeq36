
def make_dir_pyeq_output( odir ):

    # import
    from os import mkdir, path

    # results directory
    if not path.exists(odir):mkdir(odir)
    
    # info directory
    if not path.exists(odir+'/info'):mkdir(odir+'/info')
    
    # GPS time series results directory
    if not path.exists(odir+'/time_series'):mkdir(odir+'/time_series')
    if not path.exists(odir+'/time_series/obs'):mkdir(odir+'/time_series/obs')
    if not path.exists(odir+'/time_series/model'):mkdir(odir+'/time_series/model')
    if not path.exists(odir+'/time_series/model_all_dates'):mkdir(odir+'/time_series/model_all_dates')
    if not path.exists(odir+'/time_series/res'):mkdir(odir+'/time_series/res')
    if not path.exists(odir+'/time_series/plot'):mkdir(odir+'/time_series/plot')

    # slip results directory
    if not path.exists(odir+'/slip'):mkdir(odir+'/slip')
    if not path.exists(odir+'/slip/incremental'):mkdir(odir+'/slip/incremental')
    if not path.exists(odir+'/slip/cumulative'):mkdir(odir+'/slip/cumulative')
    if not path.exists(odir+'/slip/rate'):mkdir(odir+'/slip/rate')
    
    # slip time series results directory
    if not path.exists(odir+'/slip_time_series'):mkdir(odir+'/slip_time_series')
    if not path.exists(odir+'/slip_time_series/rate'):mkdir(odir+'/slip_time_series/rate')
    if not path.exists(odir+'/slip_time_series/incremental'):mkdir(odir+'/slip_time_series/incremental')
    if not path.exists(odir+'/slip_time_series/cumulative'):mkdir(odir+'/slip_time_series/cumulative')

    # conf
    if not path.exists(odir+'/conf'):mkdir(odir+'/conf')

    # stf results directory
    if not path.exists(odir+'/stf'):mkdir(odir+'/stf')

    # sum results directory
    if not path.exists(odir+'/summary'):mkdir(odir+'/summary')

    # stats results directory
    if not path.exists(odir+'/stats'):mkdir(odir+'/stats')

    # plots results directory
    if not path.exists(odir+'/plots'):mkdir(odir+'/plots')

    # stf plots directory
    if not path.exists(odir+'/plots/stf'):mkdir(odir+'/plots/stf')

    # time series plots directory
    if not path.exists(odir+'/plots/ts'):mkdir(odir+'/plots/ts')

    # model cumulative plots directory
    if not path.exists(odir+'/plots/model_cumulative'):mkdir(odir+'/plots/model_cumulative')

    # model rate plots directory
    if not path.exists(odir+'/plots/model_rate'):mkdir(odir+'/plots/model_rate')

    # displacements directory
    if not path.exists(odir+'/displacement'):mkdir(odir+'/displacement')
    if not path.exists(odir+'/displacement/cumulative'):mkdir(odir+'/displacement/cumulative')
    if not path.exists(odir+'/displacement/cumulative/model'):mkdir(odir+'/displacement/cumulative/model')
    if not path.exists(odir+'/displacement/cumulative/obs'):mkdir(odir+'/displacement/cumulative/obs')
    if not path.exists(odir+'/displacement/cumulative/res'):mkdir(odir+'/displacement/cumulative/res')

    # npz directory
    if not path.exists(odir+'/npy'):mkdir(odir+'/npy')

    # shapefile
    if not path.exists(odir+'/shapefile'):mkdir(odir+'/shapefile')
    if not path.exists(odir+'/shapefile/slip_cumulative'):mkdir(odir+'/shapefile/slip_cumulative')
    if not path.exists(odir+'/shapefile/slip_rate'):mkdir(odir+'/shapefile/slip_rate')
    if not path.exists(odir+'/shapefile/disp_cumulative'):mkdir(odir+'/shapefile/disp_cumulative')
    if not path.exists(odir+'/shapefile/disp_incremental'):mkdir(odir+'/shapefile/disp_incremental')
    
    # gmt
    if not path.exists(odir+'/gmt'):mkdir(odir+'/gmt')
    if not path.exists(odir+'/gmt/slip_cumulative'):mkdir(odir+'/gmt/slip_cumulative')
    if not path.exists(odir+'/gmt/slip_rate'):mkdir(odir+'/gmt/slip_rate')

    