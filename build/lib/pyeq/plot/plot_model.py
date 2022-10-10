def plot_model(model, interpolation=False):
    """
    wrapper to plot_model_shp for plotting various models
    """

    ###################################################################
    # IMPORT
    ###################################################################

    import numpy as np
    import multiprocessing
    from joblib import Parallel, delayed
    from tqdm import tqdm
    from os import mkdir, path
    from str2bool import str2bool
    from pyeq.plot import plot_model_shp
    import pyacs.lib.units
    import pyeq.message.message as MESSAGE
    import pyeq.message.verbose_message as VERBOSE
    import pyeq.message.error as ERROR
    import pyeq.message.warning as WARNING
    import pyeq.message.debug_message as DEBUG

    #TODO: not sure this works when model time step are different from observation dates

    ###################################################################
    # NUMBER OF CPU FOR PARALLEL JOBS
    ###################################################################
    num_cores = np.max([multiprocessing.cpu_count() - 10, int(model.plot_settings.ncpu)])

    ###################################################################
    # DEFINE THE UPPER BOUNDS VALUE FOR COLOR SCALE
    ###################################################################

    def get_max_scale(mx):
        import numpy as np
        scale = np.arange(10 ** int(np.log10(mx)), 10 ** (int(np.log10(mx)) + 1) + 10 ** int(np.log10(mx)),
                          10 ** int(np.log10(mx)))
        return scale[np.where(scale > mx)][0]


    mxslip = model.cumulative_slip_max * 1.1

    # color bounds
    if model.plot_settings.slip_bounds in ['None','NONE']:
        bound_color = [0., mxslip]
    else:
        # user provided bounds [min_slip, max_slip]
        lbounds = model.plot_settings.slip_bounds.replace('[','').replace(']','').split(',')
        # deal with upper bound
        if lbounds[1] in ['None','NONE']:
            ub = mxslip
        else:
            ub = float( lbounds[1] )
        # deal with lower bound
        # lower bound provided as percent
        if '%' in lbounds[0]:
            lb = float(lbounds[0].split('%')[0])/100. * mxslip
        else:
        # lower bound as true value
            lb = float( lbounds[0] )

        bound_color = [lb,ub]

    # ncolor
    if model.plot_settings.ncolor in ['None','NONE'] or model.plot_settings.ncolor is None:
        ncolor = None
    else:
        ncolor = int(model.plot_settings.ncolor)

    # size displacement arrows
    sda = 1./float(model.plot_settings.disp_scale)


    ###################################################################
    # CUMULATIVE MODELS
    ###################################################################


    # title
    # load cstf
    H_title = {}
    H_shp = {}
    H_ishp = {}
    H_disp = {}
    H_obs = {}
    H_pred = {}
    H_res = {}

    H_obs_up = {}
    H_pred_up = {}
    H_res_up = {}

    cstf = np.genfromtxt(model.odir + '/stf/cstf.dat', dtype=str)

    # INITIALIZE DATA TO BE PLOTTED
    for i in np.arange(cstf.shape[0]):
        moment = float(cstf[i, 1])
        if moment > 0:
            magnitude = pyacs.lib.units.moment_to_magnitude(moment)
        else:
            magnitude = 0.

        H_title[i] = ("%s time step %04d - day %04d - %s %s\nMo=%s N.m (Mw %4.2lf)" % (
        model.name, i, int(float((cstf[i, 0]))), cstf[i, 2], cstf[i, 3], cstf[i, 1], magnitude))
        H_shp[i] = ("%s/shapefile/slip_cumulative/%04d_cumulative_slip.shp" % (model.odir, i))
        H_ishp[i] = ("%s/shapefile/i_slip_cumulative/%04d_cumulative_slip.shp" % (model.odir, i))
        H_disp[i] = ("%s/shapefile/disp_cumulative/%04d_model_cum_disp.shp" % (model.odir, i))
        H_obs[i] = ("%s/displacement/cumulative/obs/%04d_obs_cum_disp.dat" % (model.odir, i))
        H_pred[i] = ("%s/displacement/cumulative/model/%04d_model_cum_disp.dat" % (model.odir, i))
        H_res[i] = ("%s/displacement/cumulative/res/%04d_res_cum_disp.dat" % (model.odir, i))

        H_obs_up[i] = ("%s/displacement/cumulative/obs/%04d_obs_cum_disp_up.dat" % (model.odir, i))
        H_pred_up[i] = ("%s/displacement/cumulative/model/%04d_model_cum_disp_up.dat" % (model.odir, i))
        H_res_up[i] = ("%s/displacement/cumulative/res/%04d_res_cum_disp_up.dat" % (model.odir, i))


    ###########################################################################
    # cumulative slip plot
    ###########################################################################
    if str2bool(model.plot_settings.cum_slip):
        MESSAGE("plotting cumulative slip")

        if interpolation:
            if not path.exists(model.odir + '/plots/i_model_cumulative'): mkdir(
                model.odir + '/plots/i_model_cumulative')
            outdir = model.odir + '/plots/i_model_cumulative'
            MESSAGE("making plot for interpolated cumulative slip in %s" % outdir)
        else:
            if not path.exists(model.odir + '/plots/model_cumulative'): mkdir(model.odir + '/plots/model_cumulative')
            outdir = model.odir + '/plots/model_cumulative'
            if not path.exists(outdir): mkdir(outdir)
            MESSAGE("making plot for cumulative slip in %s " % outdir)

        processed_list = Parallel(n_jobs=num_cores) \
            (delayed(plot_model_shp) \
                 (H_shp[i],
                  cmap=model.plot_settings.colormap,
                  ncolor=ncolor,
                  contour=None,
                  crs=None,
                  log=False,
                  title=H_title[i],
                  bounds=bound_color,
                  outdir=outdir,
                  outfile=None,
                  shp_poly=model.external_shapefile_poly,
                  shp_line=model.external_shapefile_line,
                  shp_point=H_disp[i])

             for i in tqdm(np.arange(cstf.shape[0])))


    ###########################################################################
    # cumulative slip plot with observed/modeled displacements
    ###########################################################################

    if str2bool(model.plot_settings.disp_obs) or str2bool(model.plot_settings.disp_model):
        MESSAGE("plotting model/observed GPS displacements with cumulative slip")
        outdir = model.odir + '/plots/disp'
        if not path.exists(outdir): mkdir(outdir)

        H_GPS = {}
        H_GPS_up = {}
        if str2bool(model.plot_settings.disp_obs) and str2bool(model.plot_settings.disp_model):
            for i in np.arange(len(H_obs)):
                H_GPS[i] = list((H_obs[i], H_pred[i]))
                H_GPS_up[i] = list((H_obs_up[i], H_pred_up[i]))
        else:
            if str2bool(model.plot_settings.disp_obs):
                H_GPS = H_obs
                H_GPS_up = H_obs_up
            if str2bool(model.plot_settings.disp_model):
                H_GPS = H_pred
                H_GPS_up = H_pred_up

        # horizontal
        processed_list = Parallel(n_jobs=num_cores) \
            (delayed(plot_model_shp) \
                 (H_shp[i],
                  cmap=model.plot_settings.colormap,
                  ncolor=ncolor,
                  contour=None,
                  crs=None,
                  log=False,
                  title=H_title[i],
                  bounds=bound_color,
                  outdir=outdir,
                  outfile=None,
                  shp_poly=model.external_shapefile_poly,
                  shp_line=model.external_shapefile_line,
                  shp_point=H_disp[i],
#                  disp=[H_obs[i], H_pred[i]],
                  disp=H_GPS[i],
                  disp_scale=sda)

             for i in tqdm(np.arange(cstf.shape[0])))
        # up
        processed_list = Parallel(n_jobs=num_cores) \
            (delayed(plot_model_shp) \
                 (H_shp[i],
                  cmap=model.plot_settings.colormap,
                  ncolor=ncolor,
                  contour=None,
                  crs=None,
                  log=False,
                  title=H_title[i],
                  bounds=bound_color,
                  outdir=outdir,
                  outfile=None,
                  shp_poly=model.external_shapefile_poly,
                  shp_line=model.external_shapefile_line,
                  shp_point=H_disp[i],
                  #                  disp=[H_obs[i], H_pred[i]],
                  disp=H_GPS_up[i],
                  disp_scale=sda)

             for i in tqdm(np.arange(cstf.shape[0])))

    ###########################################################################
    # plot slip_dir_en
    ###########################################################################

    slip_dir_file = model.odir + '/info/slip_dir_en.dat'
    outdir = model.odir + '/plots/map'

    MESSAGE("plotting slip direction")

    plot_model_shp(H_shp[cstf.shape[0]-1],
              cmap=model.plot_settings.colormap,
              ncolor=ncolor,
              contour=model.plot_settings.ncontour,
              interpolate=model.plot_settings.interpolate.lower(),
              crs=None,
              log=False,
              title=H_title[i],
              bounds=bound_color,
              outdir=outdir,
              outfile='slip_dir_en.pdf',
              shp_poly=model.external_shapefile_poly,
              shp_line=model.external_shapefile_line,
              shp_point=H_disp[i],
              disp=slip_dir_file,
              disp_scale=10.)

    ###########################################################################
    # cumulative slip plot with res
    ###########################################################################
    if str2bool(model.plot_settings.disp_residual):
        MESSAGE("plotting residuals displacement with cumulative slip")

        outdir = model.odir + '/plots/disp_res'
        if not path.exists(outdir): mkdir(outdir)

        processed_list = Parallel(n_jobs=num_cores) \
            (delayed(plot_model_shp) \
                 (H_shp[i],
                  cmap=model.plot_settings.colormap,
                  ncolor=ncolor,
                  contour=None,
                  crs=None,
                  log=False,
                  title=H_title[i],
                  bounds=bound_color,
                  outdir=outdir,
                  outfile=None,
                  shp_poly=model.external_shapefile_poly,
                  shp_line=model.external_shapefile_line,
                  shp_point=H_disp[i],
                  disp=[H_res[i]],
                  disp_scale = 5*sda)

             for i in tqdm(np.arange(cstf.shape[0])))

    ###########################################################################
    # cumulative slip contour plot
    ###########################################################################
    if str2bool(model.plot_settings.ccum_slip):
        MESSAGE("plotting contour cumulative slip")
        outdir = model.odir + '/plots/model_cumulative_contour'
        if not path.exists(outdir): mkdir(outdir)

        processed_list = Parallel(n_jobs=num_cores) \
            (delayed(plot_model_shp) \
                 (H_shp[i],
                  cmap=model.plot_settings.colormap,
                  ncolor=ncolor,
                  contour=model.plot_settings.ncontour,
                  interpolate=model.plot_settings.interpolate.lower(),
                  crs=None,
                  log=False,
                  title=H_title[i],
                  bounds=bound_color,
                  outdir=outdir,
                  outfile=None,
                  shp_poly=model.external_shapefile_poly,
                  shp_line=model.external_shapefile_line,
                  shp_point=H_disp[i])
             for i in tqdm(np.arange(cstf.shape[0])))


    ###################################################################
    # PLOT MODEL RATE
    ###################################################################

    ###################################################################
    # VARIABLES FOR LOOP OF RATE MODELS
    ###################################################################

    mxslip = model.rate_slip_max * 1.1

    # color bounds
    if model.plot_settings.rate_bounds in ['None','NONE']:
        bound_color = [0., mxslip]
    else:
        # user provided bounds [min_slip, max_slip]
        lbounds = model.plot_settings.rate_bounds.replace('[','').replace(']','').split(',')
        # deal with upper bound
        if lbounds[1] in ['None','NONE']:
            ub = mxslip
        else:
            ub = float( lbounds[1] )
        # deal with lower bound
        # lower bound provided as percent
        if '%' in lbounds[0]:
            lb = float(lbounds[0].split('%')[0])/100. * mxslip
        else:
        # lower bound as true value
            lb = float( lbounds[0] )

        bound_color = [lb,ub]

    H_title = {}
    H_shp = {}
    H_ishp = {}
    stf = np.atleast_2d(np.genfromtxt(model.odir + '/stf/stf.dat', dtype=str))


    for i in np.arange(stf.shape[0]):
        H_title[i] = ("%s -- time step %04d \n day %6.1lf - %s %s\nMo=%s N.m/day" % (
        model.name, i, float(stf[i, 0]), stf[i, 2], stf[i, 3], stf[i, 1]))
        H_shp[i] = ("%s/shapefile/slip_rate/%04d_slip_rate.shp" % (model.odir, i))
        H_ishp[i] = ("%s/shapefile/i_slip_rate/%04d_slip_rate.shp" % (model.odir, i))


    ###########################################################################
    # model rate plot (linear scale)
    ###########################################################################
    if str2bool(model.plot_settings.rate):

        if interpolation:
            if not path.exists(model.odir + '/plots/i_model_rate'): mkdir(model.odir + '/plots/i_model_rate')
            outdir = model.odir + '/plots/i_model_rate'
            MESSAGE("making plot for interpolated slip rate in %s" % outdir)
        else:
            if not path.exists(model.odir + '/plots/model_rate'): mkdir(model.odir + '/plots/model_rate')
            outdir = model.odir + '/plots/model_rate'
            MESSAGE("making plot for slip rate in %s" % outdir)

        processed_list = Parallel(n_jobs=num_cores) \
            (delayed(plot_model_shp) \
                 (H_shp[i],
                  cmap=model.plot_settings.colormap,
                  ncolor=ncolor,
                  contour=None,
                  crs=None,
                  log=False,
                  title=H_title[i],
                  bounds=bound_color,
                  outdir=outdir,
                  outfile=None,
                  shp_poly=model.external_shapefile_poly,
                  shp_line=model.external_shapefile_line,
                  shp_point=H_disp[i + 1])
             for i in tqdm(np.arange(stf.shape[0])))

    ###########################################################################
    # model rate contour plot (log scale)
    ###########################################################################
    if str2bool(model.plot_settings.crate):

        # output directory
        if interpolation:
            if not path.exists(model.odir + '/plots/i_model_rate'): mkdir(model.odir + '/plots/i_model_rate_contour')
            outdir = model.odir + '/plots/i_model_rate_contour'
            MESSAGE("making plot for contour interpolated slip rate in %s" % outdir)
        else:
            if not path.exists(model.odir + '/plots/model_rate_contour'): mkdir(model.odir + '/plots/model_rate_contour')
            outdir = model.odir + '/plots/model_rate_contour'
            MESSAGE("making plot for contour slip rate in %s" % outdir)

        # color bounds
        bounds = np.linspace(np.log10(0.1), np.log10(get_max_scale(mxslip)), 1001)

        # color bounds
        if model.plot_settings.rate_bounds in ['None', 'NONE']:
            bound_color = [np.log10(0.1), np.log10(get_max_scale(mxslip))]
        else:
            # user provided bounds [min_slip, max_slip]
            lbounds = model.plot_settings.rate_bounds.replace('[', '').replace(']', '').split(',')
            # deal with upper bound
            if lbounds[1] in ['None', 'NONE']:
                ub = np.log10(get_max_scale(mxslip))
            else:
                ub = np.log10(get_max_scale(float(lbounds[1]) ))
            # deal with lower bound
            # lower bound provided as percent
            if '%' in lbounds[0]:
                lb = np.log10( float(lbounds[0].split('%')[0]) / 100. * mxslip )
            else:
                # lower bound as true value
                if float(lbounds[0])==0:
                    lb = np.log10(0.1)
                else:
                    lb = np.log10(float(lbounds[0]))

            bound_color = [lb, ub]
#        print(lb,ub)

        processed_list = Parallel(n_jobs=num_cores) \
            (delayed(plot_model_shp) \
                 (H_shp[i],
                  cmap=model.plot_settings.colormap,
                  ncolor=ncolor,
                  contour=model.plot_settings.ncontour,
                  interpolate=model.plot_settings.interpolate.lower(),
                  crs=None,
                  log=True,
                  title=H_title[i],
                  bounds=bound_color,
                  outdir=outdir,
                  outfile=None,
                  shp_poly=model.external_shapefile_poly,
                  shp_line=model.external_shapefile_line,
                  shp_point=H_disp[i + 1])
             for i in tqdm(np.arange(stf.shape[0])))
