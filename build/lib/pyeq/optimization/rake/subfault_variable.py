
def subfault_variable( model ):
    """
    Runs a non-linear optimization with variable rake.
    variable rake is implemented as a correction to the main slip azimuth determined
    from the model.rake_type and model.rake_value.
    The estimated correction is controled by model.rake_constraint.
    model.rake_contraint = 'X/Y/Z'.
    X = 0: no correctionn estimated
    X = 1: a single value correction is estimated
    X > 1: variable rake, X**2 is the number of points over the geometry provided. Slip azimuth will be interpolated
    over all subfaults. X = 3 is usually a good guess.
    Y : float for absolute value bounds in slip direction correction. For instance, 20 means that slip azimuth will
    be bounded by -20 to + 20 degrees.
    Z: is the maximum iteration for the optimization algorithm. 10 is usually a good guess to start.
    Rake variable is loggued in model.rake_variable_log attribute.
    """


    ###########################################################################
    # IMPORT
    ###########################################################################

    # GENERAL

    import sys
    import numpy as np
    from time import time
    import copy
    import logging

    # PYACS

    import pyeq.date
    import pyeq.regularization
    import pyeq.regularization.damping
    import pyacs.lib.astrotime as at
    import pyacs.lib.faultslip

    # PYEQ
    import pyeq.forward_model
    import pyeq.lib.objects
    import pyeq.log
    import pyeq.conf
    import pyeq.elastic_tensor
    import pyeq.lib.green_tensor
    # change 27/03/2021
    import pyeq.elastic_tensor.green_ds_ss_to_main_rake_conjugate
    import pyeq.gps_time_series
    import pyeq.obs_tensor.set_zero_at_first_obs
    import pyeq.optimization.wrapper.make_inversion

    import pyeq.message.message as MESSAGE
    import pyeq.message.verbose_message as VERBOSE
    import pyeq.message.error as ERROR
    import pyeq.message.warning as WARNING
    import pyeq.message.debug_message as DEBUG


    MESSAGE("SUBFAULT VARIABLE SLIP DIRECTION INVERSION" , level=1)

    ###########################################################################
    # DECIPHER model.rake_constraint
    ###########################################################################

    if model.rake_constraint.count('/') != 2:
        ERROR(("Could not decipher rake_constraint argument: %s " % model.rake_constraint),exit=True)

    try:
        [rv_ndelta, rv_bound, rv_maxiter ] = list( map(int , model.rake_constraint.split('/')))
    except:
        ERROR(("Could not decipher rake_constraint argument: %s " % model.rake_constraint),exit=True)

    MESSAGE("Number of points used to estimate rake/slip azimuth delta: %03d " % rv_ndelta)
    MESSAGE("Bounds for rake/slip azimuth delta: %03d " % rv_bound)
    MESSAGE("Maxiter for optimization: %03d " % rv_maxiter)


    ###########################################################################
    # ROUTINE FOR SUBFAULT-FIXED RAKE INVERSION
    ###########################################################################

    def mini_pyaks( mmodel, rake_subfault ):

        # import
        import copy

        # computes the Green tensor according to new RAKE

        mmodel.green, _ = pyeq.elastic_tensor.green_ds_ss_to_main_rake_conjugate(  \
            mmodel.green_ds_ss, mmodel.geometry, mmodel.sgeometry, 'vector', rake_subfault )

        # now these are firmly defined

        mmodel.normal_system_build = 'pyeq.forward_model.build5'
        (mmodel.N, mmodel.Nd) = pyeq.forward_model.build5(mmodel)

        mmodel.Nobs = copy.deepcopy( mmodel.N )
        mmodel.Ndobs = copy.deepcopy( mmodel.Nd )

        mmodel = pyeq.regularization.laplace.add_laplace_cons(mmodel)
        mmodel = pyeq.regularization.damping.add_damping_cons(mmodel)

        # computes solution

        if mmodel.regularization in ['covariance', 'laplacian']:
            mmodel.parameters, time_inversion = pyeq.optimization.wrapper.make_inversion.pyeq_nnls(mmodel.N, mmodel.Nd, mmodel.nnls,
                                                                                                  verbose=False)

        return mmodel

    ###########################################################################
    # COST FUNCTION
    ###########################################################################

    def cost_function( DELTA_SLIPAZGRD , model , az_points ):

        import pyacs.lib.faultslip
        from scipy.interpolate import griddata

        VERBOSE("Computing cost function for variable rake")

        # make a deep copy of model
        mmodel = copy.deepcopy( model )
        VERBOSE( ("delta_slip_azimuth grid min - max - mean : %10.2lf %10.2lf %10.2lf" % (np.min(DELTA_SLIPAZGRD),np.max(DELTA_SLIPAZGRD),np.mean(DELTA_SLIPAZGRD)) ) )

        # interpolate slip azimuth from points located at az_points with values DELTA_SLIPAZGRD

        DELTA_SLIPAZ =  griddata( az_points, DELTA_SLIPAZGRD, ( mmodel.sgeometry.centroid_long, mmodel.sgeometry.centroid_lat ), method='linear')

        VERBOSE( ("delta_slip_azimuth subfaults min - max - mean : %10.2lf %10.2lf %10.2lf" % (np.min(DELTA_SLIPAZ),np.max(DELTA_SLIPAZ),np.mean(DELTA_SLIPAZ)) ) )

        # update rake
        updated_slipdir =  mmodel.slip_az + DELTA_SLIPAZ
        if mmodel.rake_type.lower() == 'euler':
            np_motion_type  =  np.array([mmodel.rake_value.split('/')[-1]]*DELTA_SLIPAZ.shape[0])
        else:
            if mmodel.rake_type in ['fixed','constant']:
                if float (mmodel.rake_value ) >0: np_motion_type = np.array(['inverse']*DELTA_SLIPAZ.shape[0])
                else: np_motion_type = np.array(['normal']*DELTA_SLIPAZ.shape[0])


        VERBOSE(("original rake min - max - mean : %10.2lf %10.2lf %10.2lf" % (np.min(mmodel.rake_subfault),np.max(mmodel.rake_subfault),np.mean(mmodel.rake_subfault))))
        RAKE = pyacs.lib.faultslip.rake_from_slip_az( mmodel.sgeometry.strike, mmodel.sgeometry.dip, updated_slipdir, np_motion_type)
        VERBOSE(("test rake     min - max - mean : %10.2lf %10.2lf %10.2lf" % (np.min(RAKE),np.max(RAKE),np.mean(RAKE))))

        # min pyaks - runs inversion

        mmodel = mini_pyaks( mmodel , RAKE )

        # observation + regularization

        C = np.copy( mmodel.N ) @ mmodel.parameters
        C = C.flatten()

        # cost inversion
        O_C = mmodel.Nd - C
        cost = np.sum( O_C**2 )
        VERBOSE("cost Normal system  %.12E" % cost )

        # adds cost of damping on correction
        cost = cost + (mmodel.rake_subfault.shape[0] / az_points.shape[0] / rv_bound)**2 * np.sum( DELTA_SLIPAZGRD**2 )
        #+ 1.E2 * cost_smoothing_rake

        VERBOSE("cost Normal System + cost delta azimuth %.12E" % cost )

        return cost
    # end cost function

    ###########################################################################
    # MAIN
    ###########################################################################

    from scipy.optimize import minimize
    from scipy.optimize import Bounds
    from scipy.interpolate import griddata

    # get spatial DLO
    #VERBOSE("Making Discrete Laplace Operator for the triangular mesh")
    #dlos = pyeq.regularization.laplace.make_dlo_trimesh(model.geometry, stencil=4, verbose=model.verbose)
    # normalize DLO
    #VERBOSE("Normalizing Discrete Laplace Operator for the triangular mesh")
    #ndlos = np.copy( pyeq.regularization.laplace.normalize_dlo_trimesh(dlos, 2).todense() )

    MESSAGE("MAKING GRID FOR SLIP AZIMUTH")

    min_lon = np.min( model.sgeometry.centroid_long )
    max_lon = np.max( model.sgeometry.centroid_long )

    min_lat = np.min( model.sgeometry.centroid_lat )
    max_lat = np.max( model.sgeometry.centroid_lat )

    np_lon = np.linspace( min_lon, max_lon, rv_ndelta )
    np_lat = np.linspace( min_lat, max_lat, rv_ndelta )

    (lon,lat) = np.meshgrid(np_lon,np_lat)
    az_points = np.array([ lon.flatten(), lat.flatten()]).T

    # BOUNDS
    # TODO bounds on rake delta
    lb = np.zeros( az_points.shape[0] ) - rv_bound
    ub = np.zeros( az_points.shape[0] ) + rv_bound
    bounds = Bounds( lb,ub )

    MESSAGE("MAKING OPTIMIZATION USING L-BFGS-B")

    X0 = np.zeros( az_points.shape[0] )

    #from scipy.optimize import dual_annealing
    #RESULT = dual_annealing( cost, bounds=list(zip(lb, ub)) , args=(model, ndlos , az_points )),

    RESULT = minimize(cost_function, X0, args=(model, az_points ),
                          method='L-BFGS-B' , bounds=bounds , options={'maxiter':rv_maxiter, 'disp':True} )


    MESSAGE("RESULTS FROM OPTIMIZATION" , level=1)
    print(RESULT)

    best_slip_dir = griddata(az_points, RESULT.x , (model.sgeometry.centroid_long, model.sgeometry.centroid_lat),
             method='linear')

    cost_original  = cost_function( X0 , model , az_points )
    cost_optimized = cost_function( RESULT.x , model , az_points )

    cost_gain = (cost_original - cost_optimized)/cost_original * 100.

    MESSAGE("cost initial   rake: %15.10E" % cost_original)
    MESSAGE("cost optimized rake: %15.10E" % cost_optimized)
    MESSAGE("gain (percent)     : %15.4lf" % cost_gain)



    # info on slip azimuth change
    from scipy.interpolate import griddata
    DELTA_SLIPAZ = griddata(az_points, RESULT.x, (model.sgeometry.centroid_long, model.sgeometry.centroid_lat),
                            method='linear')

    MESSAGE("min max mean slip azimuth change: %6.2lf %6.2lf %6.2lf" % (np.min(DELTA_SLIPAZ),np.max(DELTA_SLIPAZ), np.mean(DELTA_SLIPAZ)))

    # info on rake change

    updated_slipdir = model.slip_az + DELTA_SLIPAZ
    # TODO handle case rake_type is constant; currently Euler only implemented
    np_motion_type = np.array([model.rake_value.split('/')[-1]] * DELTA_SLIPAZ.shape[0])

    new_rake = pyacs.lib.faultslip.rake_from_slip_az(model.sgeometry.strike, model.sgeometry.dip, updated_slipdir,
                                             np_motion_type)

    delta_rake = new_rake - model.rake_subfault
    MESSAGE("min max mean rake change        : %6.2lf %6.2lf %6.2lf" % (np.min(delta_rake),np.max(delta_rake), np.mean(delta_rake)))


    # chi2 obs

    mmodel = copy.deepcopy(model)
    mmodel = mini_pyaks( mmodel , mmodel.rake_subfault )

    cost_obs_original             = np.sum( (mmodel.Nobs @ mmodel.parameters - mmodel.Ndobs)**2 )
    cost_obs_regularized_original = np.sum( ( np.asarray( mmodel.N )    @ mmodel.parameters - mmodel.Nd)**2 )


    mmodel = copy.deepcopy(model)

    mmodel = mini_pyaks( mmodel , new_rake )
    mmodel.rake_subfault =  new_rake
    mmodel.slip_az =  updated_slipdir
    mmodel.slip_dir_en = np.array( [np.sin( np.radians( mmodel.slip_az )) ,  np.cos( np.radians( mmodel.slip_az ) ) ] ).T

    cost_obs_optimized = np.sum( (mmodel.Nobs @ mmodel.parameters - mmodel.Ndobs)**2 )
    cost_obs_gain = (cost_obs_original - cost_obs_optimized)/cost_obs_original * 100.

    cost_obs_regularized_optimized = np.sum( ( np.asarray( mmodel.N ) @ mmodel.parameters - mmodel.Nd)**2 )
    cost_obs_regularized_gain = (cost_obs_regularized_original - cost_obs_regularized_optimized)/cost_obs_regularized_original * 100.


    mmodel.log_variable_rake = ''

    mmodel.log_variable_rake += ("rake_type           : %s \n" %  mmodel.rake_type)
    mmodel.log_variable_rake += ("rake_value          : %s \n" %  mmodel.rake_value)
    mmodel.log_variable_rake += ("rake_constraint     : %s \n" %  mmodel.rake_constraint)

    for i in np.arange( az_points.shape[0] ):
        mmodel.log_variable_rake += ("az_point #%02d  : %10.4lf  %10.4lf  %10.2lf \n" % ( i, az_points[i,0], az_points[i,1] , RESULT.x[i]))

    mmodel.log_variable_rake += ("delta slip azimuth : min %5.2lf max %5.2lf mean %5.2lf  \n" %  ( np.min( DELTA_SLIPAZ ) , np.max( DELTA_SLIPAZ ) , np.mean( DELTA_SLIPAZ ) ) )
    mmodel.log_variable_rake += ("delta rake         : min %5.2lf max %5.2lf mean %5.2lf  \n" %  ( np.min( delta_rake ) , np.max( delta_rake ) , np.mean( delta_rake ) ) )

    mmodel.log_variable_rake += ("mean slip azimuth  : old %5.2lf optimized %5.2lf \n" %  ( np.mean( model.slip_az ) , np.mean( mmodel.slip_az ) ) )
    mmodel.log_variable_rake += ("mean rake          : old %5.2lf optimized %5.2lf \n" %  ( np.mean( model.rake_subfault ) , np.mean( mmodel.rake_subfault ) ) )

    mmodel.log_variable_rake += ("optimization msg   : %s \n" % RESULT.message)
    mmodel.log_variable_rake += ("opt.     success   : %s \n" % str(RESULT.success))

    mmodel.log_variable_rake += ("cost initial   rake: %15.10E \n" % cost_original)
    mmodel.log_variable_rake += ("cost optimized rake: %15.10E \n" % cost_optimized)
    mmodel.log_variable_rake += ("gain (percent)     : %15.4lf \n" % cost_gain)

    mmodel.log_variable_rake += ("cost obs initial  : %15.10E \n" % cost_obs_original)
    mmodel.log_variable_rake += ("cost obs optimized: %15.10E \n" % cost_obs_optimized)
    mmodel.log_variable_rake += ("gain obs (percent): %15.4lf \n" % cost_obs_gain)

    mmodel.log_variable_rake += ("cost obs regularized initial  : %15.10E \n" % cost_obs_regularized_original)
    mmodel.log_variable_rake += ("cost obs regularized optimized: %15.10E \n" % cost_obs_regularized_optimized)
    mmodel.log_variable_rake += ("gain obs regularized (percent): %15.4lf \n" % cost_obs_regularized_gain)

    for line in mmodel.log_variable_rake.split('\n'):
        MESSAGE(line)

    # remove the normal matrix to save space

    try:
        delattr(mmodel, 'Nobs')
        delattr(mmodel, 'Ndobs')
    except:
        ERROR("deleting normal observation system attributes")

    return mmodel