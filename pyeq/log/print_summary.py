def print_summary( model ):

    ###################################################################
    # IMPORT
    ###################################################################
    
    import numpy as np
    import pyacs.lib.astrotime as at
    from shutil import copyfile
    import copy
    import os
    import shapefile
    from datetime import datetime,timedelta
    import time
    
    import pyacs
    from pyeq.lib import eq_disloc_3d as DL
    from pyacs.lib.vel_field import Velocity_Field as VF
    from pyacs.lib.gmtpoint import GMT_Point
    import pyacs.lib.utils
    import pyacs.lib.glinalg
    from pyacs.gts.Sgts import Sgts
    from pyacs.gts.Gts import Gts
    import pyacs.lib.coordinates
    import pyeq.lib.lib_inversion
    import pyeq.lib.green_tensor
    import pyeq.lib.geometry
#    import pyeq.lib.log.make_dir_pyeq_output
#    import pyeq.lib.log.model2shp_gmt

    import pyeq.message.message as MESSAGE
    import pyeq.message.verbose_message as VERBOSE
    import pyeq.message.error as ERROR
    import pyeq.message.warning as WARNING
    import pyeq.message.debug_message as DEBUG

    fsumn = model.odir+'/summary/sum.dat'


    fsum=open(fsumn,'w')
    fsum.write('#--------------------------------------------------------------------------\n')
    fsum.write('# User and System information \n')
    fsum.write('#--------------------------------------------------------------------------\n')
    
    fsum.write("    user                             : %s\n" % model.username) 
    fsum.write("    OS                               : %s\n" % model.os_version)
    fsum.write("    hostname                         : %s\n" % model.hostname)
    fsum.write("    directory                        : %s\n" % model.wdir )
    fsum.write("    number of cpu                    : %d\n" % model.n_cpu)   
    fsum.write("    number of threads                : %d\n" % model.n_thread )
    fsum.write("    available memory in Gb           : %.2f\n" % model.memory)
    fsum.write("    python version                   : %s\n" % model.python_version)
    fsum.write("    pyacs  version                   : %s\n" % model.pyacs_version)
    fsum.write("    pyeq   version                   : %s\n" % model.pyeq_version)

    
    fsum.write('#--------------------------------------------------------------------------\n')
    fsum.write('# Inversion speed and memory usage: \n')
    fsum.write('#--------------------------------------------------------------------------\n')
    fsum.write("    number of parameters estimated   : %d ( %d faults %d steps) \n" % (model.nparameters , model.nfaults, model.nstep ))
    fsum.write("    run start date & time            : %s \n" % time.strftime("%Y/%m/%d  %H:%M:%S", time.localtime(model.start_time) ) )
    end_time = time.time()
    fsum.write("    run end   date & time            : %s \n" % time.strftime("%Y/%m/%d  %H:%M:%S", time.localtime( end_time ) ) )
    duration = end_time - model.start_time
    fsum.write("    run duration                     : %.1lf s = %.1lf mn = %.1lf h \n" % ( duration , duration/60., duration/3600.) )
    fsum.write("    build obs. linear system (s)     : %.1lf \n" % model.time_build)
    fsum.write("    build regularization     (s)     : %.1lf \n" % model.time_regularization)
    fsum.write("    merging obs & regularization (s) : %.1lf \n" % model.time_merging_obs_regularization)
    fsum.write("    time nnls (s)                    : %.1lf \n" % model.time_inversion)
    fsum.write("    nnls algorithm                   : %s \n" % model.nnls)
    fsum.write("    total memory used in Gb          : %.1lf\n" %  (model.memory_usage) )
    fsum.write("    memory after  loading data       : %.1lf\n" %  (model.memory_before_ons) )
    fsum.write("    memory after  build N_obs        : %.1lf\n" %  (model.memory_before_rns ) )
    fsum.write("    memory after  build N_reg        : %.1lf\n" %  (model.memory_before_mns ) )
    fsum.write("    memory after  N_obs + N_reg      : %.1lf\n" %  (model.memory_before_inv ) )
    fsum.write("    Normal matrix in Gb              : %.1lf\n" %  (model.N_size ) )

    
    fsum.write('#--------------------------------------------------------------------------\n')
    fsum.write('# Inversion argument: \n')
    fsum.write('#--------------------------------------------------------------------------\n')
    fsum.write("    command line                     : %s \n" % (model.odir+'/info/command_line.dat'))
    fsum.write("    input npz file                   : %s \n" % model.input_npz)

    fsum.write('#--------------------------------------------------------------------------\n')
    fsum.write('# Geometry information: \n')
    fsum.write('#--------------------------------------------------------------------------\n')
    fsum.write("    number of subfaults              : %d \n" % model.nfaults)

    fsum.write('#--------------------------------------------------------------------------\n')
    fsum.write('# Regularization: \n')
    fsum.write('#--------------------------------------------------------------------------\n')
    fsum.write("    regularization                   : %s \n" % model.regularization)

    # case regularization covariance
    if model.regularization == 'covariance':
        fsum.write("    user input sigma                 : %s \n" % model.sigma)
        fsum.write("    pyeq input sigma                 : %s \n" % model.sigma_type)
        if 'time_variable' in model.sigma_type:
            fsum.write("    sigma_time min - max             : %.1f - %.1f mm/sqrt(day)\n" % ( np.min(model.sigma_time) , np.max(model.sigma_time) ) )
        else:
            fsum.write("    sigma_time                       : %.1f mm/sqrt(day)\n" % ( model.sigma_time ) )

        if 'fault_variable' in model.sigma_type:
            fsum.write("    sigma_fault min - max            : %.1f - %.1f mm/sqrt(day)\n" % ( np.min(model.sigma_fault) , np.max(model.sigma_fault) ) )
        else:
            fsum.write("    sigma_fault                      : %.1f mm/sqrt(day)\n" % ( model.sigma_fault ) )

        fsum.write("    dc (km)                          : %.1f \n" % model.dc)
        # TAU
        if model.tau == 0:
            fsum.write("    tau (days) no temporal smoothing : %.1f \n" % model.tau)
        else:
            fsum.write("    tau (days) temporal smoothing    : %.1f \n" % model.tau)
        fsum.write("    cm_type                          : %s \n" % model.cm_type )
        if model.cm_norm == 'd0':
            fsum.write("    cm_norm                          : normalized \n" )
        else:
            fsum.write("    cm_norm                          : not normalized \n" )

    if model.regularization == 'laplacian':
        if isinstance(model.lambda_spatial_smoothing,np.ndarray):
            fsum.write("    lambda spatial  smoothing        : %.1f - %.1f\n" % ( np.min(model.lambda_spatial_smoothing) , np.max(model.lambda_spatial_smoothing) ) )
        else:
            fsum.write("    lambda spatial  smoothing        : %.1f \n" % model.lambda_spatial_smoothing )
        fsum.write("    lambda temporal smoothing        : %.1f \n" % model.lambda_temporal_smoothing )
        fsum.write("    user input sigma                 : %s \n" % model.sigma)
        fsum.write("    pyeq input sigma                 : %s \n" % model.sigma_type)
        if 'time_variable' in model.sigma_type:
            fsum.write("    sigma_time min - max             : %.1f - %.1f mm/sqrt(day)\n" % ( np.min(model.sigma_time) , np.max(model.sigma_time) ) )
        else:
            fsum.write("    sigma_time                       : %.1f mm/sqrt(day)\n" % ( model.sigma_time ) )

        if 'fault_variable' in model.sigma_type:
            fsum.write("    sigma_fault min - max            : %.1f - %.1f mm/sqrt(day)\n" % ( np.min(model.sigma_fault) , np.max(model.sigma_fault) ) )
        else:
            fsum.write("    sigma_fault                      : %.1f mm/sqrt(day)\n" % ( model.sigma_fault ) )


    fsum.write('#--------------------------------------------------------------------------\n')
    fsum.write('# Rake information: \n')
    fsum.write('#--------------------------------------------------------------------------\n')
    fsum.write("    rake_type (Euler or constant)    : %s \n" % model.rake_type)
    fsum.write("    rake_value (Euler pole or value) : %s \n" % model.rake_value)
    fsum.write("    rake_constraint (0=fixed, float) : %s \n" % model.rake_constraint)
#    fsum.write("    max_slip (0=from Euler,else user): %s \n" % model.max_slip)

    fsum.write('#--------------------------------------------------------------------------\n')
    fsum.write('# Observation information: \n')
    fsum.write('#--------------------------------------------------------------------------\n')
    fsum.write("    time series are from             : %s \n" % model.dir_ts)
    fsum.write("    number of input GPS time series  : %d \n" %  model.n_site_input_ts )
    fsum.write("    number of input time series epoch: %d \n" %  model.n_date_input_ts )
    fsum.write("    number of GPS sites without data : %d \n" %  model.n_site_no_data )
#    fsum.write("    number of input Green functions  : %d \n" %  model.n_site_in_green )
    fsum.write("    number of GPS sites used         : %d \n" %  model.t_obs.shape[1] )
    fsum.write("    number of epochs used            : %d \n" %  model.t_obs.shape[0] )
    fsum.write("    uncertainties E & N              : %s \n" %  model.h_uncertainty )
    fsum.write("    up component                     : %s \n" %  model.up )
    fsum.write("    uncertainty up rescaled by       : %s \n" %  model.s_up )

    fsum.write('#--------------------------------------------------------------------------\n')
    fsum.write('# Green tensor information: \n')
    fsum.write('#--------------------------------------------------------------------------\n')
    fsum.write("    Green tensor used for inversion  : %s \n" %  str( model.green.shape ) )


    fsum.write('#--------------------------------------------------------------------------\n')
    fsum.write('# Date information: \n')
    fsum.write('#--------------------------------------------------------------------------\n')
    fsum.write("    date rounding option             : %s \n" % ( model.rounding ))
    fsum.write("    number of model dates            : %d \n" % ( model.np_model_date_s.shape[0] ))
    fsum.write("    number of model time steps       : %d \n" % ( model.np_model_step_duration_s.shape[0] ))
    fsum.write("    model start date                 : %s doy %03d %13.8lf\n" % \
               ( model.np_model_date_isoformat[0] , at.datetime2dayno( model.np_model_datetime[0] )[1] , at.datetime2decyear( model.np_model_datetime[0])))
    fsum.write("    model end  date                  : %s doy %03d %13.8lf\n" % \
               ( model.np_model_date_isoformat[-1] , at.datetime2dayno( model.np_model_datetime[-1] )[1] , at.datetime2decyear( model.np_model_datetime[-1])))
    fsum.write("    median,shortest,longest time step: %.2lf %.2lf %.2lf days\n" % \
               (model.median_time_step_days,model.shortest_time_step_days,model.longest_time_step_days))
        
    
    fsum.write('#--------------------------------------------------------------------------\n')
    fsum.write('#   Inversion results: \n')
    fsum.write('#--------------------------------------------------------------------------\n')
    fsum.write("    Cumulative Moment (N.m)          : %8.2E \n" % model.M0)
    fsum.write("    Equivalent Mw                    : %8.1f \n" % model.magnitude)

    i_mr_max = np.argmax( model.STF ) 
    
    fsum.write("    max moment rate (N.m / day )     : %8.2E #%05d %s-%s\n" % \
               ( model.STF[i_mr_max], i_mr_max, model.np_model_date_isoformat[i_mr_max], model.np_model_date_isoformat[i_mr_max+1]))

    magnitude_mr= 2./3.*(np.log10( model.STF[ i_mr_max ])-9.1)
    fsum.write("    Equivalent Mw / day              : %8.1f \n" % magnitude_mr)

    i_cs_max = np.unravel_index(np.argmax( model.NORM_CUMULATIVE_SLIP_PER_TIME_STEP), model.NORM_CUMULATIVE_SLIP_PER_TIME_STEP.shape)
    slip_max = model.NORM_CUMULATIVE_SLIP_PER_TIME_STEP[ i_cs_max ]
    
    model.cumulative_slip_max = slip_max
    model.idx_cumulative_slip_max = i_cs_max[1]
    
    fsum.write("    cumulative slip max  (mm)        : %8.1f fault #%d %s\n" % ( slip_max , i_cs_max[1] , model.np_model_date_isoformat[ i_cs_max[0] ] ) )


    i_is_max =  np.unravel_index(np.argmax (model.NORM_INC_SLIP_PER_TIME_STEP), model.NORM_INC_SLIP_PER_TIME_STEP.shape)
    inc_slip_max = model.NORM_INC_SLIP_PER_TIME_STEP[ i_is_max ]
    
    fsum.write("    incremental slip max  (mm)       : %8.1f fault #%04d time_step %04d %s -> %s\n" % \
               (inc_slip_max,i_is_max[1],i_is_max[0],model.np_model_date_isoformat[i_is_max[0]],model.np_model_date_isoformat[i_is_max[0]+1]))



    i_rs_max = lindex = np.unravel_index(np.argmax(model.NORM_RATE_SLIP_PER_TIME_STEP), model.NORM_RATE_SLIP_PER_TIME_STEP.shape)
    rate_slip_max = model.NORM_RATE_SLIP_PER_TIME_STEP[ i_rs_max ]
    
    model.rate_slip_max = rate_slip_max
    

    
    fsum.write("    slip rate max  (mm/day)          : %8.1f fault #%04d time_step %04d %s -> %s\n" % \
               (rate_slip_max,i_rs_max[1],i_rs_max[0],model.np_model_date_isoformat[i_rs_max[0]],model.np_model_date_isoformat[i_rs_max[0]+1]))
    
    
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
    
    # PRINT RESULTS ON SCREEN
    
    f = open(fsumn, "r")
    text = f.read()
    VERBOSE(text)
    f.close()


    return model 