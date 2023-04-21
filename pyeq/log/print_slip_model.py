def print_slip_model( model ):
    """
    Print slip model as text files
    """

    ###########################################################################
    # IMPORT
    ###########################################################################

    import numpy as np
    import pyacs.lib.astrotime as at

    import pyeq.message.message as MESSAGE
    import pyeq.message.verbose_message as VERBOSE
    import pyeq.message.error as ERROR
    import pyeq.message.warning as WARNING
    import pyeq.message.debug_message as DEBUG


    ###########################################################################
    # INITIALIZE COORDINATES FROM GEOMETRY FOR PRINTING SLIP
    ###########################################################################

    RATE_SLIP_TIME_STEP = np.zeros(( model.nfaults, 2 + model.np_rake.shape[0] ))
    RATE_SLIP_TIME_STEP[:,0] = model.sgeometry.centroid_long
    RATE_SLIP_TIME_STEP[:,1] = model.sgeometry.centroid_lat

    INC_SLIP_TIME_STEP = np.zeros(( model.nfaults, 2 + model.np_rake.shape[0] ))
    INC_SLIP_TIME_STEP[:,0] = model.sgeometry.centroid_long
    INC_SLIP_TIME_STEP[:,1] = model.sgeometry.centroid_lat


    CUMULATIVE_SLIP_TIME_STEP=np.zeros(( model.nfaults, 2 + model.np_rake.shape[0] ))
    CUMULATIVE_SLIP_TIME_STEP[:,0] = model.sgeometry.centroid_long
    CUMULATIVE_SLIP_TIME_STEP[:,1] = model.sgeometry.centroid_lat

    
    ###########################################################################
    # LOOP ON MODEL TIME STEP FOR SLIP RATE & INCREMENTAL SLIP
    ###########################################################################
    for i in np.arange( model.np_model_step_duration_s.shape[0] ):
        DEBUG("calculating incremental slip, slip rate for model time step: %d" % (i) )
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
        RATE_SLIP_TIME_STEP[ : , 2: ]  =  model.RATE_SLIP_PER_TIME_STEP[i].reshape( model.nfaults, -1 )
        fname=model.odir+"/slip/rate/"+("%04d_slip_rate.dat" % ( i ))
        # 1 rake
        if RATE_SLIP_TIME_STEP.shape[1] == 3:
            RATE_SLIP_TIME_STEP[:,2] = RATE_SLIP_TIME_STEP[:,2] 
#            np.savetxt(fname, RATE_SLIP_TIME_STEP, fmt="%10.5lf %10.5lf %10.3lf", header=date_info)
            np.savetxt(fname, RATE_SLIP_TIME_STEP, fmt="%10.5lf %10.5lf %10.3G", header=date_info)
        # 2 rake
        else:
            np.savetxt(fname, RATE_SLIP_TIME_STEP, fmt="%10.5lf %10.5lf %10.3lf %10.3lf ", header=date_info)
            
        # INCREMENTAL SLIP
        INC_SLIP_TIME_STEP[:,2:] = model.INC_SLIP_PER_TIME_STEP[i].reshape( model.nfaults, -1 )
        fname=model.odir+"/slip/incremental/"+("%04d_delta_slip.dat" % ( i ))
        # 1 rake
        if INC_SLIP_TIME_STEP.shape[1] == 3:
            np.savetxt(fname, INC_SLIP_TIME_STEP, fmt="%10.5lf %10.5lf %10.3lf", header=date_info)
        # 2 rake
        else:
            np.savetxt(fname, INC_SLIP_TIME_STEP, fmt="%10.5lf %10.5lf %10.3lf %10.3lf ", header=date_info)
        
    # LOOP ON MODEL DATES FOR THE CUMULATIVE SLIP
    for i in np.arange( model.np_model_date_s.shape[0] ):

        isodate = model.np_model_date_isoformat[i]

        DEBUG("writing cumulative slip for model time : %04d %s" % (i ,  isodate ) )

        date_info=("model date  %04d %s" % (i ,  isodate ) )

        # CUMULATIVE SLIP
        CUMULATIVE_SLIP_TIME_STEP[:,2:] = model.CUMULATIVE_SLIP_PER_TIME_STEP[i].reshape( model.nfaults, -1 )
        fname=model.odir+"/slip/cumulative/"+("%04d_cumulative_slip.dat" % (i))
        # 1 rake
        if CUMULATIVE_SLIP_TIME_STEP.shape[1] == 3:
            np.savetxt(fname, CUMULATIVE_SLIP_TIME_STEP, fmt="%10.5lf %10.5lf %10.3lf", header=date_info)
        # 2 rake
        else:
            np.savetxt(fname, CUMULATIVE_SLIP_TIME_STEP, fmt="%10.5lf %10.5lf %10.3lf %10.3lf ", header=date_info)
