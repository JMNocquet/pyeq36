def print_sol_to_slip( model ):
    """
    Converts inversion solution into various slip information and print them as npy files
    """
    
    # import
    import numpy as np
    
    # initialization
    np_rake = model.np_rake
    nfaults = model.green.shape[1]
    n_constant = model.nconstant
    
    SLIP = model.slip
    
    # TENSOR RATE_SLIP_PER_TIME_STEP
    #if n_constant != 0:
    #    model.RATE_SLIP_PER_TIME_STEP = SLIP[:-n_constant].reshape( -1, nfaults , np_rake.shape[0]  )
    #else:
    model.RATE_SLIP_PER_TIME_STEP = SLIP.reshape( -1, nfaults , np_rake.shape[0] )

    # TENSOR INC_SLIP_PER_TIME_STEP
    model.INC_SLIP_PER_TIME_STEP = ( model.RATE_SLIP_PER_TIME_STEP.T * model.np_model_step_duration_days ).T
    
    # TENSOR CUMULATIVE_SLIP_PER_TIME_STEP
    model.CUMULATIVE_SLIP_PER_TIME_STEP = np.zeros( ( model.INC_SLIP_PER_TIME_STEP.shape[0]+1, model.INC_SLIP_PER_TIME_STEP.shape[1], model.INC_SLIP_PER_TIME_STEP.shape[2] ) )
    model.CUMULATIVE_SLIP_PER_TIME_STEP[1:,:,:] = np.cumsum( model.INC_SLIP_PER_TIME_STEP , axis=0 )
    
    print("-- saving slip_rate.npy [index_date,fault,rake] tensor associated slip_rate_dates.npy (datetime format) in %s " % (model.odir+'/npy') )

    print("-- slip_rate.npy shape " , model.RATE_SLIP_PER_TIME_STEP.shape )
    np.save( model.odir+'/npy/slip_rate.npy' , model.RATE_SLIP_PER_TIME_STEP )
    np.save(  model.odir+'/npy/slip_rate_datetime.npy' , model.np_mid_model_datetime )
    np.save(  model.odir+'/npy/slip_rate_delta_d.npy' , model.np_mid_model_delta_d )

    print("-- incremental_slip.npy shape " , model.INC_SLIP_PER_TIME_STEP.shape )
    np.save( model.odir+'/npy/incremental_slip.npy' , model.INC_SLIP_PER_TIME_STEP )
    np.save(  model.odir+'/npy/incremental_slip_datetime.npy' , model.np_mid_model_datetime )
    np.save(  model.odir+'/npy/incremental_slip_delta_d.npy' , model.np_mid_model_delta_d )

    print("-- cumulative_slip.npy shape " , model.CUMULATIVE_SLIP_PER_TIME_STEP.shape )
    np.save( model.odir+'/npy/cumulative_slip.npy' , model.CUMULATIVE_SLIP_PER_TIME_STEP )
    np.save(  model.odir+'/npy/cumulative_slip_datetime.npy' , model.np_model_datetime )
    np.save(  model.odir+'/npy/cumulative_slip_delta_d.npy' ,  model.np_model_delta_d )

    return model