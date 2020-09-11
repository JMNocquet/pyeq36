def print_stf( model ):
    """
    Print stf, inc_stf, cstf, fault_average slio as npy and text file
    """

    ###########################################################################
    # import
    ###########################################################################
    import numpy as np
    import pyacs.lib.utils
    import pyacs.lib.astrotime as at
    
    ###########################################################################
    # shear elastic modulus for moment calculation
    ###########################################################################
    mu=3.0E10
    
    ###########################################################################
    # STF, INC_STF, CSTF as npy
    ###########################################################################
    # area in meter-square
    model.AREA = model.sgeometry.tdis_area * 1.0E6
    model.SAREA = np.sum( model.AREA )
    
    # STF
    # norm of slip if variable rake
    model.NORM_RATE_SLIP_PER_TIME_STEP = np.sqrt( np.sum( model.RATE_SLIP_PER_TIME_STEP**2, axis=2) )
    model.STF = mu * 1.E-3 * np.sum( ( model.NORM_RATE_SLIP_PER_TIME_STEP * model.AREA ) , axis= 1)
    model.FAULT_AVERAGE_SLIP_RATE = model.STF / mu * 1.E3 / model.SAREA
    np.save( model.odir+'/npy/stf.npy' , model.STF )
    np.save( model.odir+'/npy/fault_average_slip_rate.npy' , model.FAULT_AVERAGE_SLIP_RATE )

    # INC_STF    
    model.NORM_INC_SLIP_PER_TIME_STEP = np.sqrt( np.sum( model.INC_SLIP_PER_TIME_STEP**2, axis=2) )
    model.INC_STF = mu * 1.E-3 * np.sum( ( model.NORM_INC_SLIP_PER_TIME_STEP * model.AREA ) , axis= 1)
    np.save( model.odir+'/npy/incremental_stf.npy' , model.INC_STF )

    # CSTF
    model.NORM_CUMULATIVE_SLIP_PER_TIME_STEP = np.sqrt( np.sum( model.CUMULATIVE_SLIP_PER_TIME_STEP**2, axis=2) )
    model.CSTF = mu * 1.E-3 * np.sum( ( model.NORM_CUMULATIVE_SLIP_PER_TIME_STEP * model.AREA ) , axis= 1)
    np.save( model.odir+'/npy/cstf.npy' , model.CSTF )

    ###########################################################################
    # STF, INC_STF, CSTF as text files
    ###########################################################################
    
    # STF
    model.DELTA_D_AND_STF = np.vstack( ( model.np_mid_model_delta_d , model.STF) ).T
    format = "%5.2lf  %.4E  %s"
    comment = ("dates are at the middle of model time steps. Decimal days since model date start: %s" % \
               ( model.np_model_date_isoformat[0] ) )

    pyacs.lib.utils.save_np_array_with_string( model.DELTA_D_AND_STF , \
                                               model.np_mid_model_isoformat ,\
                                              format, \
                                              model.odir+'/stf/stf.dat', \
                                              comment )

    comment = ("dates are at the middle of model time steps. Decimal days since model date start: %s. Displacement values in mm" % \
               ( model.np_model_date_isoformat[0] ) )
    model.DELTA_D_AND_FAULT_AVERAGE_SLIP_RATE = np.vstack( ( model.np_mid_model_delta_d , model.FAULT_AVERAGE_SLIP_RATE) ).T
    pyacs.lib.utils.save_np_array_with_string( model.DELTA_D_AND_FAULT_AVERAGE_SLIP_RATE , \
                                               model.np_mid_model_isoformat ,\
                                              format, \
                                              model.odir+'/stf/fault_average_slip_rate.dat', \
                                              comment )
    
    
    # INC_STF
    model.DELTA_D_AND_INC_STF = np.vstack( ( model.np_mid_model_delta_d , model.INC_STF) ).T
    format = "%5.2lf  %.4E  %s"
    comment = ("dates are the end date of model time steps. Decimal days since model date start: %s" % ( model.np_model_date_isoformat[0] ) )

    pyacs.lib.utils.save_np_array_with_string( model.DELTA_D_AND_INC_STF , \
                                               model.np_mid_model_isoformat ,\
                                              format, \
                                              model.odir+'/stf/incremental_stf.dat', \
                                              comment )

    # CSTF
    model.DELTA_D_AND_CSTF = np.vstack( ( model.np_model_delta_d , model.CSTF) ).T
    format = "%5.2lf  %.4E  %s"
    comment = ("Decimal days since model date start: %s" % ( at.seconds2datetime( model.np_model_date_s[0] ).isoformat(' ')))

    pyacs.lib.utils.save_np_array_with_string( model.DELTA_D_AND_CSTF , \
                                               model.np_model_date_isoformat ,\
                                              format, \
                                              model.odir+'/stf/cstf.dat', \
                                              comment )

    return model 
