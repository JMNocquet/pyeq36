###############################################################################
def print_results( model ):
        
###############################################################################

    """
    Print kinematic inversion results for the new pyeq >= 0.50.8
    """

    ###################################################################
    # IMPORT
    ###################################################################
    
    import numpy as np
    from shutil import copyfile
    import pickle

    import pyeq.log.make_dir_pyeq_output
    import pyeq.message.message as MESSAGE
    import pyeq.message.verbose_message as VERBOSE
    import pyeq.message.error as ERROR
    import pyeq.message.warning as WARNING
    import pyeq.message.debug_message as DEBUG


    ###################################################################
    # PARAMETERS
    ###################################################################

    if model.nconstant > 0:
        model.slip = model.parameters[ :-model.nconstant ]
        model.estimated_offsets = model.parameters[ -model.nconstant: ]
    else:
        model.slip = model.parameters
        model.estimated_offsets = None

    ###################################################################
    # OUTPUT PRINT_RESULTS DIRECTORY
    ###################################################################
    model.odir = model.name
    odir = model.odir
    MESSAGE("all inversion results in: %s" % model.odir )
    pyeq.log.make_dir_pyeq_output.make_dir_pyeq_output(model.odir)

    ###################################################################
    # PRINT CONF FILE
    ###################################################################
    VERBOSE("printing conf files: %s" % (model.odir+'/conf' ) )
    pyeq.log.print_conf(model)


    ###################################################################
    # SAVE MODEL AND OBSERVATION DATES
    ###################################################################
    VERBOSE("printing model and observation dates in: %s" % model.odir )
    model = pyeq.log.print_dates(model)

    ###################################################################
    # SAVE PYEQ_KINEMATICS COMMAND LINE IN INFO DIR
    ###################################################################

    VERBOSE("saving pyeq command line in %s" % model.odir+'/info/command_line.dat')
    
    fcmd=open(model.odir+'/info/command_line.dat','w')
    fcmd.write(model.cmd_line)
    fcmd.close()

    ###################################################################
    # PRINT ESTIMATED OFFSETS IN INFO DIR
    ###################################################################

    VERBOSE("printing estimated time series offsets in %s" % model.odir+'/info/ts_offset.dat' )
    model = pyeq.log.print_offset(model)

    ###########################################################################
    # SAVE GEOMETRY IN INFO DIR
    ###########################################################################

    VERBOSE("saving geometry file in %s/info/geometry.dat" % model.odir )
    model = pyeq.log.print_geometry(model)
    
    ###########################################################################
    # SAVE WARNING FILE
    ###########################################################################


    VERBOSE("saving warning file in %s/info/warning.dat" % model.odir)
    fwarning=open(model.odir+'/info/warning.dat','w')
    fwarning.write(model.warning)
    fwarning.close()
    
    ###########################################################################
    # NPY TENSORS
    # RATE_SLIP, DELTA_SLIP & SLIP
    # HANDLES ORIGIN TIME CONSTANT
    # HANDLES VARIABLES RAKE CASE
    # SAVE SOME NPY FILES in NPY DIR
    ###########################################################################

    VERBOSE("saving slip, t_obs, green tensors as npy in dir: %s" % (model.odir+'/npy') )
    
    # SOLUTION
    np.save(model.odir+'/npy/slip.npy',model.slip)
    
    # COPY INPUT NPZ
    # 13/01/2021
    # this takes a lot of space and is barely used
    # furthermore input_npz might no be avaialble if print results is run on a computer different
    # from the one used for the run
    #copyfile(model.input_npz, model.odir+'/npy/input.npz')

    # OBS TENSOR
    ###########################################################################
    np.save(model.odir+'/npy/t_obs.npy', model.t_obs)
    
    # GREEN TENSOR
    ###########################################################################
    np.save(model.odir+'/npy/green.npy', model.green)

    # INCREMENTAL SLIP, SLIP RATE
    ###########################################################################
    model = pyeq.log.print_sol_to_slip(model)

    # Geometry
    ###########################################################################
    np.save(  model.odir+'/npy/geometry.npy' , model.geometry )

    ###########################################################################
    # SAVE SLIP TIME SERIES AS TEXT FILE
    ###########################################################################

    VERBOSE("saving text files: rate, incremental and cumulative slip time series in %s/slip_time_series" % (model.odir) )
    model = pyeq.log.print_slip_time_series(model)

    ###########################################################################
    # STF, INC_STF, CSTF as npy and text file
    ###########################################################################

    VERBOSE("print stf files in %s and %s" % ( model.odir+'/npy' , model.odir+'/stf' ) )
    model = pyeq.log.print_stf(model)

    ###########################################################################
    # PRINT SLIP MODEL
    ###########################################################################

    VERBOSE("print slip models in %s" % ( model.odir+'/slip' ) )
    # no model here
    pyeq.log.print_slip_model(model)

    ###########################################################################
    # PRINT MODEL PREDICTED TIME SERIES
    ###########################################################################

    VERBOSE("print observed, modeled ans residual time series in %s" % ( model.odir+'/time_series' ))
    model = pyeq.log.print_modeled_time_series(model)
            
    ###########################################################################
    # OBS TIME SERIES REALLY USED IN THE INVERSION
    ###########################################################################

    model = pyeq.log.print_observed_time_series(model)

    ###########################################################################
    # RESIDUAL TIME SERIES
    ###########################################################################
    # no return here
    
    model = pyeq.log.print_residual_time_series(model)

    ###########################################################################
    # STATS
    ###########################################################################

    VERBOSE("print inversion statistics")
    model = pyeq.log.print_stats(model)

    ###########################################################################
    # MODEL/OBS/RESIDUAL DISPLACEMENT FIELDS FOR EACH OBSERVATION DATE
    ###########################################################################

    pyeq.log.print_displacement_fields(model)

    ###########################################################################
    # WRITE SUMMARY
    ###########################################################################

    model.M0 = model.CSTF[-1]
    model.magnitude = 2./3.*(np.log10( model.M0 )-9.1)

    model = pyeq.log.print_summary(model)

    ###################################################################
    # SAVE MODEL PCK (MODEL AS A PICKLE)
    ###################################################################
    VERBOSE("writting model with results attributes as pickle in %s" % ( odir+'/npy/model.mpck' ))
    ofile = open( odir+'/npy/model.mpck', 'wb')
    pickle.dump( model , ofile , pickle.HIGHEST_PROTOCOL)
    ofile.close()



   