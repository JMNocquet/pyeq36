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
    import pyacs.lib.astrotime as at
    from shutil import copyfile
    import copy
    import os
    import shapefile
    from datetime import datetime,timedelta
    import time
    import matplotlib.pyplot as plt
    import pickle
    
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
    import pyeq.lib.log.make_dir_pyeq_output


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
    print("-- printing inversion results in: %s" % model.odir )
    pyeq.lib.log.make_dir_pyeq_output.make_dir_pyeq_output( model.odir )

    ###################################################################
    # PRINT CONF FILE
    ###################################################################
    print("-- printing conf files: %s" % (model.odir+'/conf' ) )
    pyeq.lib.log.print_conf( model )


    ###################################################################
    # SAVE MODEL AND OBSERVATION DATES
    ###################################################################
    print("-- printing model and observation dates in: %s" % model.odir )
    model = pyeq.lib.log.print_dates( model )


    ###################################################################
    # SAVE PYEQ_KINEMATICS COMMAND LINE IN INFO DIR
    ###################################################################

    print('-- saving pyeq command line in ', model.odir+'/info/command_line.dat')
    
    fcmd=open(model.odir+'/info/command_line.dat','w')
    fcmd.write(model.cmd_line)
    fcmd.close()

    ###################################################################
    # PRINT ESTIMATED OFFSETS IN INFO DIR
    ###################################################################

    print('-- printing estimated time series offsets in ' , model.odir+'/info/ts_offset.dat' )
    model = pyeq.lib.log.print_offset( model )

    ###########################################################################
    # SAVE GEOMETRY IN INFO DIR
    ###########################################################################

    print("-- saving geometry file in %s/info/geometry.dat" % model.odir )
    model = pyeq.lib.log.print_geometry( model )
    
    ###########################################################################
    # SAVE WARNING FILE
    ###########################################################################
    
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

    print("-- saving slip, input_npz, t_obs, green tensors as npy in dir: %s" % (model.odir+'/npy') )
    
    # SOLUTION
    np.save(model.odir+'/npy/slip.npy',model.slip)
    
    # COPY INPUT NPZ
    copyfile(model.input_npz, model.odir+'/npy/input.npz')

    # OBS TENSOR
    ###########################################################################
    np.save(model.odir+'/npy/t_obs.npy', model.t_obs)
    
    # GREEN TENSOR
    ###########################################################################
    np.save(model.odir+'/npy/green.npy', model.green)

    # INCREMENTAL SLIP, SLIP RATE
    ###########################################################################
    model = pyeq.lib.log.print_sol_to_slip( model )


    # Geometry
    ###########################################################################
    np.save(  model.odir+'/npy/geometry.npy' ,  model.geometry )

    ###########################################################################
    # SAVE SLIP TIME SERIES AS TEXT FILE
    ###########################################################################

    print("-- saving text files: rate, incremental and cumulative slip time series in %s/slip_time_series" % (model.odir) )
    model = pyeq.lib.log.print_slip_time_series( model )


    ###########################################################################
    # STF, INC_STF, CSTF as npy and text file
    ###########################################################################

    print("-- print stf files in %s and %s" % ( model.odir+'/npy' , model.odir+'/stf' ) )    
    model = pyeq.lib.log.print_stf( model )


    ###########################################################################
    # PRINT SLIP MODEL
    ###########################################################################

    print("-- print slip models in %s" % ( model.odir+'/slip' ) )    
    # no model here
    pyeq.lib.log.print_slip_model( model )

    ###########################################################################
    # PRINT MODEL PREDICTED TIME SERIES
    ###########################################################################
    
    model = pyeq.lib.log.print_modeled_time_series( model )
            
    ###########################################################################
    # OBS TIME SERIES REALLY USED IN THE INVERSION
    ###########################################################################

    model = pyeq.lib.log.print_observed_time_series( model )

    ###########################################################################
    # RESIDUAL TIME SERIES
    ###########################################################################
    # no return here
    
    model = pyeq.lib.log.print_residual_time_series( model )

    ###########################################################################
    # STATS
    ###########################################################################

    model = pyeq.lib.log.print_stats( model )

    ###########################################################################
    # MODEL/OBS/RESIDUAL DISPLACEMENT FIELDS FOR EACH OBSERVATION DATE
    ###########################################################################

    pyeq.lib.log.print_displacement_fields( model )


    ###########################################################################
    # WRITE SUMMARY
    ###########################################################################

    model.M0 = model.CSTF[-1]
    model.magnitude = 2./3.*(np.log10( model.M0 )-9.1)

    model = pyeq.lib.log.print_summary( model )
    
    print("-- all results written in %s " % model.odir )

    ###################################################################
    # SAVE MODEL PCK (MODEL AS A PICKLE)
    ###################################################################
    print("-- writting model as pickle in %s" % ( odir+'/npy/model.pck' ))
    ofile = open( odir+'/npy/model.pck', 'wb') 
    pickle.dump( model , ofile , pickle.HIGHEST_PROTOCOL)
    ofile.close()



   