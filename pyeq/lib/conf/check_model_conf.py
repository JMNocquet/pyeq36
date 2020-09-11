def check_model_conf( model ):
    """
    Check all model attributes.
    """
    
    import pyacs.lib.utils
    import str2bool
    from colors import red
    import sys
    import numpy as np
    
    ###############################################################################
    # BEHAVIOUR CONTROL 
    ###############################################################################
    
    # VERBOSE
    model.verbose = str2bool.str2bool( model.verbose )
    # DEBUG
    model.debug = str2bool.str2bool( model.debug )
    # NO_OPT
    model.no_opt = str2bool.str2bool( model.no_opt )
    # SAVE
    model.save = str2bool.str2bool( model.save)
    # PLOT
    model.plot = str2bool.str2bool( model.plot)
    # TAR
    model.tar = str2bool.str2bool( model.tar)

    ###############################################################################
    # GEOMETRY
    ###############################################################################
    if model.geometry_range_lon is not None:
        model.geometry_range_lon = range_lon = pyacs.lib.utils.str2list_float( model.geometry_range_lon )

    if model.geometry_range_lat is not None:
        model.geometry_range_lat = range_lat = pyacs.lib.utils.str2list_float( model.geometry_range_lat )

    if model.geometry_range_depth is not None:
        model.geometry_range_depth = range_lat = pyacs.lib.utils.str2list_float( model.geometry_range_depth )

    if model.geometry_remove_idx is not None:
        model.geometry_remove_idx = np.array( pyacs.lib.utils.str2list_float( model.geometry_remove_idx ) , dtype=int ) 
    
    ###############################################################################
    # DATA
    ###############################################################################
    # UP
    if model.up == 'YES':model.up = True
    if model.up == 'NO' :model.up = False
    
    # OFFSET
    try:
        model.offset = float( model.offset )
    except:
        print("[ERROR]: offset must be float.")
        import sys
        sys.exit()

    # INTERSEISMIC
    try:
        model.interseismic = float( model.interseismic )
    except:
        print("[ERROR]: interseismic must be float.")
        import sys
        sys.exit()

    # EXCLUDE_GPS
    if model.lexclude_gps is not None:
        if type( model.lexclude_gps ) != list:
            model.lexclude_gps = model.lexclude_gps.split()
    
    # CONVERT TO FLOAT
    lopt_float = [ 's_up','s_h','dc','tau','rake_constraint' ]
    
    for option in lopt_float:
        model.__dict__[ option ] = float( model.__dict__[ option ] )

    ###############################################################################
    # ALGORITHM
    ###############################################################################
    # CONVERT TO INT
    lopt_int = [ 'build' ]
    
    for option in lopt_int:
        model.__dict__[ option ] = float( model.__dict__[ option ] )

    ###############################################################################
    # MANDATORY ARGUMENTS
    ###############################################################################

    mandatory_missing = False

    lmandatory_options = [ \
                          # INPUT
                          'input_npz' ,\
                          # DATA
                           'dir_ts' ,'up','s_up','s_h','offset',\
                          # REGULARIZATION 
                           'dc','sigma','cm_type','cm_norm','tau','m0',\
                          # RAKE
                          'rake_type','rake_value',\
                          # DATES
                          'dates','rounding',\
                          #ALGO
                          'nnls','build',\
                          # GENERAL
                          'verbose','debug','no_opt','save'
                           ]

    for option in lmandatory_options:
        option_value = getattr(  model , option )
        if option_value is None:
            print( red("[PYEQ ERROR] option not understood or missing: %s" % option ) )
            mandatory_missing = True
    
    if mandatory_missing:
        print(red("[PYEQ ERROR]: STOP."))
        sys.exit()

    return model