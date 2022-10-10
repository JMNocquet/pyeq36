def check_model_conf( model ):
    """
    Check all model attributes.
    """
    
    import pyacs.lib.utils
    from str2bool import str2bool
    import sys
    import numpy as np

    import pyeq.message.message as MESSAGE
    import pyeq.message.verbose_message as VERBOSE
    import pyeq.message.error as ERROR
    import pyeq.message.warning as WARNING
    import pyeq.message.debug_message as DEBUG

    ###############################################################################
    # BEHAVIOUR CONTROL (BOOLEAN OPTIONS)
    ###############################################################################

    lboolean_opt = ['verbose','debug','no_opt','save','plot','tar','print_result', 'norm_resolution']

    for boolean_opt in lboolean_opt:
        if boolean_opt not in model.__dict__:
            WARNING("option %s not provided. Set to False." % boolean_opt)
            model.__dict__[boolean_opt] = False
        else:
            if isinstance(model.__dict__[boolean_opt],str):
                model.__dict__[boolean_opt] = str2bool(model.__dict__[boolean_opt])

        if not isinstance(model.__dict__[boolean_opt], bool):
            ERROR("option %s must be boolean" % boolean_opt,exit=True)


    ###############################################################################
    # BEHAVIOUR CONTROL (None as default)
    ###############################################################################

    lnone_opt = ['mpck','name']

    for none_opt in lnone_opt:
        if none_opt not in model.__dict__:
            WARNING("option %s not provided. Set to None." % none_opt)
            model.__dict__[none_opt] = None


    ###############################################################################
    # GEOMETRY
    ###############################################################################

    lopt_geometry = ['geometry_range_lon','geometry_range_lat','geometry_range_depth','geometry_remove_idx']
    for opt_geometry in lopt_geometry:
        if opt_geometry not in model.__dict__:
            WARNING("option %s not provided. Set to None" % opt_geometry)
            model.__dict__[opt_geometry] = None


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
    if 'up' not in model.__dict__:
        WARNING("option up not provided. Set to False")
        model.__dict__['up'] = False
    if isinstance(model.__dict__['up'], str):
        model.__dict__['up'] = str2bool(model.__dict__['up'])

    # OFFSET
    try:
        model.offset = float( model.offset )
    except:
        ERROR("option offset must be float",exit=True)

    # INTERSEISMIC
    try:
        model.interseismic = float( model.interseismic )
    except:
        ERROR("option interseismic must be float",exit=True)
        import sys
        sys.exit()

    # EXCLUDE_GPS
    if 'lexclude_gps' not in model.__dict__:
        model.lexclude_gps = []
    else:
        if type( model.lexclude_gps ) != list:
            model.lexclude_gps = model.lexclude_gps.split()

    # UNWEIGHT_GPS
    if 'lunweight_gps' not in model.__dict__:
        model.lunweight_gps = []
    else:
        if type( model.lunweight_gps ) != list:
            model.lunweight_gps = model.lunweight_gps.split()


    ###############################################################################
    # CONVERT TO FLOAT
    ###############################################################################
    lopt_float = [ 's_up','s_h','dc','tau' ]
    
    for option in lopt_float:
        try:
            model.__dict__[ option ] = float( model.__dict__[ option ] )
        except:
            DEBUG("option %s is not float or was not provided, but it might irrelevant for your run.")
            pass

    ###############################################################################
    # OPTIONS THAT ARE FLOAT WITH ZERO TO BE IGNORED. SET TO ZERO IF NOT PROVIDED
    ###############################################################################

    lzero = ['offset','interseismic']
    for option in lzero:
        try:
            model.__dict__[ option ] = float( model.__dict__[ option ] )
        except:
            DEBUG("option %s is not float or was not provided. Set to zero")
            model.__dict__[option] = 0.

    ###############################################################################
    # OPTIONS THAT ARE FLOAT WITH 1 HAVING NO EFFECT. SET TO 1 IF NOT PROVIDED
    ###############################################################################

    lone = ['s_up']
    for option in lone:
        try:
            model.__dict__[ option ] = float( model.__dict__[ option ] )
        except:
            DEBUG("option %s is not float or was not provided. Set to 1")
            model.__dict__[option] = 1.


    ###############################################################################
    # ALGORITHM
    ###############################################################################
    # CONVERT TO INT
    lopt_int = [ 'build' ]
    
    for option in lopt_int:
        model.__dict__[ option ] = int( model.__dict__[ option ] )

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
                          # 'dc','sigma','cm_type','cm_norm','tau','m0',\
                          # RAKE
                          'rake_type','rake_value','rake_constraint',\
                          # DATES
                          'dates','rounding',\
                          #ALGO
                          'nnls','build',\
                          # GENERAL
                          'verbose','debug','no_opt','save','plot','tar','print_result','mpck'
                           ]

    for option in lmandatory_options:
        if option not in lmandatory_options:
            ERROR("option not understood or missing: %s" % option ,exit=True )
            mandatory_missing = True
    
    if mandatory_missing:
        ERROR("One mandatory option missing",exit=True)

    return model