def default( model ):
    """
    Fills pyeq default values when None is not a suitable default value
    """

    # TAR
    model.tar = True
    
    # PLOT
    model.plot = True
    
    # PCK
    model.pck = None
    
    # GEOMETRY
    [ model.geometry_remove_idx , model.geometry_range_lon , model.geometry_range_lat , model.geometry_range_depth ] = [None,None,None,None]
    
    # NO TEMPORAL SMOOTHING
    model.tau = 0.
    
    # NO GEOMETRY CRITICAL SIZE RENORMALIZATION
    model.cm_norm = 1.
    
    # SPATIAL CORRELATION KERNEL
    model.cm_type = 'exponential'
    
    # PRIOR
    model.m0 = 0.
    
    # RAKE_CONSTRAINT : FIXED
    model.rake_constraint = 0.
    
    # ROUNDING FOR GPS DATA - DEFAULT IS DAY
    model.rounding = 'day'
    
    # NO UP BY DEFAULT
    model.up = False
    
    # S_UP - DEFAULT 5.
    model.s_up = 5.
    
    # S_H - DEFAULT - NO RESCALING, USE ORIGINAL VALUE FROM TIME SERIES
    model.s_h = 0.
    
    # EXCLUDE_GPS - USE ALL GPS
    model.lexclude_gps = []
    
    # NNLS ALGORITHM
    model.nnls = 'nnlsm_block_pivot'
    
    # BUILD ALGORITHM
    model.build = 5
    
    # SAVE
    model.save = False
    
    # NO_OPT
    model.no_opt = False
    
    # VERBOSE
    model.verbose = False
    
    # DEBUG
    model.debug = False
    
    return model 
    
    