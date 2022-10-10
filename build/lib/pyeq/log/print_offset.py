def print_offset( model ):
    """
    print offsets.
    """
    
    ###########################################################################
    # IMPORT
    ###########################################################################
    
    import pyacs.lib.utils
    import numpy as np
    
    
    ###########################################################################
    # RESHAPE
    ###########################################################################
    
    OFFSETS = np.zeros( (  model.np_gps_site.shape[0] , 3  ) )
    
    if model.nconstant > 0:
        TMP = model.estimated_offsets.reshape( model.np_gps_site.shape[0] , -1 ) + model.shift_constant
        OFFSETS[ : , :TMP.shape[1] ] = TMP

    ###########################################################################
    # SAVE OFFSET FILE
    ###########################################################################
    
    fmt = " %10.2f %10.2f %10.2f %s"
    if model.nconstant > 0:
        comment=" Offsets in mm East North Up"
    else:
        comment=" Offsets not estimated."
        
    pyacs.lib.utils.save_np_array_with_string( OFFSETS , model.np_gps_site , fmt,  model.odir+'/info/ts_offset.dat', comment)
    
    ###########################################################################
    # SAVE MODEL.ESTIMATED_OFFSETS_PER_SITE
    ###########################################################################

    model.estimated_offset_per_site = {}
    
    for i in np.arange( model.np_gps_site.shape[0] ):
        model.estimated_offset_per_site[ model.np_gps_site[i] ] = OFFSETS[i,:]
    
    return model