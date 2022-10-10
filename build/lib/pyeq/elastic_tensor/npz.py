"""
manipulate some npy format used by pyeq
"""
###############################################################################
def read_pyeq_input_npz(input_npz):
###############################################################################
    """
    reads an input_npz file used by pyeq and created by pyeq_make_green.py
    """

    # import
    
    import numpy as np
    import pyeq.lib.lib_inversion
    from pyeq.lib.errors import PyeqNpzReadError , PyeqNpzFormatError


    try:
        INPUT = np.load(input_npz)
    except:
        raise PyeqNpzReadError( __name__ , input_npz )

    try:
        GEOMETRY    = INPUT["GEOMETRY"]
        Dm          = INPUT["Dm"]
        GREEN       = INPUT["GREEN"]
        GREEN_UP    = INPUT["GREEN_UP"]
        OBS         = INPUT["OBS"]
        NAME_OBS    = INPUT["NAME_OBS"]
        OBS_UP      = INPUT["OBS_UP"]
        NAME_OBS_UP = INPUT["NAME_OBS_UP"]
    except:
        raise PyeqNpzFormatError( __name__ , input_npz )
 
    names=['rdis_long','rdis_lat','rdis_depth','rdis_length','rdis_width',\
           'rdis_area','ratio_rdis_tdis','strike','dip',\
           'centroid_long','centroid_lat','centroid_depth',\
           'tdis_long1','tdis_lat1','tdis_depth1',\
           'tdis_long2','tdis_lat2','tdis_depth2',\
           'tdis_long3','tdis_lat3','tdis_depth3','tdis_area']
    
    SGEOMETRY=pyeq.lib.lib_inversion.numpy_array_2_numpy_recarray(GEOMETRY,names)
 
        
    
    return(SGEOMETRY , GEOMETRY, Dm, GREEN, GREEN_UP, OBS, NAME_OBS, OBS_UP, NAME_OBS_UP)