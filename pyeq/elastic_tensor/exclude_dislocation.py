###############################################################################
def exclude_dislocation(G , geometry, Dm, \
                        range_lon=None, \
                        range_lat=None, \
                        range_depth=None, \
                        exclude_idx=None, \
                        verbose=False):
###############################################################################
    """
    
    From a 4-dimension elastic tensor G, make a selection of subfaults using some criterion
    
    The input G tensor is assumed to be organized as follows:
    
    G(i,j,k,l) is the prediction for dislocation j at site i component k for rake l
    k=0,1,2 = east, north, up
    l=0,1 : rake_principal & rake_conjugate

    geometry has the following columns:
    0 : rdis_long
    1 : rdis_lat
    2 : rdis_depth
    3 : rdis_length
    4 : rdis_width
    5 : rdis_area
    6 : ratio_rdis_tdis
    7 : strike
    8 : dip
    9 : centroid_long
    10: centroid_lat
    11: centroid_depth
    12: tdis_long1
    13: tdis_lat1
    14: tdis_depth1
    15: tdis_long2
    16: tdis_lat2
    17: tdis_depth2
    18: tdis_long3
    19: tdis_lat3
    20: tdis_depth3
    21: tdis_area
    
    :param range_lon: removes all subfaults with centroid strictly outside the provided longitude range (unit dec.deg)
    :param range_lat: removes all subfault with centroid strictly outside the provided latitude range (unit dec.deg)
    :param range_depth: removes all subfaults with centroid strictly outside the provided depth range (unit dec.deg)
    :param exclude_idx: removes a list of subfault by index 
     
    """

    ###########################################################################
    # IMPORT 
    ###########################################################################
    
    import numpy as np
    import pyeq.lib.lib_inversion
    import pyeq.message.message as MESSAGE
    import pyeq.message.verbose_message as VERBOSE
    import pyeq.message.error as ERROR

    ###########################################################################
    # CHECK ARGUMENT
    ###########################################################################
    
    # Dimension of input G
    if G.ndim != 4:
        ERROR(("Input Green array MUST be of dimension 4. G.ndim=%d" % G.ndim))

    ###########################################################################
    # IDX CASE
    # MUST BE FIRST, OTHERWISE IDX ARE WRONG
    ###########################################################################

    
    if exclude_idx is not None:
        VERBOSE(("Removing dislocations for index %d " % exclude_idx ))
        np_idx = np.delete( np.arange( G.shape[0] ) , np.array(exclude_idx) )
        G = G[np_idx,:,:,:]
        geometry = geometry[ np_idx, : ]
        Dm = Dm[np_idx].T[np_idx].T  
    
    ###########################################################################
    # MIN/MAX CASE
    ###########################################################################

    if range_lon is not None:
        VERBOSE(("Removing dislocations outside longitude range: %.3lf, %.3lf " % tuple(range_lon) ))
        
        [ min_lon , max_lon ] = range_lon
        np_idx =  np.where( ( geometry[:,9] >= min_lon ) & ( geometry[:,9] <= max_lon ) )[0]
        G = G[np_idx,:,:,:]
        geometry = geometry[ np_idx, : ]
        Dm = Dm[np_idx].T[np_idx].T  
        
    if range_lat is not None:
        VERBOSE(("Removing dislocations outside latitude range: %.3lf, %.3lf " % tuple(range_lat) ))

        [ min_lat , max_lat ] = range_lat
        np_idx =  np.where( ( geometry[:,10] >= min_lat ) & ( geometry[:,10] <= max_lat ) )[0]
        G = G[np_idx,:,:,:]
        geometry = geometry[ np_idx, : ]
        Dm = Dm[np_idx].T[np_idx].T  
        
    if range_depth is not None:
        VERBOSE(("Removing dislocations outside depth range: %.1lf, %.1lf " % tuple(range_depth) ))

        [ min_depth , max_depth ] = range_depth
        np_idx =  np.where( ( np.fabs( geometry[:,11] ) >= min_depth ) & ( np.fabs( geometry[:,11] ) <= max_depth ) )[0]
        G = G[np_idx]
        geometry = geometry[ np_idx ]
        Dm = Dm[np_idx].T[np_idx].T  

    names=['rdis_long','rdis_lat','rdis_depth','rdis_length','rdis_width',\
           'rdis_area','ratio_rdis_tdis','strike','dip',\
           'centroid_long','centroid_lat','centroid_depth',\
           'tdis_long1','tdis_lat1','tdis_depth1',\
           'tdis_long2','tdis_lat2','tdis_depth2',\
           'tdis_long3','tdis_lat3','tdis_depth3','tdis_area']
    
    SGEOMETRY = pyeq.lib.lib_inversion.numpy_array_2_numpy_recarray(geometry,names)

    return G , geometry, SGEOMETRY , Dm
        
