
def npy_geometry_to_array_and_recarray( npy_geometry , verbose=False ):
    """
    Converts a npy geometry to array and recarray

    :param npy_geometry: geometry as npy file
    :param verbose: verbose mode

    """
    # import
    
    import sys
    import numpy as np
    from pyeq.lib import lib_inversion

    # start
    if verbose:
        print("-- Loading geometry file ", npy_geometry )

    # reads geometry npy file
    try:
        GEOMETRY=np.load( npy_geometry )
     
    except:
        print("!!! Could not load ", npy_geometry )
        sys.exit()
    
    # converts array to recarray
    names=[\
           'rdis_long','rdis_lat','rdis_depth','rdis_length','rdis_width',\
           'rdis_area','ratio_rdis_tdis','strike','dip',\
           'centroid_long','centroid_lat','centroid_depth',\
           'tdis_long1','tdis_lat1','tdis_depth1',\
           'tdis_long2','tdis_lat2','tdis_depth2',\
           'tdis_long3','tdis_lat3','tdis_depth3','tdis_area']
    
    try:
        SGEOMETRY = lib_inversion.numpy_array_2_numpy_recarray(GEOMETRY,names)
    except:
        print(GEOMETRY.shape)
        print(len(names))
        print( red("[PYEQ ERROR] Error converting GEOMETRY array to recarray. GEOMETRY.shape=(%d,%d)" % (GEOMETRY.shape[0],GEOMETRY.shape[1])) )
        
    # print some informations about the geometry
    
    if verbose:
    
        print(("  -- geometry includes %04d subfaults" % SGEOMETRY.shape[0]))
    
        min_lon=np.min(SGEOMETRY.centroid_long)
        max_lon=np.max(SGEOMETRY.centroid_long)
        
        min_lat=np.min(SGEOMETRY.centroid_lat)
        max_lat=np.max(SGEOMETRY.centroid_lat)
        
        min_depth=np.min(SGEOMETRY.centroid_depth)
        max_depth=np.max(SGEOMETRY.centroid_depth)
        
        print((" -- geometry bounds read from centroids %.5lf/%.5lf/%.5lf/%.5lf and depth %.5lf/%.5lf " % (min_lon,max_lon,min_lat,max_lat,min_depth,max_depth)))
    
    return( GEOMETRY , SGEOMETRY )

def dat_geometry_to_array_and_recarray( dat_geometry , verbose=False ):
    """
    Converts a dat geometry to array and recarray

    :param dat_geometry: geometry as dat file
    :param verbose: verbose mode

    """
    # import
    
    import sys
    import numpy as np
    from pyeq.lib import lib_inversion

    # start
    if verbose:
        print("-- Loading geometry file ", dat_geometry )

    # reads geometry npy file
    try:
        GEOMETRY= np.atleast_2d( np.genfromtxt( dat_geometry ) )[:,1:]
     
    except:
        print("!!! Could not read ", dat_geometry )
        sys.exit()
    
    # converts array to recarray
    names=[\
           'rdis_long','rdis_lat','rdis_depth','rdis_length','rdis_width',\
           'rdis_area','ratio_rdis_tdis','strike','dip',\
           'centroid_long','centroid_lat','centroid_depth',\
           'tdis_long1','tdis_lat1','tdis_depth1',\
           'tdis_long2','tdis_lat2','tdis_depth2',\
           'tdis_long3','tdis_lat3','tdis_depth3','tdis_area']
    
    try:
        SGEOMETRY = lib_inversion.numpy_array_2_numpy_recarray(GEOMETRY,names)
    except:
        print(GEOMETRY.shape)
        print(len(names))
        print( red("[PYEQ ERROR] Error converting GEOMETRY array to recarray. GEOMETRY.shape=(%d,%d)" % (GEOMETRY.shape[0],GEOMETRY.shape[1])) )
        
    # print some informations about the geometry
    
    if verbose:
    
        print(("  -- geometry includes %04d subfaults" % SGEOMETRY.shape[0]))
    
        min_lon=np.min(SGEOMETRY.centroid_long)
        max_lon=np.max(SGEOMETRY.centroid_long)
        
        min_lat=np.min(SGEOMETRY.centroid_lat)
        max_lat=np.max(SGEOMETRY.centroid_lat)
        
        min_depth=np.min(SGEOMETRY.centroid_depth)
        max_depth=np.max(SGEOMETRY.centroid_depth)
        
        print((" -- geometry bounds read from centroids %.5lf/%.5lf/%.5lf/%.5lf and depth %.5lf/%.5lf " % (min_lon,max_lon,min_lat,max_lat,min_depth,max_depth)))
    
    return( GEOMETRY , SGEOMETRY )



