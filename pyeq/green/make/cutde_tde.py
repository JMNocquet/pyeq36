"""
Python wrapper of CUTDE Half space Triangular Dislocation Element (TDE) code
source: https://github.com/tbenthompson/cutde
"""


def cutde_tde( geometry ,
                 array_gps,
                 coor_type='geo' ,
                 tensile=False,
                 disp=True,
                 strain=False,
                 stress= False,
                 poisson_ratio=0.25,
                 lam= 0.28758E+11,
                 mu=0.29353E+11,
                 verbose=False ):
    """

    Compute the green function for a triangular dislocation.

    :param geometry: a geometry numpy array of shape (nfaults,22).
    :param array_gps: a coordinates x,y or lon, lat of observation sites. Can be of shape (n, 2) or (n,3) for z<0 observations
    :param coor_type: 'geo' or 'xy'
    :param tensile: boolean, whether tensile component is added to the Green tensor
    :param disp: boolean, True will return displacement Green's matrix at array_gps coordinates
    :param strain: boolean, True will return strain Green's matrix at array_gps coordinates
    :param stress: boolean, True will return stress Green's matrix at dislocation centroids coordinates
    :param poisson_ratio: Poisson coefficient, default 0.25
    :param lam,mu: Lamé parameters used for strain-stress conversion only
    :param verbose: boolean, verbose mode

    :return: the Green tensor of dim = 4

    Displacement Green tensor is:
    G[i,j,k,l] gives the displacement prediction for dislocation i at site j component k for rake l
    k=0,1,2 = east, north, up
    l=0,1 : rake_00 & rake_90; ; if tensile is True, tensile component is added for l=2

    Strain Green tensor is:
    G[i,j,k,l] gives the displacement prediction for dislocation i at site j component k for rake l
    k=0,1,2,3,4,5 = east-east, north-north, up-up, east-north, east-up, north-up
    l=0,1 : rake_00 & rake_90; if tensile is True, tensile component is added for l=2

    Unlike Displacement and Strain, Stress Green tensor is calculated at the centroid of faults:
    G[i,j,k,l] gives the displacement prediction for dislocation i at the center of dislocation j component k for rake l
    k=0,1,2,3,4,5 = east-east, north-north, up-up, east-north, east-up, north-up
    l=0,1 : rake_00 & rake_90: ; if tensile is True, tensile component is added for l=2

    Note: apparently, conventions are: ENU (positive Upward) ss>0 left-lateral, ds>0 inverse if triangle vertices are ordered upward.
    """
    #TODO check lam, mu and poisson_ratio coherence

    ###############################################################################
    # IMPORT
    ###############################################################################

    import numpy as np
    import pyacs.lib.coordinates
    from progress.bar import Bar

    import pyeq.message.message as MESSAGE
    import pyeq.message.warning as WARNING
    import pyeq.message.error as ERROR
    import pyeq.coulomb
    from cutde.halfspace import disp_matrix, strain_matrix

    ###############################################################################
    # INIT
    ###############################################################################

    n_dislocations = geometry.shape[0]
    n_gps = array_gps.shape[0]

    ###############################################################################
    # CONVERT TO FLAT EARTH APPROXIMATION IF NEEDED
    ###############################################################################

    if array_gps.shape[1] == 3:
        z = array_gps[:,2]
    else:
        z = np.zeros(n_gps)

    if coor_type == 'geo':
        if verbose: MESSAGE('converting GPS coordinates geo to flat Earth coordinates')
        # convert GPS coordinates
        [lon, lat] = [array_gps[:,0], array_gps[:,1]]
        (x,y) = pyacs.lib.coordinates.geo2flat_earth( lon, lat )
        x = np.atleast_1d( x )
        y = np.atleast_1d( y )
        array_gps = np.array([x,y,z]).T.reshape(-1,3)
    else:
        array_gps = np.array([array_gps[:,0],array_gps[:,1],z]).T

    # converts TDE coordinates
    if verbose: MESSAGE('converting TDE vertices coordinates geo to flat Earth coordinates')

    tdis_x1,tdis_y1 = pyacs.lib.coordinates.geo2flat_earth( geometry[:,12], geometry[:,13] )
    tdis_z1 = -np.sqrt( geometry[:,14]**2 )

    tdis_x2,tdis_y2 = pyacs.lib.coordinates.geo2flat_earth( geometry[:,15], geometry[:,16] )
    tdis_z2 = -np.sqrt( geometry[:,17]**2 )

    tdis_x3,tdis_y3 = pyacs.lib.coordinates.geo2flat_earth( geometry[:,18], geometry[:,19] )
    tdis_z3 = -np.sqrt( geometry[:,20]**2 )

    # case there is only one dislocation
    tdis_x1 = np.atleast_1d(tdis_x1)
    tdis_y1 = np.atleast_1d(tdis_y1)
    tdis_x2 = np.atleast_1d(tdis_x2)
    tdis_y2 = np.atleast_1d(tdis_y2)
    tdis_x3 = np.atleast_1d(tdis_x3)
    tdis_y3 = np.atleast_1d(tdis_y3)

    # so far, everything is in km, convert to m and put in proper format
    TRI = np.array([tdis_x1,tdis_y1,tdis_z1,tdis_x2,tdis_y2,tdis_z2,tdis_x3,tdis_y3,tdis_z3]).T.reshape(-1,3,3)

    # converts centroid coordinates
    if stress:
        if verbose: MESSAGE('converting TDE centroid coordinates geo to flat Earth coordinates')
        cx,cy = pyacs.lib.coordinates.geo2flat_earth( geometry[:,9], geometry[:,10] )
        cz = -np.sqrt(geometry[:,11]**2)
        C = np.array([cx,cy,cz]).T


    ###############################################################################
    # DISP GREEN TENSOR CALCULATION
    ###############################################################################

    import pyeq.lib.geometry

    if tensile:
        GREEN = np.zeros((n_dislocations, n_gps, 3, 3))
        print(GREEN.shape)
        GREEN_STRAIN = np.zeros((n_dislocations, n_gps, 6, 3))
        GREEN_STRESS = np.zeros((n_dislocations, n_gps, 6, 3))
    else:
        GREEN = np.zeros((n_dislocations, n_gps, 3, 2))
        GREEN_STRAIN = np.zeros((n_dislocations, n_gps, 6, 2))
        GREEN_STRESS = np.zeros((n_dislocations, n_gps, 6, 2))

    MESSAGE('calculating displacement Green TDE tensor using B. Thompson cutde')
    green_cutde_disp = disp_matrix(array_gps, TRI, poisson_ratio)

    # cuttde convention is (N_OBS_PTS, 3, N_SRC_TRIS, 3)
    # pyacs convention is
    # Displacement Green tensor is:
    # G[i,j,k,l] gives the displacement prediction for dislocation i at site j component k for rake l
    # k=0,1,2 = east, north, up
    # l=0,1 : rake_00 & rake_90; ; if tensile is True, tensile component is added for l=2
    # needs to re-order dimension
    MESSAGE("Swaping axis from cutde to pyacs convention")
    GREEN = np.swapaxes(np.swapaxes(green_cutde_disp, 0, 2),1,2)

    print(GREEN.shape)

    ###############################################################################
    # RETURN
    ###############################################################################

    if disp:
        if strain:
            if stress:
                return GREEN, GREEN_STRAIN, GREEN_STRESS
            else:
                return GREEN , GREEN_STRAIN
        else:
            if stress:
                return GREEN , GREEN_STRESS
            else:
                return GREEN

    else:
        if strain:
            if stress:
                return GREEN_STRAIN, GREEN_STRESS
            else:
                return GREEN_STRAIN
        else:
            return GREEN_STRESS

