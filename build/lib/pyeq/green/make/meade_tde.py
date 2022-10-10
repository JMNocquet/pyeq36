"""
Python wrapper of Meade Triangular Dislocation Element (TDE) code
"""

def meade_tde( geometry , array_gps, coor_type='geo' , disp=True, strain=True, stress= True, poisson_ratio=0.25, lam= 0.28758E+11, mu=0.29353E+11, verbose=False ):
    """
    :param geometry: a geometry numpy array of shape (nfaults,22).
    :param array_gps: a coordinates x,y or lon, lat of GPS sites (ngps, 2)
    :param coor_type: 'geo' or 'xy'
    :param disp: boolean, True will return displacement Green's matrix at array_gps coordinates
    :param strain: boolean, True will return strain Green's matrix at array_gps coordinates
    :param stress: boolean, True will return stress Green's matrix at dislocation centroids coordinates
    :param poisson_ratio: Poisson coefficient, default 0.25
    :param verbose: boolean, verbose mode

    :return: the Green tensor of dim = 4

    Displacement Green tensor is:
    G[i,j,k,l] gives the displacement prediction for dislocation i at site j component k for rake l
    k=0,1,2 = east, north, up
    l=0,1 : rake_00 & rake_90

    Strain Green tensor is:
    G[i,j,k,l] gives the displacement prediction for dislocation i at site j component k for rake l
    k=0,1,2,3,4,5 = east-east, north-north, up-up, east-north, east-up, north-up
    l=0,1 : rake_00 & rake_90

    Unlike Displacement and Strain, Stress Green tensor is calculated at the centroid of faults:
    G[i,j,k,l] gives the displacement prediction for dislocation i at the center of dislocation j component k for rake l
    k=0,1,2,3,4,5 = east-east, north-north, up-up, east-north, east-up, north-up
    l=0,1 : rake_00 & rake_90

    Note: Meade's convention seems to be EN(-U), that is depth is positive, oppositely to Nikkhoo.

    """

    ###############################################################################
    # IMPORT
    ###############################################################################

    import numpy as np
    import pyacs.lib.coordinates
    from progress.bar import Bar

    import pyeq.message.message as MESSAGE
    import pyeq.message.warning as WARNING
    import pyeq.message.error as ERROR
    from pyeq.green.meade_tde import tde

    ###############################################################################
    # INIT
    ###############################################################################

    n_dislocations = geometry.shape[0]
    n_gps = array_gps.shape[0]
    slip = 1.0

    ###############################################################################
    # CONVERT TO FLAT EARTH APPROXIMATION IF NEEDED
    ###############################################################################

    if coor_type == 'geo':
        if verbose: MESSAGE('converting GPS coordinates geo to flat Earth coordinates')
        # convert GPS coordinates
        [lon, lat] = [array_gps[:,0], array_gps[:,1]]
        (x,y) = pyacs.lib.coordinates.geo2flat_earth( lon, lat )
        x = np.atleast_1d( x )
        y = np.atleast_1d( y )
        array_gps = np.array([x,y]).T

    # converts TDE coordinates
    if verbose: MESSAGE('converting TDE vertices coordinates geo to flat Earth coordinates')

    tdis_x1,tdis_y1 = pyacs.lib.coordinates.geo2flat_earth( geometry[:,12], geometry[:,13] )
    tdis_z1 = np.sqrt( geometry[:,14]**2 )

    tdis_x2,tdis_y2 = pyacs.lib.coordinates.geo2flat_earth( geometry[:,15], geometry[:,16] )
    tdis_z2 = np.sqrt( geometry[:,17]**2 )

    tdis_x3,tdis_y3 = pyacs.lib.coordinates.geo2flat_earth( geometry[:,18], geometry[:,19] )
    tdis_z3 = np.sqrt( geometry[:,20]**2 )

    # so far, everything is in km, convert to m
    X = np.array( [tdis_x1,tdis_x2,tdis_x3] ).T
    Y = np.array( [tdis_y1,tdis_y2,tdis_y3] ).T
    Z = np.array( [tdis_z1,tdis_z2,tdis_z3] ).T

    # check dim
    if X.ndim ==1:
        X = np.array(X , ndmin=2)
        Y = np.array(Y , ndmin=2)
        Z = np.array(Z , ndmin=2)


    # converts centroid coordinates
    if stress:
        if verbose: MESSAGE('converting TDE centroid coordinates geo to flat Earth coordinates')
        cx,cy = pyacs.lib.coordinates.geo2flat_earth( geometry[:,9], geometry[:,10] )
        cz = geometry[:,11]


    ###############################################################################
    # DISP GREEN TENSOR CALCULATION
    ###############################################################################


    if disp:

        GREEN = np.zeros((n_dislocations, n_gps, 3, 2))

        # observation points
        SX = array_gps[:, 0]
        SY = array_gps[:, 1]
        SZ = array_gps[:, 0] * 0.0

        MESSAGE('calculating displacement TDE Green tensor using Meade')

        bar = Bar('', max=n_dislocations, suffix='%(percent).1f%% - %(elapsed)ds')

        for index in np.arange(n_dislocations):

            U_rake_00 = tde.calc_tri_displacements(SX, SY, SZ, X[index], Y[index], Z[index], poisson_ratio, -1.0, 0.0, 0.0)
            U_rake_90 = tde.calc_tri_displacements(SX, SY, SZ, X[index], Y[index], Z[index], poisson_ratio, 0.0, 0.0, -1.0)

            # Meade tde convention: positive Uz downward - corrected 18/02/2020 by adding -

            GREEN[index, :, :, 0] = np.array( [ U_rake_00['x'],U_rake_00['y'],-U_rake_00['z'] ]).T
            GREEN[index, :, :, 1] = np.array( [ U_rake_90['x'],U_rake_90['y'],-U_rake_90['z'] ]).T

            bar.next()

        bar.finish()

    ###############################################################################
    # STRAIN GREEN TENSOR CALCULATION
    ###############################################################################


    if strain:

        GREEN_STRAIN = np.zeros((n_dislocations, n_gps, 6, 2))

        # observation points
        SX = array_gps[:, 0] * 1.E3
        SY = array_gps[:, 1] * 1.E3
        SZ = array_gps[:, 0] * 0.0

        MESSAGE('calculating strain Green tensor using Meade')

        bar = Bar('', max=n_dislocations, suffix='%(percent).1f%% - %(elapsed)ds')

        for index in np.arange(n_dislocations):

            S_rake_00 = tde.calc_tri_strains(SX, SY, SZ, X[index], Y[index], Z[index], poisson_ratio, -1.0, 0.0, 0.0)
            S_rake_90 = tde.calc_tri_strains(SX, SY, SZ, X[index], Y[index], Z[index], poisson_ratio, 0.0, 0.0, -1.0)

            GREEN_STRAIN[index, :, :, 0] = np.array( [ S_rake_00['xx'],S_rake_00['yy'], S_rake_00['zz'],S_rake_00['xy'], S_rake_00['xz'], S_rake_00['yz'] ]).T
            GREEN_STRAIN[index, :, :, 1] = np.array( [ S_rake_90['xx'],S_rake_90['yy'], S_rake_90['zz'],S_rake_90['xy'], S_rake_90['xz'], S_rake_90['yz'] ]).T

            bar.next()

        bar.finish()

    ###############################################################################
    # STRESS GREEN TENSOR CALCULATION (AT SUBFAULTS CENTROIDS)
    ###############################################################################


    if stress:

        GREEN_STRESS = np.zeros((n_dislocations, n_dislocations, 6, 2))

        # observation points
        SX = cx * 1.E3
        SY = cy * 1.E3
        SZ = cz * 1.E3

        MESSAGE('calculating stress Green tensor using Meade')

        bar = Bar('', max=n_dislocations, suffix='%(percent).1f%% - %(elapsed)ds')


        for index in np.arange(n_dislocations):


            S_rake_00 = tde.strain_to_stress( tde.calc_tri_strains(SX, SY, SZ, X[index], Y[index], Z[index], poisson_ratio, -1.0, 0.0, 0.0) , lam, mu )
            S_rake_90 = tde.strain_to_stress( tde.calc_tri_strains(SX, SY, SZ, X[index], Y[index], Z[index], poisson_ratio, 0.0, 0.0, -1.0) , lam, mu )

            GREEN_STRESS[index, :, :, 0] = np.array( [ S_rake_00['xx'],S_rake_00['yy'], S_rake_00['zz'],S_rake_00['xy'], S_rake_00['xz'], S_rake_00['yz'] ]).T
            GREEN_STRESS[index, :, :, 1] = np.array( [ S_rake_90['xx'],S_rake_90['yy'], S_rake_90['zz'],S_rake_90['xy'], S_rake_90['xz'], S_rake_90['yz'] ]).T

            bar.next()

        bar.finish()

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

