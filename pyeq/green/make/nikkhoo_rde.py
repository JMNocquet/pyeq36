"""
Python wrapper of Nikkhoo Triangular Dislocation Element (TDE) code - Port to Python from E. Lindsey using B. Thompson
source: https://github.com/ericlindsey/tdcalc
"""


def nikkhoo_rde( geometry ,
                 array_gps,
                 coor_type='geo' ,
                 tensile=False,
                 disp=True,
                 strain=True,
                 stress= True,
                 poisson_ratio=0.25,
                 lam=0.28758E+11,
                 mu=0.29353E+11,
                 verbose=False ):
    """

    Compute the green function for a rectangular dislocation. rectangle are obtained from two adjacent triangles
    defining a rectangle.
    This routine is mainly for check, Okada's formula allowing a direct calculation in this case.

    :param geometry: a geometry numpy array of shape (nfaults,22).
    :param array_gps: a coordinates x,y or lon, lat of GPS sites (ngps, 2)
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
    from pyeq.green.nikkhoo_tde import tdcalc
    import pyeq.coulomb

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


    # converts RDE coordinates
    if verbose: MESSAGE('converting RDE origin point coordinates geo to flat Earth coordinates')

    rdis_x,rdis_y = pyacs.lib.coordinates.geo2flat_earth( geometry[:,0], geometry[:,1] )
    rdis_z = -np.sqrt( geometry[:,2]**2 )

    # case there is only one dislocation
    rdis_x = np.atleast_1d(rdis_x)
    rdis_y = np.atleast_1d(rdis_y)

    # so far, everything is in km which is fine with disp_tilt_stress_from_edcmp
    # array_faults : 2D numpy array of faults parameters [slip,xf,yf,depth,length,width,strike,dip,rake]
    array_fault = np.zeros((geometry.shape[0],9))

    array_fault[:,0] = 1.  # slip
    array_fault[:,1] = rdis_x
    array_fault[:,2] = rdis_y
    array_fault[:,3] = rdis_z
    array_fault[:,4] = geometry[:,3] # length
    array_fault[:,5] = geometry[:,4] # width
    array_fault[:,6] = geometry[:,7] # strike
    array_fault[:,7] = geometry[:,8] # dip

    # computes delta_x, delta_y, delta_z for the rectangular dislocation vertices
    alpha = np.radians(-array_fault[:,6]) + np.pi / 2.

    delta_x = np.cos(np.radians( array_fault[:,7] )) * np.cos(alpha - np.pi / 2.) * array_fault[:,5]
    delta_y = np.cos(np.radians( array_fault[:,7] )) * np.sin(alpha - np.pi / 2.) * array_fault[:,5]
    delta_z = np.sin(np.radians(array_fault[:,7])) * array_fault[:,5]
    DELTA_BOTTOM = np.array([delta_x, delta_y, -delta_z]).T

    delta_x =  np.cos(alpha) * array_fault[:,4]
    delta_y =  np.sin(alpha) * array_fault[:,4]
    delta_z = delta_x * 0.0
    DELTA_TOP = np.array([delta_x, delta_y, delta_z]).T


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
        GREEN_STRAIN = np.zeros((n_dislocations, n_gps, 6, 3))
        GREEN_STRESS = np.zeros((n_dislocations, n_gps, 6, 3))
    else:
        GREEN = np.zeros((n_dislocations, n_gps, 3, 2))
        GREEN_STRAIN = np.zeros((n_dislocations, n_gps, 6, 2))
        GREEN_STRESS = np.zeros((n_dislocations, n_gps, 6, 2))

    if verbose:
        MESSAGE('calculating displacement Green tensor using Nikkhoo for rectangular dislocations')

    slip_00 = [1., 0., 0]
    slip_90 = [0., 1., 0.]
    slip_tensile = [0., 0., 1.]

    bar = Bar('', max=n_dislocations, suffix='%(percent).1f%% - %(elapsed)ds')
    for index in np.arange(n_dislocations):

        # get rectangle coordinates
        A = array_fault[index,1:4]
        B = A + DELTA_TOP[index]
        C = B + DELTA_BOTTOM[index]
        D = A + DELTA_BOTTOM[index]

        TRI1 = np.array([ A,B,C ])
        if not pyeq.lib.geometry.is_triangle_up(A,B,C):
            TRI1 = np.array([ B,A,C ])

        TRI2 = np.array([ A,C,D ])
        if not pyeq.lib.geometry.is_triangle_up(A,C,D):
            TRI2 = np.array([ C,A,D ])


        U_rake_00 = tdcalc.TDdispHS(array_gps,TRI1,slip_00,poisson_ratio)
        U_rake_90 = tdcalc.TDdispHS(array_gps,TRI1,slip_90,poisson_ratio)

        U_rake_00 += tdcalc.TDdispHS(array_gps,TRI2,slip_00,poisson_ratio)
        U_rake_90 += tdcalc.TDdispHS(array_gps,TRI2,slip_90,poisson_ratio)

        if tensile:
            U_tensile = tdcalc.TDdispHS(array_gps,TRI1,slip_tensile,poisson_ratio)
            U_tensile += tdcalc.TDdispHS(array_gps,TRI2,slip_tensile,poisson_ratio)
            GREEN[index, :, :, 2] = U_tensile


        GREEN[index, :, :, 0] = U_rake_00
        GREEN[index, :, :, 1] = U_rake_90

        # sum green disp for TRI1 & TRI2
        GREEN[index, :, :, 0] = tdcalc.TDdispHS(array_gps,TRI1,slip_00,poisson_ratio)
        GREEN[index, :, :, 1] = tdcalc.TDdispHS(array_gps,TRI1,slip_90,poisson_ratio)
        GREEN[index, :, :, 0] += tdcalc.TDdispHS(array_gps,TRI2,slip_00,poisson_ratio)
        GREEN[index, :, :, 1] += tdcalc.TDdispHS(array_gps,TRI2,slip_90,poisson_ratio)

        # sum green strain for TRI1 & TRI2
        GREEN_STRAIN[index, :, :, 0] = tdcalc.TDstrainHS(array_gps,TRI1,slip_00,poisson_ratio)
        GREEN_STRAIN[index, :, :, 1] = tdcalc.TDstrainHS(array_gps,TRI1,slip_90,poisson_ratio)
        GREEN_STRAIN[index, :, :, 0] += tdcalc.TDstrainHS(array_gps,TRI2,slip_00,poisson_ratio)
        GREEN_STRAIN[index, :, :, 1] += tdcalc.TDstrainHS(array_gps,TRI2,slip_90,poisson_ratio)

        # strain2stress
        GREEN_STRESS[index, :, :, 0] = pyeq.coulomb.strain2stress( GREEN_STRAIN[index, :, :, 0],lam=lam,mu=mu )
        GREEN_STRESS[index, :, :, 1] = pyeq.coulomb.strain2stress( GREEN_STRAIN[index, :, :, 1],lam=lam,mu=mu )


        if tensile:
            GREEN[index, :, :, 2] = tdcalc.TDdispHS(array_gps,TRI1,slip_tensile,poisson_ratio)
            GREEN[index, :, :, 2] += tdcalc.TDdispHS(array_gps,TRI2,slip_tensile,poisson_ratio)

            GREEN_STRAIN[index, :, :, 2] = tdcalc.TDstrainHS(array_gps,TRI1,slip_tensile,poisson_ratio)
            GREEN_STRAIN[index, :, :, 2] += tdcalc.TDstrainHS(array_gps,TRI2,slip_tensile,poisson_ratio)

            GREEN_STRESS[index, :, :, 2] = pyeq.coulomb.strain2stress( GREEN_STRAIN[index, :, :, 2],lam=lam,mu=mu )


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

