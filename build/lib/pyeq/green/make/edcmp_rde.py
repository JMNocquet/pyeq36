"""
Python wrapper of edcmp code to calculate displacements, tilt, strain and stress for rectuganlar dislocations.
edcmp is available at: https://www.gfz-potsdam.de/en/section/physics-of-earthquakes-and-volcanoes/infrastructure/tool-development-lab
and must be installed.
reference: Wang, R., F. Lorenzo-Martín and F. Roth (2003), Computation of deformation induced by earthquakes in a multi-layered elastic crust - FORTRAN programs EDGRN/EDCMP, Computer and Geosciences, 29(2), 195-207.
"""


def edcmp_rde( geometry ,
               array_gps,
               coor_type='geo' ,
               disp=True,
               strain=True,
               stress= True,
               lam= 0.28758E+11,
               mu=0.29353E+11,
               verbose=False ):
    """
    :param geometry: a geometry numpy array of shape (nfaults,22).
    :param array_gps: a coordinates x,y or lon, lat of GPS sites (ngps, 2)
    :param coor_type: 'geo' or 'xy'
    :param disp: boolean, True will return displacement Green's matrix at array_gps coordinates
    :param strain: boolean, True will return strain Green's matrix at array_gps coordinates
    :param stress: boolean, True will return stress Green's matrix at dislocation centroids coordinates
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
    from pyeq.green.edcmp_rde import disp_tilt_strain_stress_from_edcmp
    from pyeq.green.edcmp_rde import check_edcmp_executable

    ###############################################################################
    # CHECK EDCMP IS INSTALLED AND PROPERLY WORKING
    ###############################################################################

    if not check_edcmp_executable():
        ERROR('edcmp not installed or not in your PATH. edcmp is available from https://www.gfz-potsdam.de/en/section/physics-of-earthquakes-and-volcanoes/infrastructure/tool-development-lab',exit=True)

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
        array_gps = np.array([x,y,np.zeros(x.shape[0])]).T.reshape(-1,3)
    else:
        array_gps = np.array([array_gps[:,0],array_gps[:,1],np.zeros(array_gps.shape[0])]).T


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
    array_fault[:,4] = geometry[:,3]   # length
    array_fault[:,5] = geometry[:,4]   # width
    array_fault[:,6] = geometry[:,7] # strike
    array_fault[:,7] = geometry[:,8] # dip

    # converts centroid coordinates
    if stress:
        if verbose: MESSAGE('converting RDE centroid coordinates geo to flat Earth coordinates')
        cx,cy = pyacs.lib.coordinates.geo2flat_earth( geometry[:,9], geometry[:,10] )
        cz = -np.sqrt(geometry[:,11]**2)
        array_centroid = np.array([cx,cy,cz]).T


    ###############################################################################
    # DISP STRAIN TILT GREEN TENSOR CALCULATION
    ###############################################################################

    GREEN = np.zeros((n_dislocations, n_gps, 3, 2))
    GREEN_STRAIN = np.zeros((n_dislocations, n_gps, 6, 2))
    GREEN_STRESS = np.zeros((n_dislocations, n_gps, 6, 2))
    GREEN_TILT = np.zeros((n_dislocations, n_gps, 2, 2))


    MESSAGE('calculating displacement, strain, and tilt Green tensor using edcmp')

    bar = Bar('', max=n_dislocations, suffix='%(percent).1f%% - %(elapsed)ds')

    for index in np.arange(n_dislocations):
        single_array_fault = array_fault[index,:].reshape(1,9)
        single_array_fault[0,-1] = 0.

        U_rake_00, STR_rake_00, STS_rake_00, TILT_rake_00 = disp_tilt_strain_stress_from_edcmp(single_array_fault, array_gps, verbose=verbose )
        single_array_fault[0,-1] = 90.
        U_rake_90, STR_rake_90, STS_rake_90, TILT_rake_90 = disp_tilt_strain_stress_from_edcmp(single_array_fault, array_gps, verbose=verbose )

        GREEN[index, :, :, 0] = U_rake_00[:,2:5]
        GREEN[index, :, :, 1] = U_rake_90[:,2:5]

        GREEN_STRAIN[index, :, :, 0] = STR_rake_00[:,2:8]
        GREEN_STRAIN[index, :, :, 1] = STR_rake_90[:,2:8]

        GREEN_STRESS[index, :, :, 0] = STS_rake_00[:,2:8]
        GREEN_STRESS[index, :, :, 1] = STS_rake_90[:,2:8]

        GREEN_TILT[index, :, :, 0] = TILT_rake_00[:,2:4]
        GREEN_TILT[index, :, :, 1] = TILT_rake_90[:,2:4]


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

