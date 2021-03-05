"""
Python wrapper for Okada dc3d code to calculate displacements, tilt, strain and stress for rectangular dislocations.
ref: Okada (1992) [Bull. Seism. Soc. Am., 82, 1018-1040]
https://www.bosai.go.jp/e/dc3d.html
"""


def okada_rde( geometry , array_gps, coor_type='geo' , tensile=False, disp=True, strain=True, stress= True, lam= 0.28758E+11, mu=0.29353E+11, verbose=False ):
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
    l=0,1 : rake_00 & rake_90; optionally provides tensile component for l=2 if requested

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
    from pyeq.green.okada_rde.dc3d import dc3d
    import pyeq.coulomb

    ###############################################################################
    # INIT
    ###############################################################################

    n_dislocations = geometry.shape[0]
    n_gps = array_gps.shape[0]

    if tensile:
        lr = [0,1,2]
        n_slip_component = 3
    else:
        lr = [0,1]
        n_slip_component = 2

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


    # medium constant
    ALPHA = (lam + mu) / (lam + 2.*mu)

    ###############################################################################
    # DISP STRAIN TILT GREEN TENSOR CALCULATION
    ###############################################################################

    GREEN = np.zeros((n_dislocations, n_gps, 3, n_slip_component))
    GREEN_STRAIN = np.zeros((n_dislocations, n_gps, 6, n_slip_component))
    GREEN_STRESS = np.zeros((n_dislocations, n_gps, 6, n_slip_component))
    GREEN_TILT = np.zeros((n_dislocations, n_gps, 2, n_slip_component))


    MESSAGE('calculating displacement, strain, and tilt Green tensor using Okada')

    bar = Bar('', max=n_dislocations*n_gps, suffix='%(percent).1f%% - %(elapsed)ds')

    # we take the fault origin as the coordinates origin
    AL1 = 0.
    AW2 = 0.

    # vectors of constants allowing Okada fault coordinates system to ENU
    np_st = np.radians(array_fault[:,6])
    np_csst = np.cos(np_st)
    np_ssst = np.sin(np_st)
    np_cs2st = np.cos(2. * np_st)
    np_ss2st = np.sin(2. * np_st)

    # DIS for SS, DS, TENSILE
    DIS = {}
    DIS[0] = [1.,0.,0.]
    DIS[1] = [0.,1.,0.]
    DIS[2] = [0.,0.,1.]

    for index_fault in np.arange(n_dislocations):
        AL2 = array_fault[index_fault, 4]
        AW1 = -array_fault[index_fault, 5]
        DEPTH = array_fault[index_fault, 3]
        DIP = array_fault[index_fault, 7]


        # converts the gps coordinates into the fault coordinates
        #array_gps_in_fault_coor = np.copy( array_gps )
        #array_gps_in_fault_coor[:,0] = array_gps_in_fault_coor[:,0] - array_fault[:,1]
        #array_gps_in_fault_coor[:,1] = array_gps_in_fault_coor[:,1] - array_fault[:,2]

        # coefficient for Okada to ENU coordinate conversion
        st = np_st[index_fault]
        csst = np_csst[index_fault]
        ssst = np_ssst[index_fault]
        cs2st = np_cs2st[index_fault]
        ss2st = np_ss2st[index_fault]


        for index_obs in np.arange(n_gps):

            # observation pyacs to Aki
            Y = array_gps[index_obs,0] - array_fault[index_fault, 1]
            X = array_gps[index_obs,1] - array_fault[index_fault, 2]
            Z = 0.

            # observation transform from Aki's to Okada's coordinate system

            Xfault = X * csst + Y * ssst
            Yfault = X * ssst - Y * csst

            for r in lr:
                [DISL1, DISL2, DISL3] = DIS[r]
                [UX, UY, UZ, UXX, UYX, UZX, UXY, UYY, UZY, UXZ, UYZ, UZZ] = \
                dc3d(ALPHA, Xfault, Yfault, Z, DEPTH, DIP, AL1, AL2, AW1, AW2, DISL1, DISL2, DISL3)

                # converts disp, tilt, strain for fault coordinates to Aki coordinates
                disp_X_aki = UX * csst + UY * ssst
                disp_Y_aki = UX * ssst - UY * csst
                disp_Z_aki = -UZ

                tilt_X_aki = -(UXZ * csst + UYZ * ssst)
                tilt_Y_aki = -(UXZ * ssst - UYZ * csst)

                strain_X_aki = UXX * csst * csst + UYY * ssst * ssst + 0.5 * (UXY + UYX) * ss2st
                strain_Y_aki = UXX * ssst * ssst + UYY * csst * csst - 0.5 * (UXY + UYX) * ss2st
                strain_Z_aki = UZZ
                strain_XY_aki =  0.5 * (UXX - UYY) * ss2st - 0.5 * (UXY + UYX) * cs2st
                strain_XZ_aki = -0.5 * (UZX + UXZ) * ssst  + 0.5 * (UYZ + UZY) * csst
                strain_YZ_aki = -0.5 * (UZX + UXZ) * csst  - 0.5 * (UYZ + UZY) * ssst

                # Aki to ENU
                disp_X = disp_Y_aki
                disp_Y = disp_X_aki
                disp_Z = -disp_Z_aki

                tilt_X = tilt_Y_aki
                tilt_Y = tilt_X_aki

                strain_X = strain_Y_aki
                strain_Y = strain_X_aki
                strain_Z = strain_Z_aki
                strain_XY = strain_XY_aki
                strain_XZ = strain_YZ_aki
                strain_YZ = strain_XZ_aki

                GREEN[index_fault, index_obs, :, r] = np.array([disp_X,disp_Y,disp_Z])
                GREEN_STRAIN[index_fault, index_obs, :, r] = np.array([strain_X,strain_Y,strain_Z,strain_XY,strain_XZ,strain_YZ])
                GREEN_TILT[index_fault, index_obs, :, r] = np.array([tilt_X,tilt_Y])
                GREEN_STRESS[index_fault, index_obs, :, r] = pyeq.coulomb.strain2stress( GREEN_STRAIN[index_fault, index_obs, :, r],lam=lam,mu=mu )

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

