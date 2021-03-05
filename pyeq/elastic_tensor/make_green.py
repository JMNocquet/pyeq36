def make_green( geometry, sgeometry , array_gps, type='tde', verbose=False):
    """
    Creates a Green tensor

    GREEN IS A TENSOR OF DIM 4
    GREEN(i,j,k,l) is the prediction for dislocation i at site j component k for rake l
    k=0,1,2 = east, north, up
    l=0,1 : rake_00 & rake_90

    """

    ###############################################################################
    # IMPORT
    ###############################################################################

    import numpy as np
    import pyacs.lib.coordinates
    from pyeq.lib import eq_disloc_3d as Dislocation

    # TDE
    if type == 'tde':
        from pyeq.lib.meade_tde import tde
        Poisson_Ratio=0.25
        TDE=True
    else:
        TDE=False


    ###############################################################################
    # H_FAULTS INITALIZATION
    ###############################################################################

    H_fault_lp = {}
    H_fault_xy = {}

    print("-- Building dislocations")

    # TRIANGULAR DISLOCATONS
    if TDE:
        X = np.zeros(3)
        Y = np.zeros(3)
        Z = np.zeros(3)

        XYZ = np.zeros((3, 3))

        H_TDE = {}
        print("-- Triangular dislocations will be used")

    for i in range( sgeometry.shape[0]):

        if verbose:
            print('  -- ', i, ' / ', sgeometry.shape[0])

        [rdis_long, rdis_lat, rdis_depth, rdis_length, rdis_width, rdis_area, ratio_rdis_tdis, strike, dip, \
         centroid_long, centroid_lat, centroid_depth, \
         tdis_long1, tdis_lat1, tdis_depth1, tdis_long2, tdis_lat2, tdis_depth2, tdis_long3, tdis_lat3, tdis_depth3, \
         tdis_area] = np.array(list( sgeometry[i]))

        # triangular dislocation
        if TDE:
            (X[0], Y[0]) = pyacs.lib.coordinates.geo2flat_earth(tdis_long1, tdis_lat1);
            Z[0] = tdis_depth1
            (X[1], Y[1]) = pyacs.lib.coordinates.geo2flat_earth(tdis_long2, tdis_lat2);
            Z[1] = tdis_depth2
            (X[2], Y[2]) = pyacs.lib.coordinates.geo2flat_earth(tdis_long3, tdis_lat3);
            Z[2] = tdis_depth3

            XYZ[:, 0] = X
            XYZ[:, 1] = Y
            XYZ[:, 2] = Z

            H_TDE[i] = np.copy(XYZ)

        # rectangular dislocations
        else:

            depth = np.sqrt(rdis_depth ** 2)
            index_fault = i
            if (rdis_long > 180.): rdis_long = rdis_long - 360.

            lon = rdis_long
            lat = rdis_lat
            length = rdis_length
            width = rdis_width
            area = rdis_area

            # fake values

            rake = 0.0
            max_slip = 0.0

            (x, y) = pyacs.lib.coordinates.geo2flat_earth(lon, lat)
            dislocation_lp = Dislocation.Dislocation(index_fault, lon, lat, depth, strike, dip, length, width, area, rake,
                                                     max_slip)
            dislocation_xy = Dislocation.Dislocation(index_fault, x, y, depth, strike, dip, length, width, area, rake,
                                                     max_slip)
            H_fault_xy[index_fault] = dislocation_xy
            H_fault_lp[index_fault] = dislocation_lp


    ###############################################################################
    # CREATES GREEN FUNCTIONS FOR HORIZONTAL COMPONENTS
    ###############################################################################

    n_dislocations = sgeometry.shape[0]
    n_gps = array_gps.shape[0]
    slip = 1.0

    print("-- Creating Green's functions matrix for horizontal components")

    # GREEN IS A TENSOR OF DIM 4
    # GREEN(i,j,k,l) is the prediction for dislocation i at site j component k for rake l
    # k=0,1,2 = east, north, up
    # l=0,1 : rake_00 & rake_90

    GREEN = np.zeros((n_dislocations, n_gps, 3, 2))

    if TDE:
        # observation points
        SX = array_gps[:, 0] * 1.E3
        SY = array_gps[:, 1] * 1.E3
        SZ = array_gps[:, 0] * 0.0

        for index in range(sgeometry.shape[0]):

            if verbose:
                print('  -- ', index, ' / ', sgeometry.shape[0])

            X = H_TDE[index][:, 0] * 1.E3
            Y = H_TDE[index][:, 1] * 1.E3
            Z = H_TDE[index][:, 2] * 1.E3
            Z = np.sqrt(Z ** 2)

            U_rake_00 = tde.calc_tri_displacements(SX, SY, SZ, X, Y, Z, Poisson_Ratio, -1.0, 0.0, 0.0)
            U_rake_90 = tde.calc_tri_displacements(SX, SY, SZ, X, Y, Z, Poisson_Ratio, 0.0, 0.0, -1.0)

            green_rake_00 = np.zeros((array_gps.shape[0], 5))
            green_rake_90 = np.zeros((array_gps.shape[0], 5))

            green_rake_00[:, 0:2] = array_gps
            green_rake_90[:, 0:2] = array_gps

            green_rake_00[:, 2] = U_rake_00['x']
            green_rake_00[:, 3] = U_rake_00['y']
            # Meade tde convention: positive Uz downward - corrected 18/02/2020 by adding -
            green_rake_00[:, 4] = -U_rake_00['z']

            green_rake_90[:, 2] = U_rake_90['x']
            green_rake_90[:, 3] = U_rake_90['y']
            # Meade tde convention: positive Uz downward - corrected 18/02/2020 by adding -
            green_rake_90[:, 4] = -U_rake_90['z']

            GREEN[index, :, :, 0] = green_rake_00[:, 2:5]
            GREEN[index, :, :, 1] = green_rake_90[:, 2:5]


    else:

        slip = 1.0

        for index in H_fault_xy.keys():
            fault_xy = H_fault_xy[index]
            fault_lp = H_fault_lp[index]

            #    ARRAY_SOURCES[index,:]=SOURCES[index,:]

            #     # add 24/03/2016 for the weird geometry associated with the seamount from JY. Collot (JGR, 2017)
            #     print '!!!! JY',index,fault_lp.x,fault_lp.y,fault_lp.strike,rake
            # #
            #     if fault_lp.strike > 90.0:
            #         print '!!!! ',index,fault_lp.x,fault_lp.y,fault_lp.strike,rake
            #         rake=rake+180.0
            #         print '!!!! corrected to ',index,fault_lp.x,fault_lp.y,fault_lp.strike,rake
            #     if fault_lp.strike < -90.0:
            #         print '!!!! ',index,fault_lp.x,fault_lp.y,fault_lp.strike,rake
            #         rake=rake-180.0
            #         print '!!!! corrected to ',index,fault_lp.x,fault_lp.y,fault_lp.strike,rake

            green_rake_00 = fault_xy.disp_slip_rake(slip, 0.0, array_gps)
            green_rake_90 = fault_xy.disp_slip_rake(slip, 90.0, array_gps)

            #        green_rake_00=fault_xy.disp_slip_rake_no_edcmp(slip,0.0,array_gps)
            #        green_rake_90=fault_xy.disp_slip_rake_no_edcmp(slip,90.0,array_gps)

            green_en = green_rake_00[:, 2:4]

            #    print 'GREEN[index,:,:,0] ',GREEN[index,:,:,0].shape
            #    print 'green_rake_00[:,2:5] ',green_rake_00[:,2:5].shape
            GREEN[index, :, :, 0] = green_rake_00[:, 2:5]
            GREEN[index, :, :, 1] = green_rake_90[:, 2:5]

    return GREEN