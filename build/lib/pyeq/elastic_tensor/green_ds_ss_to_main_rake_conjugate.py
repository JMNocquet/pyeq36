"""
Converts a Green tensor of ndim=4 having dip-slip and strike-slip components to main rake and its conjugate
"""


def green_ds_ss_to_main_rake_conjugate(GREEN, GEOMETRY, SGEOMETRY, rake_type, rake_or_pole):
    """

    :param GREEN: original Green tensor (GREEN.ndim=4). See below for format
    :param GEOMETRY: 2D numpy array of pyacs format fault geometry
    :param SGEOMETRY: same as geometry except that geometry is a recarray
    :param rake_type: 'euler','fixed','vector'.
            'euler' means rake is calculated at every subfault from an euler pole
            'fixed' means a single fixed rake for all subfaults
            'vector' means that rake will be provided for every subfault as a 1D numpy array
    :param rake_or_pole: value for rake.
            if rake_type == 'euler' then a string 'lon/lat/w/style' is expected
            if rake_type == 'fixed' a single rake in the range of -180 -- 180 is expected
            if rake_type == 'vector' a 1D numpy array of length GEOMETRY.shape[0] having rake in the range of -180 -- 180 is expected

    :return: the new green tensor, but fault and gps have been switched
    """

    ###########################################################################
    # IMPORT
    ###########################################################################

    import numpy as np
    import pyacs.lib.faultslip
    import pyacs.lib.euler

    import pyeq.message.message as MESSAGE
    import pyeq.message.verbose_message as VERBOSE
    import pyeq.message.error as ERROR
    import pyeq.message.warning as WARNING
    import pyeq.message.debug_message as DEBUG

    # dealing with the principal rake

    RAKE = np.zeros(GEOMETRY.shape[0])
    VEL_FROM_EULER = np.zeros(GEOMETRY.shape[0])

    ## EULER CASE
    if rake_type.lower() == 'euler':

        _tmp, elon, elat, ew, _motion_type = rake_or_pole.split('/')
        #        pole=("%s/%s/%s" % (elon,elat,ew))

        for i in range(GEOMETRY.shape[0]):
            [x, y, strike, dip] = [SGEOMETRY.centroid_long[i], SGEOMETRY.centroid_lat[i], SGEOMETRY.strike[i],
                                   SGEOMETRY.dip[i]]

            RAKE[i] = pyacs.lib.faultslip.rake_from_euler(x, y, strike, dip, rake_or_pole)
            (ve, vn) = pyacs.lib.euler.vel_from_euler(x, y, float(elon), float(elat), float(ew))
            VEL_FROM_EULER[i] = np.sqrt(ve ** 2 + vn ** 2)
        #print("RAKE EULER min max mean %lf %lf %lf " % (np.min(RAKE),np.max(RAKE),np.mean(RAKE)))

    ## CONSTANT FIXED CASE
    elif rake_type.lower() in ('fixed','constant'):
        RAKE = RAKE + float(rake_or_pole)


    ## SUBFAULT VARIABLE RAKE CASE
    else:
        if not isinstance( rake_or_pole , np.ndarray):
            ERROR("fault variable rake option required a 1D numpy array", exit=True)
        if rake_or_pole.shape[0] != GEOMETRY.shape[0]:
            ERROR("rake_or_pole.shape[0] != GEOMETRY.shape[0]", exit=True)

        RAKE = rake_or_pole
        #print("RAKE VECTOR min max mean %lf %lf %lf " % (np.min(RAKE),np.max(RAKE),np.mean(RAKE)))

    RAKE_RADIANS = np.radians(RAKE)
    CONJUGATE_RAKE_RADIANS = np.radians(RAKE + 90.0)

    NEW_GREEN = np.zeros((GREEN.shape[0], GREEN.shape[1], GREEN.shape[2], GREEN.shape[3]))

    GREEN_4GPS_EAST_RAKE_00 = GREEN[:, :, 0, 0]
    GREEN_4GPS_EAST_RAKE_90 = GREEN[:, :, 0, 1]
    GREEN_4GPS_NORTH_RAKE_00 = GREEN[:, :, 1, 0]
    GREEN_4GPS_NORTH_RAKE_90 = GREEN[:, :, 1, 1]
    GREEN_4UP_RAKE_00 = GREEN[:, :, 2, 0]
    GREEN_4UP_RAKE_90 = GREEN[:, :, 2, 1]

    # Now calculating the Green's functions in the principal rake direction

    GREEN_4GPS_EAST_RAKE_PRINCIPAL = np.cos(RAKE_RADIANS) * GREEN_4GPS_EAST_RAKE_00 + np.sin(
        RAKE_RADIANS) * GREEN_4GPS_EAST_RAKE_90
    GREEN_4GPS_NORTH_RAKE_PRINCIPAL = np.cos(RAKE_RADIANS) * GREEN_4GPS_NORTH_RAKE_00 + np.sin(
        RAKE_RADIANS) * GREEN_4GPS_NORTH_RAKE_90

    GREEN_4GPS_EAST_RAKE_CONJUGATE = np.cos(CONJUGATE_RAKE_RADIANS) * GREEN_4GPS_EAST_RAKE_00 + np.sin(
        CONJUGATE_RAKE_RADIANS) * GREEN_4GPS_EAST_RAKE_90
    GREEN_4GPS_NORTH_RAKE_CONJUGATE = np.cos(CONJUGATE_RAKE_RADIANS) * GREEN_4GPS_NORTH_RAKE_00 + np.sin(
        CONJUGATE_RAKE_RADIANS) * GREEN_4GPS_NORTH_RAKE_90

    GREEN_4GPS_UP_RAKE_PRINCIPAL = np.cos(RAKE_RADIANS) * GREEN_4UP_RAKE_00 + np.sin(
        RAKE_RADIANS) * GREEN_4UP_RAKE_90
    GREEN_4GPS_UP_RAKE_CONJUGATE = np.cos(CONJUGATE_RAKE_RADIANS) * GREEN_4UP_RAKE_00 + np.sin(
        CONJUGATE_RAKE_RADIANS) * GREEN_4UP_RAKE_90

    NEW_GREEN[:, :, 0, 0] = GREEN_4GPS_EAST_RAKE_PRINCIPAL
    NEW_GREEN[:, :, 1, 0] = GREEN_4GPS_NORTH_RAKE_PRINCIPAL
    NEW_GREEN[:, :, 2, 0] = GREEN_4GPS_UP_RAKE_PRINCIPAL

    NEW_GREEN[:, :, 0, 1] = GREEN_4GPS_EAST_RAKE_CONJUGATE
    NEW_GREEN[:, :, 1, 1] = GREEN_4GPS_NORTH_RAKE_CONJUGATE
    NEW_GREEN[:, :, 2, 1] = GREEN_4GPS_UP_RAKE_CONJUGATE

    return (NEW_GREEN , RAKE )