"""
Module for Coulomb and stress calculation
"""

import numpy as np
from pyeq.lib import edcmp_from_27 as edcmp


###############################################################################
def fault_normal_vector(strike, dip):
    ###############################################################################
    """
    Calculates the vector normal to a plane defined by its strike and dip
    Coordinates system is x (Easting), y (Northing), z (positive downward)
    Returns a 3-component vector as a 1-D numpy array

    :param strike: fault strike in decimal degrees
    :param dip: fault dip in decimal degrees
    """

    rstrike = np.radians(strike)
    rdip = np.radians(dip)

    n = np.zeros((3))

    n[0] = -np.sin(rstrike) * np.sin(rdip)
    n[1] = np.cos(rstrike) * np.sin(rdip)
    n[2] = -np.cos(rdip)

    return (n)


###############################################################################
def unit_strike_vector(strike):
    ###############################################################################
    """
    Calculates the unit strike vector of a fault
    Coordinates system is x (Easting), y (Northing), z (positive downward)
    Returns a 3-component vector as a 1-D numpy array

    :param strike: fault strike in decimal degrees
    """

    rstrike = np.radians(strike)

    us = np.zeros((3))

    us[0] = np.cos(rstrike)
    us[1] = np.sin(rstrike)
    us[2] = 0.0

    return (us)


###############################################################################
def unit_dip_vector(strike, dip):
    ###############################################################################
    """
    Calculates the unit dip vector of a fault
    Coordinates system is x (Easting), y (Northing), z (positive downward)
    Returns a 3-component vector as a 1-D numpy array

    :param strike: fault strike in decimal degrees
    :param dip: fault dip in decimal degrees
    """

    rstrike = np.radians(strike)
    rdip = np.radians(dip)

    ud = np.zeros((3))

    ud[0] = np.sin(rstrike) * np.cos(rdip)
    ud[1] = -np.cos(rstrike) * np.cos(rdip)
    ud[2] = -np.sin(rdip)

    return (ud)


###############################################################################
def stress2traction(S, strike, dip):
    ###############################################################################
    """
    Calculates the traction along a plane defined by its strike and dip from a stress tensor

    :param S: stress tensor as 3x3  2D numpy array
    :param strike: fault strike in decimal degrees
    :param dip: fault dip in decimal degrees

    Returns a 3-component vector as a 1-D numpy array

    """

    n = fault_normal_vector(strike, dip)
    t = np.dot(S, n.reshape(3, 1))

    return (t)


###############################################################################
def burger(strike, dip, rake):
    ###############################################################################
    """
    Calculates the unit slip or Burger vector defined by strike, dip, rake

    :param strike: fault strike in decimal degrees
    :param dip: fault dip in decimal degrees
    :param rake: slip rake in decimal degrees

    Returns a 3-component vector as a 1-D numpy array

    """
    rrake = np.radians(rake)

    s = unit_strike_vector(strike)
    d = unit_dip_vector(strike, dip)

    b = np.cos(rrake) * s + np.sin(rrake) * d

    return (b)


###############################################################################
def stress2shear(S, strike, dip, rake):
    ###############################################################################
    """
    Calculates the shear stress from a stress tensor for slip in strike,dip,rake direction

    :param S: stress tensor as 3x3  2D numpy array
    :param strike: fault strike in decimal degrees
    :param dip: fault dip in decimal degrees
    :param rake: slip rake in decimal degrees

    Returns a float
    """

    t = stress2traction(S, strike, dip)
    b = burger(strike, dip, rake)
    shear = np.dot(t.T, b)

    return (shear)


###############################################################################
def stress2normal(S, strike, dip, rake):
    ###############################################################################
    """
    Calculates the normal stress from a stress tensor for slip in strike,dip,rake direction

    :param S: stress tensor as 3x3  2D numpy array
    :param strike: fault strike in decimal degrees
    :param dip: fault dip in decimal degrees
    :param rake: slip rake in decimal degrees

    Returns a float
    """

    t = stress2traction(S, strike, dip)
    n = fault_normal_vector(strike, dip)

    tn = np.dot(t.T, n)

    return (tn)


###############################################################################
def stress2coulomb(S, strike, dip, rake, mu):
    ###############################################################################
    """
    Calculates the Coulomb stress for slip in strike,dip,rake direction

    :param S: stress tensor as 3x3  2D numpy array
    :param strike: fault strike in decimal degrees
    :param dip: fault dip in decimal degrees
    :param rake: slip rake in decimal degrees
    :param mu: friction coefficient (usually in the range 0.1-0.8)

    Returns a float
    """

    shear = stress2shear(S, strike, dip, rake)
    normal = stress2normal(S, strike, dip, rake)

    coulomb = shear + mu * normal

    return (coulomb, shear, normal)


###############################################################################
def coulomb_for_single_receiver_edcmp(rec_long, rec_lat, rec_depth, rec_strike, rec_dip, rec_rake, \
                                      fault_long, fault_lat, fault_depth, fault_length, fault_width, fault_strike,
                                      fault_dip, fault_rake, fault_slip, mu, coor_type='geo'):
    ###############################################################################
    """Calculates the Coulomb stress tensor at a given point induced by a single dislocation
       Uses edcmp

    if coor_typpe='geo' (default), parameters units are decimal degrees for long, lat coordinates and dip/rake and km for depth, length, and width
    if coor_typpe='xyz' parameters units are km for long, lat, decimal degrees for dip/rake and km for depth, length, and width.

    :param rec_long:
    :param rec_lat:
    :param rec_depth:
    :param rec_strike:
    :param rec_dip:
    :param rec_rake:
    :param fault_long:
    :param fault_lat:
    :param fault_depth:
    :param fault_length:
    :param fault_dip:
    :param fault_strike:
    :param fault_slip:
    :param fault_rake:

    returns the Coulomb stress

    """
    # print "--- in Coulomb"
    # print 'rec_long,rec_lat,rec_depth,rec_strike,rec_dip,rec_rake'
    # print rec_long,rec_lat,rec_depth,rec_strike,rec_dip,rec_rake
    # print '---------------------'
    # print 'fault_long,fault_lat,fault_depth,fault_length,fault_width,fault_strike,fault_dip,fault_rake, fault_slip'
    # print fault_long,fault_lat,fault_depth,fault_length,fault_width,fault_strike,fault_dip,fault_rake, fault_slip

    STRESS_REC = stress_tensor_for_single_receiver_edcmp(rec_long, rec_lat, rec_depth, rec_strike, rec_dip, rec_rake, \
                                                         fault_long, fault_lat, fault_depth, fault_length, fault_width,
                                                         fault_strike, fault_dip, fault_rake, fault_slip,
                                                         coor_type=coor_type)

    STRESS = edcmp.get_stress_tensor(STRESS_REC, 0)

    coulomb, shear, normal = stress2coulomb(STRESS, rec_strike, rec_dip, rec_rake, mu)
    # print 'Coulomb in coulomb_from_single_disloc_edcmp ', coulomb
    return (coulomb, shear, normal)


###############################################################################
def stress_tensor_for_single_receiver_edcmp(rec_long, rec_lat, rec_depth, rec_strike, rec_dip, rec_rake, \
                                            fault_long, fault_lat, fault_depth, fault_length, fault_width, fault_strike,
                                            fault_dip, fault_rake, fault_slip, coor_type='geo'):
    ###############################################################################
    """Calculates the stress tensor at a given point induced by a single dislocation
       Uses edcmp

    if coor_typpe='geo' (default), parameters units are decimal degrees for long, lat coordinates and dip/rake and km for depth, length, and width
    if coor_typpe='xyz' parameters units are km for long, lat, decimal degrees for dip/rake and km for depth, length, and width.

    :param rec_long:
    :param rec_lat:
    :param rec_depth:
    :param rec_strike:
    :param rec_dip:
    :param rec_rake:
    :param fault_long:
    :param fault_lat:
    :param fault_depth:
    :param fault_length:
    :param fault_dip:
    :param fault_strike:
    :param fault_slip:
    :param fault_rake:

    returns a 3x3 2D numpy array.
    The stress components are ordered to be used with the edcmp convention:
    (x=Northing,y=Easting,z=depth aka positive downward)
    component are assumed to be in Pa if slip has been provided in meters.

    """

    # converts geographical coordinates to local cartesian coordinates with the edcmp convention
    # (x=Northing,y=Easting,z=depth aka positive downward), all units in meters

    # print '  -- coor_type ',coor_type

    if coor_type == 'geo':

        from pyacs.lib import coordinates as Coordinates

        x, y = Coordinates.geo2flat_earth(rec_long, rec_lat)

        (rec_x, rec_y) = (y, x)
        rec_depth = rec_depth

        (x, y) = Coordinates.geo2flat_earth(fault_long, fault_lat)
        (fx, fy) = (y, x)
        fslip = fault_slip
        fdepth = fault_depth
        flength = fault_length
        fwidth = fault_width
        fstrike = fault_strike
        fdip = fault_dip
        frake = fault_rake

    elif coor_type == 'xyz':

        (rec_x, rec_y) = (rec_lat, rec_long)
        rdepth = rec_depth
        (fx, fy) = (fault_lat, fault_long)
        fslip = fault_slip
        fdepth = fault_depth
        flength = fault_length
        fwidth = fault_width
        fstrike = fault_strike
        fdip = fault_dip
        frake = fault_rake
    else:
        import sys
        print('!!! ERROR: coor_type MUST be either geo or xyz')
        sys.exit()

    # converts all entries to np.array

    (fault_long, fault_lat, fault_depth, fault_length, fault_width, fault_strike, fault_dip, fault_rake, fault_slip) = \
        map(edcmp.to_1Dnparray, \
            [fault_long, fault_lat, fault_depth, fault_length, fault_width, fault_strike, fault_dip, fault_rake,
             fault_slip])

    edcmp.make_edcmp_input_file(rec_x, rec_y, rec_depth, fslip, fx, fy, fdepth, flength, fwidth, fstrike, fdip, frake, \
                                DISP=False, STRAIN=False, STRESS=True, TILT=False, elambda=0.32000E+11, emu=0.32000E+11)

    edcmp.run()

    STRESS_TENSOR_REC = edcmp.get_stress_from_edcmp_output()

    return (STRESS_TENSOR_REC)

















