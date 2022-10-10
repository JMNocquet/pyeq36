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

    import numpy as np

    rstrike = np.radians(strike)
    rdip = np.radians(dip)

    ud = np.zeros((3))

    ud[0] = np.sin(rstrike) * np.cos(rdip)
    ud[1] = -np.cos(rstrike) * np.cos(rdip)
    ud[2] = -np.sin(rdip)

    return (ud)
