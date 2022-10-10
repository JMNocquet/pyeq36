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

    import numpy as np

    rstrike = np.radians(strike)
    rdip = np.radians(dip)

    n = np.zeros((3))

    n[0] = -np.sin(rstrike) * np.sin(rdip)
    n[1] = np.cos(rstrike) * np.sin(rdip)
    n[2] = -np.cos(rdip)

    return (n)
