###############################################################################
def unit_strike_vector(strike):
###############################################################################
    """
    Calculates the unit strike vector of a fault
    Coordinates system is x (Easting), y (Northing), z (positive downward)
    Returns a 3-component vector as a 1-D numpy array

    :param strike: fault strike in decimal degrees
    """

    import numpy as np

    rstrike = np.radians(strike)

    us = np.zeros((3))

    us[0] = np.cos(rstrike)
    us[1] = np.sin(rstrike)
    us[2] = 0.0

    return (us)
