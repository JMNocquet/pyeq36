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

    import numpy as np
    import pyeq.coulomb

    rrake = np.radians(rake)

    s = pyeq.coulomb.unit_strike_vector(strike)
    d = pyeq.coulomb.unit_dip_vector(strike, dip)

    b = np.cos(rrake) * s + np.sin(rrake) * d

    return (b)
