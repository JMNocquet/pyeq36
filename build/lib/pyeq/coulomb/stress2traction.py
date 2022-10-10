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

    import pyeq.coulomb
    import numpy as np

    n = pyeq.coulomb.unit_normal_vector(strike, dip)
    t = np.dot(S, n.reshape(3, 1))

    return (t)
