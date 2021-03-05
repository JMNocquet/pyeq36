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


    import pyeq.coulomb
    import numpy as np

    t = pyeq.coulomb.stress2traction(S, strike, dip)
    b = pyeq.coulomb.burger(strike, dip, rake)
    shear = np.dot(t.T, b)

    return (shear)
