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

    import pyeq.coulomb
    import numpy as np


    shear  = pyeq.coulomb.stress2shear(S, strike, dip, rake)
    normal = pyeq.coulomb.stress2normal(S, strike, dip, rake)

    coulomb = shear + mu * normal

    return (coulomb, shear, normal)
