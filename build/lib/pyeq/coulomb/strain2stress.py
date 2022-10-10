###############################################################################
def strain2stress(strain, lam= 0.28758E+11, mu=0.29353E+11):
###############################################################################
    """
    Strain to stress conversion

    strain can be provided as:
    [Exx,Eyy,Ezz,Exy,Exz,Eyz] (np.array(strain).shape=(1,6))

    [[Exx,Eyy,Ezz,Exy,Exz,Eyz]
    [....                   ]
    [Exx,Eyy,Ezz,Exy,Exz,Eyz]] (np.array(strain).shape=(n,6))

    [[Exx,Exy,Exz],
     [Exy,Eyy,Eyz],
     [Exz,Eyz,Ezz]] np.array(strain).shape=(3,3)

     or an array of 3x3 tensor np.array(strain).shape(n,3,3)

    :param strain: stress tensor as 3x3  2D numpy array

    Returns the associated stress tensor
    """

    import pyeq.coulomb
    import numpy as np
    from pyeq.message import error

    E = np.atleast_2d(np.array(strain))

    if E.shape[-1] == 6:
        [Exx, Eyy, Ezz, Exy, Exz, Eyz] = E.T

        Sxx = 2*mu*Exx+lam*(Exx+Eyy+Ezz);
        Syy = 2*mu*Eyy+lam*(Exx+Eyy+Ezz);
        Szz = 2*mu*Ezz+lam*(Exx+Eyy+Ezz);
        Sxy = 2*mu*Exy;
        Sxz = 2*mu*Exz;
        Syz = 2*mu*Eyz;

        return  np.array([Sxx, Syy, Szz, Sxy, Sxz, Syz]).T

    if E.shape[-1] == 3:
        lS = []
        E = E.reshape(-1,3,3)
        for i in np.arange(E.shape[0]):
            S = 2 * mu * E[i]
            np.fill_diagonal(S,S+lam*np.sum(np.diag(E[i])))
            lS.append(S)

        return np.array([lS])

    error(("Bad dimension/shape for strain: %s" % str(strain.shape)),exit=True)
