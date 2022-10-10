def make_dlo_temporal(nfault, nstep):
    """
    :nfault: number of sub-faults
    :nstep: number of model time steps

    :note: returns a sparse scipy.sparse.spdiags matrix
    """
    # import

    import scipy as scipy
    import numpy as np

    # Form second difference matrix.
    #e = np.ones((1, nfault * nstep))
    #D = scipy.sparse.spdiags(np.vstack((e, e)), [-nfault, nfault], nfault * nstep, nfault * nstep)

    e1 = np.ones((1, nfault * nstep))
    e2 = np.ones((1, nfault * nstep))

    D = scipy.sparse.spdiags(np.vstack((e1, e2)), [-nfault, nfault], nfault * nstep, nfault * nstep)

    # Fills the diagonal
    diag = -D.sum(axis=1).A1
    D.setdiag(diag)
    # return as sparse matrix
    return D