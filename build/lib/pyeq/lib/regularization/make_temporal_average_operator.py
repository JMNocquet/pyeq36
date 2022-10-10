def make_temporal_average_operator(nfault, nstep):
    """
    Note Returns a sparse matrix
    """
    # import

    import scipy as scipy
    import numpy as np

    # Form second difference matrix.
    e = np.ones((1, nfault * nstep))
    D = scipy.sparse.spdiags(np.vstack((e, -2 * e, e)), [0, nfault, 2 * nfault], nfault * nstep, nfault * nstep)

    # return as sparse matrix
    return D