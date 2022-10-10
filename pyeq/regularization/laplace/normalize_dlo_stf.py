def normalize_dlo_temporal( dlo , lam_stf ):
    """
    Normalize a temporal Discrete Laplace Operator matrix for stf

    :param dlo: scipy.sparse.lil 2D numpy array with dimension (nstep x nfault) x (nstep x nfault)
    :param lam_s: weight for spatial smoothing
    """


    # normalize

    return dlo / 2 * 0.5 * lam_stf
