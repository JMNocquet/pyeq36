def normalize_dlo_temporal( dlo , lam_t ):
    """
    Normalize a temporal Discrete Laplace Operator matrix so that
    sum j !=i dlo[i,j] = 1/2 lam_t
    Additionally set dlo[i,i] = 0

    :param dlo: scipy.sparse.lil 2D numpy array with dimension (nstep x nfault) x (nstep x nfault)
    :param lam_s: weight for spatial smoothing
    """

    # set diagonal terms as 0
    dlo.setdiag(0.)

    # normalize

    return dlo / 2 * 0.5 * lam_t