def normalize_dlo_trimesh( dlo , lam_s ):
    """
    Normalize a spatial Discrete Laplace Operator matrix so that
    sum j !=i dlo[i,j] = 1/2 lam_s
    Additionally set dlo[i,i] = 0

    :param dlo: scipy.sparse.lil 2D numpy array with dimension nfault x nfault
    :param lam_s: weight for spatial smoothing
    """
    # import
    import numpy as np

    # set diagonal terms as 0
    dlo.setdiag(0.)

    # normalize
    sum = dlo.sum(axis=1)
    norm_coeff = np.max( np.fabs(sum))

    return dlo / norm_coeff * 0.5 * lam_s