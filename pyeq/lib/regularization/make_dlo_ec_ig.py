def dlo_ec_ig( Dm , d0, cutoff=0.01, sparse=False, verbose=True, debug=False ):
    """
    Computes the discrete Laplace operator on an irregular grid.

    Uses the relationship between a correlation matrix with exponential kernel and the laplacian.

    :param Dm: distance matrix
    :param d0: characteristic lengthscale. Should reflect the mean inter-points distance
    :param cutoff: cutoff value for setting coefficients to zero.
    :param sparse: if True, return the result as a sparse matrix
    :param verbose: verbose mode
    :param debug: print information for debugginh
    :return DLO: the discrete Laplace operator as a 2D numpy array
    """

    # import
    import numpy as np
    import scipy.linalg
    import pyacs.lib.glinalg

    # inverse the correlation matrix (Choleski for improved stability)
    icorr = pyacs.lib.glinalg.syminv( np.exp(-Dm/d0 ))

    # compute the eigen values and vectors
    e,v = scipy.linalg.eigh(icorr, eigvals_only=False)

    # square root of eigen values vector
    sqrt_e = np.sqrt(e)

    # reuilt the true inverse square root
    sqrt_icorr = np.dot( np.dot( v, np.diag(sqrt_e) ) , v.T )

    # check if debug
    if debug:
        print('check')
        print(np.allclose(icorr,np.dot(sqrt_icorr,sqrt_icorr)))

    # computes the sum of elements, row-wise
    diag = np.diagonal( sqrt_icorr )
    S = np.sum(sqrt_icorr,axis=1)

    alpha = S
    best_alpha = np.median(alpha)

    if verbose:
        print("-- expected alpha from analytical solution: %.2lf " % ( 1./np.sqrt(2*np.pi) ) )
        print("-- median estimate for this grid          : %.2lf " % ( best_alpha ) )

    if debug:
        print('alpha')
        print(alpha)

    # compute the Discrete Laplace Operator
    DL = -np.eye(Dm.shape[0]) + 1./best_alpha * sqrt_icorr

    # set negligible elements to zero
    if verbose:
        print("-- number of non-null elements: %d / %d" % (np.count_nonzero(DL),DL.size))
        print("-- max: %.2lf " % ( np.max( DL ) ))
        print("-- max: %.2lf " % (np.min(DL)))

    return DL
