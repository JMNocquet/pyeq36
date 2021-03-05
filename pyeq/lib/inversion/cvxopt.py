"""
Performs inversions using CVXOPT and CVXPY.
"""

def inversion_cvxopt( Nm,
                      Nd,
                      DLS,
                      DLT,
                      D,
                      lambda_spatial_smoothing,
                      lambda_temporal_smoothing,
                      lambda_damping,
                      norm_N = 'L2',
                      norm_DLS='L2',
                      norm_DLT='L2',
                      norm_D='L2',
                      m0=0.):
    """
    Performs inversion using cvxopt and cvxpy

    :param Nm: Normal matrix of the direct problem
    :param Nd: Normal observation vector
    :param DLS: Discrete Laplace Operator for spatial smoothing as a scipy sparse matrix
    :param DLT: Discrete Laplace Operator for temporal smoothing as a scipy sparse matrix
    :param D: damping parameter as a scalar
    :param norm_N: norm to minimize 0.5 * ||Gm-d||
    :param norm_DLS: norm to minimize spatial smoothing lambda_spatial_smoothing * ||DLS*m||
    :param norm_DLT: norm to minimize spatial smoothing lambda_temporal_smoothing * ||DLT*m||
    :param norm_d: norm to minimize lambda_damping * ||m-m0||
    :param m0: prior model. scalar, or 1D array of model.nfaults values, or 2D model.nfaults * model.nstep
    """

    # import
    import numpy as np
    import cvxpy as cp
    import scipy as scipy
    import cvxopt as cvxopt


    # renormalize operators
    DLT = np.dot(DLT.T,DLT)
    DLT = DLT / np.median( DLT.diagonal() ) * lambda_temporal_smoothing

    DLS = np.dot(DLS.T,DLS)
    DLS = DLS / np.median( DLS.diagonal() ) * lambda_spatial_smoothing

    D = lambda_damping * D

    # model parameters
    m = cp.Variable(shape= Nm.shape[0] )

    # list of allowed option
    lnorm = [
        ['LS', 'LS', 'LS', 'LS'],
        ['LS', 'L2', 'L2', 'LS'],
        ['L1', 'L2', 'L2', 'L2']
            ]


    # Least-squares case
    if [norm_N,norm_DLS,norm_DLT,norm_D] == ['LS','LS','LS','LS']:

        obj = cp.Minimize( cp.sum_squares( Nm @ m - Nd )
                           + cp.sum_squares( DLS @ m )
                           + cp.sum_squares( DLT @ m )
                           + cp.sum_squares( D @ m ))

    # Least-square for forward model and damping, L2 for smoothing () from Kim et al., 2009
    if [norm_N,norm_DLS,norm_DLT,norm_D] == ['LS','L2','L2','LS']:
        obj = cp.Minimize( cp.sum_squares( Nm @ m - Nd )
                   + cp.norm(DLS * m, 2)
                   + cp.norm(DLT * m, 2)
                   + cp.sum_squares( D @ m ))

    # L1 for forward model and damping, L2 for smoothing and damping () for tests
    if [norm_N, norm_DLS, norm_DLT, norm_D] == ['L1', 'L2', 'L2', 'L2']:
        obj = cp.Minimize(cp.norm(Nm @ m - Nd , 1 )
                          + cp.norm(DLS * m, 2)
                          + cp.norm(DLT * m, 2)
                          + cp.sum_squares(D * m))

    # SETUP PROBLEM WITH NON-NEGATIVITY CONSTRAINTS

    prob = cp.Problem(obj, [ m>= 0 ])

    # ECOS and SCS solvers fail to converge before
    # the iteration limit. Use CVXOPT instead.
    prob.solve(solver=cp.CVXOPT, verbose=True)
    print('Solver status: {}'.format(prob.status))

    # Check for error.
    if prob.status != cp.OPTIMAL:
        raise Exception("Solver did not converge!")

    print("optimal objective value: {}".format(obj.value))

    return m.value