def nnls_predotted(A_dot_A, A_dot_b, tol=1e-8):
    """
    Solve ``argmin_x || Ax - b ||_2`` for ``x>=0``. This version may be superior to the FORTRAN implementation when ``A`` has more rows than
    columns, and especially when ``A`` is sparse.
    Note that the arguments and return values differ from the FORTRAN implementation; in particular, this implementation does not expect the actual
    design matrix ``A`` nor the RHS vector ``b``, but rather ``A.T.dot(A)`` and ``A.T.dot(b)``. These are smaller than the original ``A`` and ``b``
    iff ``A`` has more rows than columns.
    This function also does not return the residual. The squared residual ``|| Ax-b ||^2`` may be calculated efficiently as:
        ``b.dot(b) + x.dot(A_dot_A.dot(x) - 2*A_dot_b)``
    where ``x`` is the output of this function
    Parameters
    ----------
    A_dot_A : ndarray
        Square matrix corresponding to ``A.T.dot(A)``, where ``A`` is as shown above.
    A_dot_b : ndarray
        Vector corresponding to ``A.T.dot(b)``, where ``A`` and ``b`` are as shown above.
    Returns
    -------
    x : ndarray
        Solution vector.
    """

    import numpy as np

    A_dot_A = np.asarray_chkfinite(A_dot_A)
    A_dot_b = np.asarray_chkfinite(A_dot_b)

    if len(A_dot_A.shape) != 2:
        raise ValueError("expected matrix")
    if len(A_dot_b.shape) != 1:
        raise ValueError("expected vector")

    nvar = A_dot_A.shape[0]
    if nvar != A_dot_A.shape[1]:
        raise ValueError("expected square matrix")

    if nvar != A_dot_b.shape[0]:
        raise ValueError("incompatible dimensions")

    P_bool = np.zeros(nvar, np.bool)
    x = np.zeros(nvar, dtype=A_dot_A.dtype)
    s = np.empty_like(x)
    w = A_dot_b
    while not P_bool.all() and w.max() > tol:
        j_idx = w[~P_bool].argmax()
        newly_allowed = np.flatnonzero(~P_bool)[j_idx]
        P_bool[newly_allowed] = True
        s[:] = 0
        currPs = np.flatnonzero(P_bool)
        if len(currPs) > 1:
            s[currPs] = np.linalg.solve(A_dot_A[currPs[:, None], currPs[None, :]], A_dot_b[currPs])
        else:
            currP = currPs[0]
            s[currP] = A_dot_b[currP]/A_dot_A[currP, currP]
        s_P_l_0 = (s[currPs] < 0)
        while s_P_l_0.any():
            currPs_s_P_l_0 = currPs[s_P_l_0]
            alpha = (x[currPs_s_P_l_0]/(x[currPs_s_P_l_0] - s[currPs_s_P_l_0])).min()
            x += alpha*(s-x)
            P_bool[currPs] = (x[currPs] > tol)
            s[:] = 0
            currPs = np.flatnonzero(P_bool)
            if len(currPs) > 1:
                s[currPs] = np.linalg.solve(A_dot_A[currPs[:, None], currPs[None, :]], A_dot_b[currPs])
            else:
                currP = currPs[0]
                s[currP] = A_dot_b[currP]/A_dot_A[currP, currP]
            s_P_l_0 = (s[currPs] < 0)
        x[:] = s[:]
        if x[newly_allowed] == 0:
            break  # avoid infinite loop
        w = A_dot_b - A_dot_A.dot(x)
    return x
