def make_dlo_stf(nfault, nstep):
    import numpy as np
    import scipy.sparse

    dlo_stf = np.zeros((nstep, nfault * nstep))

    # create DLO pattern
    pattern = np.zeros((3 * nfault))
    pattern[:nfault] = 1.
    pattern[nfault:2 * nfault] = -2.
    pattern[2 * nfault:] = 1.

    # fill DLO pattern
    for i in np.arange(1, nstep - 1):
        dlo_stf[i, (i - 1) * nfault:(i + 2) * nfault] = pattern

    # first and last row
    dlo_stf[0, 0:nfault] = -1
    dlo_stf[0, nfault:2 * nfault] = 1

    dlo_stf[-1, (nstep - 2) * nfault:(nstep - 1) * nfault] = 1
    dlo_stf[-1, (nstep - 1) * nfault:(nstep) * nfault] = -1
    #    dlo_stf[-1,nfault:2*nfault] = 1

    # return the associated sparse matrix
    return scipy.sparse.bsr_matrix(dlo_stf)