def make_normal_dlo_space_time( dlos, dlot):
    """
    Make the normal system for Discrete Laplace Operator (DLO) for space and time from space and time DLO
    """

    # import
    import numpy as np
    import scipy.sparse

    # get nfault and nstep

    nfault = dlos.shape[0]
    nstep = int( dlot.shape[0] / nfault )

    # duplicate dlos nstep times into dlot

    # dlos lil to dia
    dlos_dia_single = dlos.todia()

    # fills the block diagonal dlos_dia_single into dlot
    dlos_dia = scipy.sparse.block_diag([dlos_dia_single for _ in range(nstep)]).todia()

    diag = -dlos_dia.sum(axis=1).A1 -dlot.sum(axis=1).A1
    dlot.setdiag(diag)

    dlot += dlos_dia

    # normal system
    import pyeq.message.message as MESSAGE
    import pyeq.message.warning as WARNING
    import pyeq.message.error as ERROR
    MESSAGE("Building normal system for regularization from Discrete Laplace Operators. This can take a few minutes.")
    N_dlots = dlot.transpose().dot(dlot)

    return N_dlots