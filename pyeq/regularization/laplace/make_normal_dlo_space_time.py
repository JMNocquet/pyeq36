def make_normal_dlo_space_time( model , dlos, dlot):
    """
    Make the normal system for Discrete Laplace Operator (DLO) space and time reularization constraints

    """

    # import
    import numpy as np
    import scipy.sparse
    from time import time
    from progress.bar import Bar

    # import
    import scipy.sparse
    import pyeq.regularization.laplace
    import pyeq.message.message as MESSAGE
    import pyeq.message.verbose_message as VERBOSE
    import pyeq.message.debug_message as DEBUG
    import pyeq.message.warning as WARNING
    import pyeq.message.error as ERROR
    import pyacs.debug

    # get nfault and nstep

    nfault = dlos.shape[0]
    nstep = int( dlot.shape[0] / nfault )

    VERBOSE("Merging space and time DLO")

    # duplicate dlos nstep times into dlot
    # dlos lil to dia
    dlos_dia_single = dlos.todia()

    # fills the block diagonal dlos_dia_single into dlot
    dlos_dia = scipy.sparse.block_diag([dlos_dia_single for _ in range(nstep)]).todia()

    VERBOSE("shape DLOS_DIA %s" % str(dlos_dia.shape))
    VERBOSE("shape DLOT %s" % str(dlot.shape))
    diag = -dlos_dia.sum(axis=1).A1 -dlot.sum(axis=1).A1
    dlot.setdiag(diag)

    dlot += dlos_dia

    if model.debug or pyacs.debug():
        DEBUG('DISPLAY DLO S T')
        import matplotlib.pyplot as plt
        plt.imshow(dlot.toarray(), cmap="jet", aspect=1)
        plt.colorbar()
        plt.show()

    # normal system

    # change 09/04/2021 - improves memory use
    # dlot is a sparse matrix and dlot.T dlot still remains in sparse format
    # previous solution was model.N += dlot.transpose().dot(dlot)
    # which converts dlot.T dlot to a dense matrix before adding it.
    # This resulted in doubling memory use
    # The new solution directly uses the indices of non-zero element from the sparse matrix
    # This also appears to be faster
    # I checked that results are model.N remains the same
    # change again 02/05/2021 now using sparse matrix in bsr format which seems to be faster

    from time import time

    OLD = False

    if OLD:

        VERBOSE("DLO matrix to normal DLO")
        t0 = time()
        Ndlot = dlot.transpose().dot(dlot)
        VERBOSE("%.3f s " % (time() - t0))

        VERBOSE("memory usage after building the sparse normal matrix: %.2lf Gb" % pyeq.log.get_process_memory_usage())
        MESSAGE("Finding normal DLO non-zero rows")
        t0 = time()
        rows = sum((m * [k] for k, m in enumerate(np.diff(Ndlot.indptr))), [])
        VERBOSE("%.3f s " % (time() - t0))

        MESSAGE("Adding regularization normal matrix to normal system")
        t0 = time()
        model.N[rows,Ndlot.indices] += Ndlot.data
        VERBOSE("%.3f s " % (time() - t0))

    else:

        VERBOSE("DLO matrix to normal DLO")
        t0 = time()
        Ndlot = dlot.transpose().dot(dlot)
        VERBOSE("%.3f s " % (time() - t0))

        VERBOSE('Sparse csc matrix to to bsr format')
        t0 = time()
        bsrNdlot = Ndlot.tobsr(blocksize=(nfault,nfault))
        VERBOSE("%.3f s " % (time() - t0))

        MESSAGE("Adding regularization normal matrix to normal system")
        t0 = time()
        for i in np.arange(nstep):
            np_step = bsrNdlot.indices[bsrNdlot.indptr[i]:bsrNdlot.indptr[i + 1]]
            for k,j in enumerate(np_step):
                model.N[i*nfault:(i+1)*nfault , j*nfault:(j+1)*nfault] += bsrNdlot.data[bsrNdlot.indptr[i]: bsrNdlot.indptr[i + 1]][k]
        VERBOSE("%.3f s " % (time() - t0))


    # added 30/04/2021
    # DLO and damping constraints on the cumulated slip

    if float( model.lambda_final_spatial_smoothing ) != 0.:
        if float( model.lambda_spatial_smoothing ) != 0.:
            norm_factor = float( model.lambda_final_spatial_smoothing) / float( model.lambda_spatial_smoothing)
        else:
            norm_factor = float( model.lambda_final_spatial_smoothing )

        # re-populate the diagonal with the sum on non-diagonal elements
        diag = -dlos_dia_single.sum(axis=1).A1
        dlos_dia_single.setdiag(diag)

        # TODO Normalization might be wrong
        np_model_step_duration_days = np.diff( model.np_model_date_s) / 86400.

        MESSAGE("Adding smoothing constraints on the cumulated slip distribution: %5.1lf " % float( model.lambda_final_spatial_smoothing) )
        # individual DLO normal matrix
        ata = dlos_dia_single.transpose().dot(dlos_dia_single)*norm_factor / nstep

        # weight to account for model step duration
        MESSAGE("Computing time matrix to account for model step duration")
        TM = np.zeros((nstep,nstep))
        for i in np.arange(nstep):
            TM[i,:(nstep-i)] = np_model_step_duration_days[:(nstep-i)]
        # convert to normal weight for DLO smoothing
        NTMS = np.dot(TM.T,TM)


        # same for damping
        if float(model.sigma_final) != 0:
            MESSAGE("Adding damping constraints on the final slip distribution: %5.1lf " % float(model.sigma_final))
            Id = np.eye(nfault) * 1. / float(model.sigma_final) ** 2
            NTMD = np.dot( np_model_step_duration_days.reshape(-1,1) , np_model_step_duration_days.reshape(1,-1))

        # print progression bar
        bar = Bar('', max=nstep*nstep/2, suffix='%(percent).1f%% - %(eta)ds')

        if float(model.sigma_final) != 0:
            for i in np.arange(nstep):
                for j in np.arange(i + 1):
                    bar.next()
                    # equal weight not accounting for different model step duration
                    #                model.N[i*nfault:(i+1)*nfault,j*nfault:(j+1)*nfault] += (nstep-i) * ata
                    #                model.N[j * nfault:(j + 1) * nfault, i * nfault:(i + 1) * nfault] += (nstep - i) * ata
                    model.N[i * nfault:(i + 1) * nfault, j * nfault:(j + 1) * nfault] += NTMS[i, j] * ata + NTMD[i, j] * Id
                    model.N[j * nfault:(j + 1) * nfault, i * nfault:(i + 1) * nfault] += NTMS[i, j] * ata + NTMD[i, j] * Id

        else:
            for i in np.arange(nstep):
                for j in np.arange(i+1):
                    bar.next()
    # equal weight not accounting for different model step duration
    #                model.N[i*nfault:(i+1)*nfault,j*nfault:(j+1)*nfault] += (nstep-i) * ata
    #                model.N[j * nfault:(j + 1) * nfault, i * nfault:(i + 1) * nfault] += (nstep - i) * ata
                    model.N[i*nfault:(i+1)*nfault,j*nfault:(j+1)*nfault] += NTMS[i,j] * ata
                    model.N[j * nfault:(j + 1) * nfault, i * nfault:(i + 1) * nfault] += NTMS[i,j] * ata

        bar.finish()

    # added 19/05/2021
    # DLO for STF

    if float( model.lambda_stf ) != 0.:
        # this regularization is directly written as a normal equation
        # this is much faster and less memory consuming
        # TODO normalization

        MESSAGE("Adding smoothing constraints on STF: %5.1lf " % float(model.lambda_stf))

        one = np.ones((nfault, nfault)) / 2 * 0.5 * float(model.lambda_stf)**2 / nstep
        pattern = np.hstack([one, -4 * one, 6 * one, -4 * one, one])

        # test for post-seismic - Laplacian constraints cause unrealistic small STF at start

        # original DLO with equality constraint at start/end
        model.N[:nfault, :3 * nfault]           += np.hstack([2 * one, -3 * one, one])
        model.N[nfault:2 * nfault, :4 * nfault] += np.hstack([-3 * one, 6 * one, -4 * one, one])

        model.N[(nstep - 2) * nfault:(nstep - 1) * nfault, -4 * nfault:] += np.hstack([one, -4 * one, 6 * one, -3 * one])
        model.N[(nstep - 1) * nfault:(nstep) * nfault, -3 * nfault:] += np.hstack([one, -3 * one, 2 * one])


        # modified DLO for no temporal constraints at start/end
        #model.N[:nfault, :3 * nfault]           += np.hstack([1 * one, -2 * one, one])
        #model.N[nfault:2 * nfault, :4 * nfault] += np.hstack([-2 * one, 5 * one, -4 * one, one])

        #model.N[(nstep - 2) * nfault:(nstep - 1) * nfault, -4 * nfault:] += np.hstack([one, -4 * one, 6 * one, -3 * one])
        #model.N[(nstep - 1) * nfault:(nstep) * nfault, -3 * nfault:] += np.hstack([one, -3 * one, 2 * one])

        # print progression bar
        bar = Bar('', max=nstep/2, suffix='%(percent).1f%% - %(eta)ds')
        for i in np.arange(2, nstep - 2):

            model.N[i * nfault:(i + 1) * nfault, (i - 2) * nfault:(i + 3) * nfault] += pattern
            bar.next()

        bar.finish()

    #return N_dlots
    return model