def normalize_dlo_trimesh_resolution( dlo , lam_s , G ):
    """
    Normalize a spatial Discrete Laplace Operator matrix with the resolution matrix diag(GtG) so that
    sum j !=i dlo[i,j] = 1/2 lam_s
    Additionally set dlo[i,i] = 0

    :param dlo: scipy.sparse.lil 2D numpy array with dimension nfault x nfault
    :param lam_s: weight for spatial smoothing
    """
    # import
    import numpy as np
    import pyeq.message.message as MESSAGE
    import pyeq.message.verbose_message as VERBOSE
    import pyeq.message.error as ERROR
    import pyeq.message.warning as WARNING
    import pyeq.message.debug_message as DEBUG

    # set diagonal terms as 0
    dlo.setdiag(0.)

    # normalize DLO so that the maximum pseudo-observation weight is 1/2 lam_s / sum dlo
    sum = dlo.sum(axis=1)
    norm_coeff = np.max( np.fabs(sum))
    dlo = dlo / norm_coeff * 0.5 * lam_s

    s613_dlo = dlo[613,:].sum()
    s067_dlo = dlo[67,:].sum()

    global_weight = np.sum( sum )

    print('max and min pseudo obs sigma ')
    print(np.min(sum))
    print(np.max(sum))
    print('sum of weight ')
    print(np.sum(sum))
    # computes the resolution weight (square root)
    raw_weight = 1./np.sqrt( np.sum( G**2 , axis=0 ) )
    # renormalize weight so that the sum of weight is 1.
    raw_weight = raw_weight / np.sum( raw_weight )

    print('worst subfault',np.argmin(raw_weight))
    print('best subfault',np.argmax(raw_weight))

    MESSAGE("min, max, median for Discrete Laplacian Operator weighting from resolution: %.1E %.1E %.1E" % (np.min(raw_weight),np.max(raw_weight),np.median(raw_weight)))

    print(G.shape)
    print(dlo.shape)
    print('raw_weight %d dlo %d' % (raw_weight.shape[0],dlo.shape[0]) )

    for i in np.arange( raw_weight.shape[0] ):
        dlo[i,:] = dlo[i,:] * raw_weight[i]

    # renormalize so that np.sum(np.sum(dlo,axis=1) remains the same as before applying the resolution matrix
    current_weight = np.sum( dlo.sum(axis=1 ) )
    dlo = dlo / current_weight * global_weight

    sum = dlo.sum(axis=1)
    print('max and min pseudo obs sigma after reso weighting')

    print(np.min(sum))
    print(np.max(sum))
    print('sum of weight ')
    print(np.sum(sum))

    s613_dlor = dlo[613,:].sum()
    s067_dlor = dlo[67,:].sum()

    print('dlo 613 ' , s613_dlo)
    print('dlor 613 ' , s613_dlor)

    print('ratio' , 1./(s613_dlo/s613_dlor))

    print('dlo 067 ' , s067_dlo)
    print('dlor 067 ' , s067_dlor)

    print('ratio' , 1./(s067_dlo/s067_dlor))


    # DEBUG PRINT SMOOTHING WEIGHT


    return dlo