def make_spatial_average_operator(Dm, Dc, weight_type='1/d'):
    """
    Create a spatial average operator from a distance matrix.

    :param Dm: distance matrix
    :param Dc: critical distance so that the average operator only accounts for neighbour cells. Dc can be obtained from get_distance_laplacian.
    :param weight_type: string. weight type. default='1/r'. TODO Additional smoothing kernels to be implemented if required.
    :return: the average operator as a 2D numpy array
    """

    # import
    import numpy as np

    # do the weight
    # 1/r
    if weight_type == '1/d':
        np.fill_diagonal(Dm, 1)
        D = 1. / Dm
    else:
        print("-- weight_type not implemented yet: %s " % (weight_type))
        import sys
        sys.exit()

    # set elements with distance larger than Dc to 0
    D[Dm > Dc] = 0

    # renormalize D so that, for each row, the sum of non diagonal elements is 1 and diagonal is 1
    norm_vector = np.sum(D, axis=1)
    D = -(D.T / norm_vector).T
    np.fill_diagonal(D, 1)

    return D
