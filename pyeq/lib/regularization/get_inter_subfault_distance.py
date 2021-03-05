def get_inter_subfault_distance(Dm, dmin=2, dmax=100, step=0.25, n_contiguous=3):
    """
    Given a distance matrix of subfaults, get a median distance allowing to select neighbour elements.

    :param Dm: distance matrix. A symetric matrix where D[i,j] is the distance between sub-fault i and j.
    :param dmin: min distance for search. dmin should be greater than step
    :param dmax: max distance for search
    :param step: step in distance search
    :param n_contiguous: number of contiguous elements used for subsequent laplacian operator. Use 3 or 12 for triangles, 4 or 8 for rectangles.
    :return: the optimal distance.
    """

    import numpy as np

    n = 0
    dc = dmin - step
    while n < (n_contiguous + 1):
        dc = dc + 0.5
        D = np.copy(Dm)
        D[D > dc] = 0
        n = int(np.median(np.count_nonzero(D, axis=0)))

    print("-- sub-fault median distance (km): %.2lf" % (dc))

    return dc