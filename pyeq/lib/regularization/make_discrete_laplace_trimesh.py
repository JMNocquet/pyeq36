

def make_discrete_laplace_trimesh( geometry ):
    """
    Build a Discrete Laplace Operator from a geometry file.
    In this routine, the Discrete Laplace Operator (DLO) is a n x n (n=number of subfaults=geometry.shape[0]) is defined as followed:
    DLO[i,j] = -1/2 * 1/3 if subfault i & j share an edge
    DLO[i,j] = -1/2 * 1/9 if subfault i & j share only a vertice
    DLO[i,j] = 0 else
    DLO[i,i] = sum -DLO[i,j] for i != j

    The code first build the topology using pyeq.lib.geometry.make_topology_trimesh

    :param geometry: a 2D numpy array with number of subfaults row and 22 columns.
    return: a nxn matrix as a scipy.sparse.lil_matrix
    """

    # import
    import numpy as np
    from scipy.sparse import lil_matrix

    import pyeq.lib.geometry.make_topology_trimesh

    # get topology
    topology = pyeq.lib.geometry.make_topology_trimesh( geometry )
    # compute the discrete_laplace operator (DLO)
    # create a sparse matrix
    DLO = lil_matrix((geometry.shape[0], geometry.shape[0]))
    # loop on cells
    for i in np.arange( DLO.shape[0] ):
        DLO[i,topology.cell_neighbours_edge[i]]   = 2./3
        DLO[i,topology.cell_neighbours_vertex[i]] = 2./9
        DLO[i,i] = - DLO[i,:].sum(axis=1)

    return DLO
