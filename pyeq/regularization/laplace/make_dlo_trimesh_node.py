def make_dlo_trimesh_node( model , verbose=True ):
    """
    Build a Discrete Laplace Operator (DLO) on nodes from a geometry file including a triangular mesh.
    In this routine, the Discrete Laplace Operator (DLO) is a n x n (n=number of nodes)
    is defined as followed:

    DLO[i,j] = 1/n if i,j are connected by an edge
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
    import pyeq.message.message as MESSAGE
    import pyeq.message.verbose_message as VERBOSE
    import pyeq.message.warning as WARNING
    import pyeq.message.error as ERROR

    # get topology
    if not hasattr( model, 'topology'):
        VERBOSE("Making topology for the triangular mesh")
        model.topology = pyeq.lib.geometry.make_topology_trimesh( model.geometry )
    # compute the discrete_laplace operator (DLO)
    # create a sparse matrix
    DLO = lil_matrix((model.nvertices, model.nvertices))

    # loop on cells
    VERBOSE("Making DLO on vertices")
    for i in np.arange( DLO.shape[0] ):
        n = len(model.topology.vertex_neighbours[i])
        for j in model.topology.vertex_neighbours[i]:
            DLO[i,j]   = 1./n
        DLO[i,i] = - DLO[i,:].sum(axis=1)
    #print(DLO)
    return DLO
