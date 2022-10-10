def make_dlo_trimesh( model , stencil=4, verbose=True ):
    """
    Build a Discrete Laplace Operator (DLO) from a geometry file including a regular triangular mesh.
    In this routine, the Discrete Laplace Operator (DLO) is a n x n (n=number of subfaults=geometry.shape[0])
    is defined as followed:

    For stencil=4, that is a 4-points stencil

    DLO[i,j] = 1/3 if subfault i & j share an edge
    DLO[i,j] = 0 else
    DLO[i,i] = sum -DLO[i,j] for i != j

    For stencil=16

    DLO[i,j] = 1/2 * 1/3 if subfault i & j share an edge
    DLO[i,j] = 1/2 * 1/9 if subfault i & j share only a vertice
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

    # check
    if stencil not in [4,16]:
        err_message = ("stencil parameter = %d. Must be either 4 or 16" % stencil )
        ERROR(err_message,exit=True)

    # get topology
    if not hasattr(model,'topology'):
        VERBOSE("Making topology for the triangular mesh")
        model.topology = pyeq.lib.geometry.make_topology_trimesh( model.geometry )

    # compute the discrete_laplace operator (DLO)
    # create a sparse matrix
    DLO = lil_matrix((model.geometry.shape[0], model.geometry.shape[0]))

    # 16-points stencil
    if stencil==16:
        if verbose:
            MESSAGE("Creating Discrete Laplace Operator (DLO) using a 16 points stencil")
        # loop on cells
        for i in np.arange( DLO.shape[0] ):
            DLO[i,model.topology.cell_neighbours_edge[i]]   = 2./3
            DLO[i,model.topology.cell_neighbours_vertex[i]] = 2./9
            DLO[i,i] = - DLO[i,:].sum(axis=1)

    # 4-points stencil
    if stencil==4:
        if verbose:
            MESSAGE("Creating Discrete Laplace Operator (DLO) using a 4 points stencil")
        # loop on cells
        for i in np.arange( DLO.shape[0] ):
            DLO[i,model.topology.cell_neighbours_edge[i]]   = 1./3
            DLO[i,i] = - DLO[i,:].sum(axis=1)

    return DLO
