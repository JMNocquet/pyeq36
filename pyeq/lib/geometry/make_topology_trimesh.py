def make_topology_trimesh( model ):
    """
    Builds a topology from a pyeq geometry array. Topology is required from computing the Discrete Laplacian Operator
    used for regularization.


    returns a topology_trimesh object with the following attributes:
    - cell_vertex_idx : 2D numpy array of integers (shape = (n,3), n=number of subfaults).
                        Each row are the index of 3 vertices
    - vertex_coor     : 2D numpy array of floats (shape = (m,2), m=number of vertices of the whole geometry).
                        Each row has longitude, latitude of vertices.
    - vertex_cell_idx: dictionary. vertex_cell_idx[i] gives a list of cells having vertex #i
    - cell_neighbours_edge: dictionary. neighbours_edge[i] gives the list of subfaults sharing an edge (= 2 vertices) with subfault #i
    - cell_neighbours_vertex = dictionary. neighbours_vertex[i] gives the list of subfaults sharing only one vertex with subfault #i
    - vertex_neighbours = dictionary. vertex_neighbours[i] gives a list of indices of neighbours vertices (connected by a triangle edge)
    :param geometry: a 2D numpy array with number of subfaults row and 22 columns.
    :return topology: pyeq.lib.geometry.topology_trimesh object.

    """

    # import
    import numpy as np
    import pyeq.lib.objects.topology_trimesh
    from icecream import ic
    import pyeq.message.message as MESSAGE
    import pyeq.message.verbose_message as VERBOSE
    import pyeq.message.error as ERROR
    import pyeq.message.warning as WARNING
    import pyeq.message.debug_message as DEBUG
    import pyacs.debug

    # tolerance
    tol = 1.E-6

    # topology object to be populated
    topology = pyeq.lib.objects.topology_trimesh()

    # get distance matrix
    np_coor_vertices_1 = model.geometry[:,12:14]
    np_coor_vertices_2 = model.geometry[:,15:17]
    np_coor_vertices_3 = model.geometry[:,18:20]
    np_coor_vertices = np.vstack((np_coor_vertices_1,np_coor_vertices_2))
    np_coor_vertices = np.vstack((np_coor_vertices,np_coor_vertices_3))

    from scipy.spatial import distance_matrix
    Dm = distance_matrix(np_coor_vertices,np_coor_vertices)

    # convert distance matrix to index

    np_cell_vertices = np.zeros(( model.geometry.shape[0],3),dtype=int) * np.nan

    idx = 0
    lvertex_coor = []

    for i in np.arange(Dm.shape[0]):
        idx_cell = np.remainder(i, model.geometry.shape[0])
        idx_vertex_in_cell = np.floor_divide(i , model.geometry.shape[0])
        lindex=np.where(Dm[i,:] < tol)[0]
        tindex = np.array([ np.remainder(lindex , model.geometry.shape[0] ), np.floor_divide(lindex , model.geometry.shape[0]) ]).T

        if [idx_cell,idx_vertex_in_cell] != tindex[0,:].tolist():
            np_cell_vertices[tindex[:,0],tindex[:,1]] = np_cell_vertices[tindex[0,0],tindex[0,1]]
        else:
            np_cell_vertices[tindex[:,0],tindex[:,1]] = idx
            lvertex_coor.append([np_coor_vertices[i,:]])
            idx=idx+1

    topology.cell_vertex_idx = np_cell_vertices.astype(np.int64)
    topology.vertex_coor = np.squeeze( np.array(lvertex_coor) )

    # vertex_cell_idx: for a vertex idx, gives a list of cells having this vertex
    topology.vertex_cell_idx = {}

    for i in np.arange( topology.vertex_coor.shape[0] ):
        topology.vertex_cell_idx[i] =  np.argwhere( topology.cell_vertex_idx == i )

    # cell neighbours
    topology.cell_neighbours_edge = {}
    topology.cell_neighbours_vertex = {}
    for i in np.arange( topology.cell_vertex_idx.shape[0] ):
        topology.cell_neighbours_edge[i]   = []
        topology.cell_neighbours_vertex[i] = []

        list_cell_idx = []
        list_cell_idx = topology.vertex_cell_idx[ topology.cell_vertex_idx[i,0]][:,0].tolist()
        list_cell_idx +=  topology.vertex_cell_idx[ topology.cell_vertex_idx[i,1]][:,0].tolist()
        list_cell_idx += topology.vertex_cell_idx[ topology.cell_vertex_idx[i,2]][:,0].tolist()
        list_cell_idx = [item for item in list_cell_idx if item != i]
        list_cell_idx = sorted( list_cell_idx )

        for idx_neighbour in list(set(list_cell_idx)):
            if list_cell_idx.count(idx_neighbour) ==2:
                topology.cell_neighbours_edge[i].append(idx_neighbour)
            if list_cell_idx.count(idx_neighbour) ==1:
                topology.cell_neighbours_vertex[i].append(idx_neighbour)

    # vertex neighbours
    vertex_neighbours = {}

    for i in np.arange( topology.vertex_coor.shape[0] ):
        vertex_neighbours[i] = []
    for i in np.arange( topology.cell_vertex_idx.shape[0] ):
        vertex_neighbours[topology.cell_vertex_idx[i,0]].append( topology.cell_vertex_idx[i,1] )
        vertex_neighbours[topology.cell_vertex_idx[i,0]].append( topology.cell_vertex_idx[i,2] )
        vertex_neighbours[topology.cell_vertex_idx[i,1]].append( topology.cell_vertex_idx[i,0] )
        vertex_neighbours[topology.cell_vertex_idx[i,1]].append( topology.cell_vertex_idx[i,2] )
        vertex_neighbours[topology.cell_vertex_idx[i,2]].append( topology.cell_vertex_idx[i,0] )
        vertex_neighbours[topology.cell_vertex_idx[i,2]].append( topology.cell_vertex_idx[i,1] )
        #print(vertex_neighbours)
    for i in np.arange( topology.vertex_coor.shape[0] ):
        vertex_neighbours[i] = list(set(vertex_neighbours[i]))

    topology.vertex_neighbours = vertex_neighbours

    # DEBUG - DISPLAY
    import matplotlib.pyplot as plt
    # plot cells
    from matplotlib.lines import Line2D
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for i in np.arange( topology.cell_vertex_idx.shape[0] ):
        line = Line2D([topology.vertex_coor[topology.cell_vertex_idx[i],0]], [topology.vertex_coor[topology.cell_vertex_idx[i],1]])
        ax.add_line(line)
    # plot nodes
    nn = 100
    if pyacs.debug():
        ic(vertex_neighbours[nn] )
        ax.scatter([topology.vertex_coor[vertex_neighbours[nn],0]],[topology.vertex_coor[vertex_neighbours[nn],1]],s=30,c='r')
        ax.scatter([topology.vertex_coor[nn,0]],[topology.vertex_coor[nn,1]],s=50,c='k')
        plt.show()



    VERBOSE("Number of subfaults: %d" %  topology.cell_vertex_idx.shape[0] )
    VERBOSE("Number of vertices : %d" %  topology.vertex_coor.shape[0] )
    import sys
    #sys.exit()
    model.topology =  topology
    model.nvertices = model.topology.vertex_coor.shape[0]
    model.nfaults = model.topology.cell_vertex_idx.shape[0]

    return model
