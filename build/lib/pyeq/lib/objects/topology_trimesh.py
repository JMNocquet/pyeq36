"""
Class for topology of a regular triangular mesh. It is populated through pyeq.lib.geometry.make_topology_trimesh.
When created, a topology_trimesh instance has the following attributes:
- cell_vertex_idx : 2D numpy array of integers (shape = (n,3), n=number of subfaults).
                    Each row are the index of 3 vertices
- vertex_coor     : 2D numpy array of floats (shape = (m,2), m=number of vertices of the whole geometry).
                    Each row has longitude, latitude of vertices.
- vertex_cell_idx: dictionary. vertex_cell_idx[i] gives a list of cells having vertex #i
- cell_neighbours_edge: dictionary. neighbours_edge[i] gives the list of subfaults sharing an edge (= 2 vertices) with subfault #i
- cell_neighbours_vertex = dictionary. neighbours_vertex[i] gives the list of subfaults sharing only one vertex with subfault #i

"""
class topology_trimesh:
  def __init__(self):
    pass
