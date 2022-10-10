def barycentric_interpolation_one_triangle( vertices, np_points ):
    """
    :param vertices: triangles vertices coordinates as a numpy array with shape = (3,2)
    :param np_points: coorrdinates of points for interpolation as numpy array with shape = (n,2)
    """

    # import
    import numpy as np

    # global barycentric matrix
    B = np.zeros( (np_points.shape[0],3 ) )

    for i in np.arange(np_points.shape[0]):
        # individual barycentric matrix
        # see https://en.wikipedia.org/wiki/Barycentric_coordinate_system
        x1,y1 = vertices[0,:]
        x2,y2 = vertices[1,:]
        x3,y3 = vertices[2,:]

        x,y = np_points[i,:]
        X = np.array([1,x,y]).reshape(-1,1)

        A_2 = x1*(y2-y3) + x2*(y3-y1)+x3*(y1-y2)

        M = np.zeros((3,3))
        M[0,0] = x2*y3-x3*y2
        M[1,0] = x3*y1-x1*y3
        M[2,0] = x1*y2-x2*y1

        M[0,1] = y2-y3
        M[1,1] = y3-y1
        M[2,1] = y1-y2

        M[0,2] = x3-x2
        M[1,2] = x1-x3
        M[2,2] = x2-x1

        [l1,l2,l3] = 1./A_2 * np.dot(M,X).flatten()
        B[i,:] = [l1,l2,l3]

    return B