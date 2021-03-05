def is_triangle_up( A,B,C ):
    """
    check whether vertices of a triangle are ordered so that the vectorial product is upward = counterclockwise seen
    from above

    :param A,B,C: vertices coordinates provided as 3-components 1D numpy array
    """

    AB = B-A
    BC = C-B

    if (AB[0] * BC[1] - AB[1] * BC[0])>0:
        return True
    else:
        return False




