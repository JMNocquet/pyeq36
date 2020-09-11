def normalization_factor( SGEOMETRY , Dc ):

    import numpy as np 
    # gets dc0 = minimum characteristic length for fault discretization
    dc0 = np.sqrt(np.min( SGEOMETRY.rdis_area))

    return  dc0 / Dc 