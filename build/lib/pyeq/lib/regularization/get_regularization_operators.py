def get_regularization_operators( model ):
    """
    Get regularization operator for cvxopt optimization
    """

    ###########################################################################
    # IMPORT
    ###########################################################################

    import pyeq.lib.regularization
    import scipy.sparse
    import numpy as np

    # Form second difference matrix.
    e = np.ones((1, model.nfaults * model.nstep))
    D = scipy.sparse.spdiags(np.vstack((e, -2 * e, e)), [0, model.nfaults, 2 * model.nfaults], model.nfaults * model.nstep, model.nfaults * model.nstep)

    ###########################################################################
    # CHECK ATTRIBUTES EXIST
    ###########################################################################

    check_attr = ['lambda_spatial_smoothing','lambda_temporal_smoothing','lambda_damping']

    for attr in check_attr:

        if not hasattr( model ,  attr ):
            print("ERROR: model has no option %s. Exiting." % (attr) )
            import sys
            exit()


    ###########################################################################
    # CONVERT TYPES AND RENORMALIZE BY NFAULTS
    ###########################################################################

    model.lambda_spatial_smoothing = float( model.lambda_spatial_smoothing )
    model.lambda_temporal_smoothing = float( model.lambda_temporal_smoothing )
    model.lambda_damping = float( model.lambda_damping )

    ###########################################################################
    # PRINT INFO
    ###########################################################################

    print("-- parameters for laplacian like regularization")
    print("   spatial smoothing : %.2lf " % model.lambda_spatial_smoothing )
    print("   temporal smoothing: %.2lf " % model.lambda_temporal_smoothing )
    print("   damping           : %.2lf " % model.lambda_damping )

    ###########################################################################
    # SPATIAL
    ###########################################################################

    print("-- Get Discrete Laplace Spatial Operator")
    DLS = pyeq.lib.regularization.make_discrete_laplace_trimesh( model.geometry )
    model.DLS = scipy.sparse.lil_matrix( scipy.linalg.block_diag(*([DLS.toarray()] * model.nstep))).todia()

    ###########################################################################
    # TEMPORAL
    ###########################################################################

    print("-- Get Discrete Laplace Temporal Operator")
    model.DLT = pyeq.lib.regularization.make_temporal_average_operator(model.nfaults, model.nstep)

    ###########################################################################
    # DAMPING
    ###########################################################################

    print("-- Get Damping Operator")
    model.D = scipy.sparse.eye( model.nfaults * model.nstep )

    #print(model.DLS.shape)
    #print(model.DLT.shape)
    #print(model.D.shape)

    return model