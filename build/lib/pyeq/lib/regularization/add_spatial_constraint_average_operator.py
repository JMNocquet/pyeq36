def add_spatial_constraint_average_operator(model , method='laplacian_trimesh'):
    """
    Add spatial regularization constraints.
    Regularization are imposed through an operator indicating that any slip rate on a sub-fault at a given\
    model time step is equal to the weighted average of slip at surrounding sub-faults.

    The current default method is through computing a Discrete Laplace Operator through make_discrete_laplace_trimesh
    Use method = 'weighted_average' for older version
    """

    # NEW APPROACH THOUGH TRUE LAPLACIAN (BEFORE 12/11/2020)
    if method == 'laplacian_trimesh':

        ###############################################################################
        # IMPORT
        ###############################################################################

        import numpy as np
        import scipy.linalg
        import pyeq.lib.regularization
        import pyeq.log

        ###########################################################################
        # INITIALIZE SOME PARAMETERS
        ###########################################################################

        # size of the block matrix associated with slip parameters
        ncm = model.nfaults * model.nstep

        ###############################################################################
        # MAKE THE DISCRETE LAPLACE OPERATOR FOR A SINGLE TIME STEP
        ###############################################################################

        #if model.verbose:
        print("-- making the Discrete Laplace Operator for a single time step")

        DLO = pyeq.lib.regularization.make_discrete_laplace_trimesh( model.geometry )

        ###############################################################################
        # BUILD THE FULL SET OF CONSTRAINTS
        ###############################################################################

        SPATIAL_CONSTRAINT = (DLO.T).dot(DLO).toarray()
        SPATIAL_CONSTRAINT = SPATIAL_CONSTRAINT / np.median( SPATIAL_CONSTRAINT.diagonal() )

        my_lambda = pyeq.lib.regularization.normalize_lambda(model)[0]
        print('my lambda: ', my_lambda)
        print("-- Filling diagonal blocks")
        #model.N[:ncm, :ncm] = model.N[:ncm, :ncm] + scipy.linalg.block_diag(*([ SPATIAL_CONSTRAINT ] * model.nstep))
        # this is an attempt to account for the step duration
        #step_duration_days = np.diff(model.np_model_date_s) / 86400.
        #model.N[:ncm, :ncm] = model.N[:ncm, :ncm] + scipy.linalg.block_diag(*( SPATIAL_CONSTRAINT * np.array( [[step_duration_days]]).T) )
        if isinstance(my_lambda,float):
            print('lambda is float')
            my_lambda = ( np.zeros((model.nstep)) + 1 ) * my_lambda

        print('my lambda: ', my_lambda)

        model.N[:ncm, :ncm] = model.N[:ncm, :ncm] + scipy.linalg.block_diag(*( SPATIAL_CONSTRAINT * np.array( [[my_lambda]]).T) )
        # based on the follwing test using broadcasting
        #A = np.zeros((6, 6))
        #B = np.arange(4).reshape(2, 2)
        #A + scipy.linalg.block_diag(*(B * np.array([[[1, 2, 3]]]).T))

        ###############################################################################
        # HANDLE ORIGIN TIME CONSTANTS
        ###############################################################################
        if model.offset > 0:
            np.fill_diagonal(model.N[-model.nconstant:, -model.nconstant:],
                             np.diagonal(model.N[-model.nconstant:, -model.nconstant:]) + 1. / model.offset ** 2)

        return model


    # OLD APPROACH (BEFORE 12/11/2020)
    if method == 'weighted_average':

        ###############################################################################
        # IMPORT
        ###############################################################################

        import numpy as np
        import scipy.linalg
        import pyeq.lib.regularization
        import pyeq.log

        ###########################################################################
        # INITIALIZE SOME PARAMETERS
        ###########################################################################

        # size of the block matrix associated with slip parameters
        ncm = model.nfaults * model.nstep

        ###############################################################################
        # GET THE INTER SUBFAULT DISTANCE AS A SINGLE SCALAR
        ###############################################################################

        if model.verbose:
            print("-- getting inter subfault distance")

        dc = pyeq.lib.regularization.get_inter_subfault_distance( model.dm, dmin=2, dmax=100, step=0.25, n_contiguous=3)

        if model.verbose:
            print("-- subfault distance (km): %.2lf" % (dc) )

        ###############################################################################
        # MAKE THE SPATIAL AVERAGE OPERATOR FOR A SINGLE TIME STEP
        ###############################################################################

        if model.verbose:
            print("-- making the spatial average operator for a single time step")

        SAO = pyeq.lib.regularization.make_spatial_average_operator( model.dm, dc)

        ###############################################################################
        # BUILD THE FULL SET OF CONSTRAINTS
        ###############################################################################

        SPATIAL_CONSTRAINT = np.dot( SAO.T , SAO )
        SPATIAL_CONSTRAINT = SPATIAL_CONSTRAINT * (model.lambda_spatial_smoothing)

        print("-- Filling diagonal blocks")
        model.N[:ncm, :ncm] = model.N[:ncm, :ncm] + scipy.linalg.block_diag(*([ SPATIAL_CONSTRAINT ] * model.nstep))

        ###############################################################################
        # HANDLE ORIGIN TIME CONSTANTS
        ###############################################################################
        if model.offset > 0:
            np.fill_diagonal(model.N[-model.nconstant:, -model.nconstant:],
                             np.diagonal(model.N[-model.nconstant:, -model.nconstant:]) + 1. / model.offset ** 2)

        return model
