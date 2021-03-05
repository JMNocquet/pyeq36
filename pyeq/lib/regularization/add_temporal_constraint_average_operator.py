def add_temporal_constraint_average_operator(model):
    """
    Add spatial regularization constraints.
    Regularization are imposed through a time-difference operator
    """



    # NOTES

    ###############################################################################
    # IMPORT
    ###############################################################################

    import numpy as np
    import pyeq.lib.regularization
    import pyeq.log

    ###########################################################################
    # INITIALIZE SOME PARAMETERS
    ###########################################################################

    # size of the block matrix associated with slip parameters
    ncm = model.nfaults * model.nstep

    ###############################################################################
    # MAKE THE DIFFERENCE OPERATOR
    ###############################################################################

    if model.verbose:
        print("-- making the time average operator")

    TDO = pyeq.lib.regularization.make_temporal_average_operator(  model.nfaults , model.nstep )

    ###############################################################################
    # BUILD THE FULL SET OF CONSTRAINTS
    ###############################################################################

    TEMPORAL_CONSTRAINT = TDO.T @ TDO
    TEMPORAL_CONSTRAINT = TEMPORAL_CONSTRAINT / np.median( TEMPORAL_CONSTRAINT.diagonal() ) * pyeq.lib.regularization.normalize_lambda(model)[1]

    print("-- Adding temporal smoothing constraints")
    model.N[:ncm, :ncm] = model.N[:ncm, :ncm] + TEMPORAL_CONSTRAINT

    return model
