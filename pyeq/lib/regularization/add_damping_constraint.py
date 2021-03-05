def add_damping_constraint(model):
    """
    Add damping constraint.
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

    print("-- Adding damping constraints")
    np.fill_diagonal(model.N, model.N.diagonal() + pyeq.lib.regularization.normalize_lambda(model)[2] )

    return model
