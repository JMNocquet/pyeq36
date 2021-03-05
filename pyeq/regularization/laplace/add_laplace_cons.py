def add_laplace_cons( model ):
    """
    Adds Discrete Laplace regularization to the normal system
    """

    DEBUG=True
    model.verbose=True

    # import
    import scipy.sparse
    import pyeq.regularization.laplace
    import pyeq.message.message as MESSAGE
    import pyeq.message.verbose_message as VERBOSE
    import pyeq.message.debug_message as DEBUG
    import pyeq.message.warning as WARNING
    import pyeq.message.error as ERROR

    if float(model.lambda_temporal_smoothing) == 0. and \
       float(model.lambda_spatial_smoothing) == 0. and \
            not model.np_sigma.all():
        WARNING("User provided parameters imply no regularization constraint at all.")
        WARNING("This is only possible with enough observations and very few faults")

        return model

    if float(model.lambda_spatial_smoothing) != 0.:

        # get spatial DLO
        VERBOSE("Making Discrete Laplace Operator for the triangular mesh")
        dlos = pyeq.regularization.laplace.make_dlo_trimesh( model.geometry, stencil=4, verbose=model.verbose)

        # normalize
        VERBOSE("Normalizing Discrete Laplace Operator for the triangular mesh")
        ndlos = pyeq.regularization.laplace.normalize_dlo_trimesh( dlos , float(model.lambda_spatial_smoothing) )
    else:
        ndlos = scipy.sparse.lil_matrix(model.N.shape)
        MESSAGE("No spatial smoothing.")

    if float(model.lambda_temporal_smoothing) != 0.:
        # get temporal DLO
        VERBOSE("Making Discrete Laplace Operator for temporal smoothing")
        dlot = pyeq.regularization.laplace.make_dlo_temporal( model.nfaults, model.nstep)

        # normalize
        VERBOSE("Normalizing Discrete Laplace Operator for temporal smoothing")
        ndlot = pyeq.regularization.laplace.normalize_dlo_temporal( dlot , float(model.lambda_temporal_smoothing) )
    else:
        ndlot = scipy.sparse.lil_matrix(model.N.shape)
        MESSAGE("No temporal smoothing.")

    # apply Laplace regularization
    # case both spatial and temporal smoothing

    VERBOSE("Merging spatial/temporal smoothing constraints and adding it to the normal system")
    model.N += pyeq.regularization.laplace.make_normal_dlo_space_time( ndlos, ndlot)

    # not sure whether this is OK
    model.lambda_spatial_smoothing = float(model.lambda_spatial_smoothing)
    model.lambda_temporal_smoothing = float(model.lambda_temporal_smoothing)

    return model