def normalize_lambda(model):
    """
    returns the normalized lambda for both laplacian_like and cvxopt regularization
    """

    import numpy as np

    nlambda_spatial_smoothing = model.lambda_spatial_smoothing / np.sqrt(model.nfaults*model.nstep)
    nlambda_temporal_smoothing = model.lambda_temporal_smoothing / np.sqrt(model.nfaults*model.nstep)

    if model.lambda_damping == 0:
        nlambda_damping = 0
    else:
        nlambda_damping = 1./model.lambda_damping**2

    return nlambda_spatial_smoothing,nlambda_temporal_smoothing,nlambda_damping