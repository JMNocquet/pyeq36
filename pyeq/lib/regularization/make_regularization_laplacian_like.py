def make_regularization_laplacian_like( model ):
    """
    Makes a laplacian like regularization, consisting in:
    - adding constraints in the normal system through a discrete spatial average operator
    - adding constraints in the normal system through a time discrete differential operator
    - adding damping constraints in the normal system through the identity equation
    """

    ###########################################################################
    # IMPORT
    ###########################################################################

    import pyeq.lib.regularization
    import os
    import numpy as np
    import pandas
    from datetime import datetime
    from scipy import interpolate

    # reference date is 1980.0
    ref_date_time = datetime(1980, 1, 1, 0, 0, 0)


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
    # PRINT INFO
    ###########################################################################

    print("-- parameters for laplacian like regularization")
    print("   spatial smoothing : %s" % model.lambda_spatial_smoothing )
    print("   temporal smoothing: %s" % model.lambda_temporal_smoothing )
    print("   damping           : %s" % model.lambda_damping )

    ###########################################################################
    # CONVERT TYPES
    ###########################################################################
    try:
        model.lambda_spatial_smoothing = float( model.lambda_spatial_smoothing )
    except:
        if os.path.isfile(model.lambda_spatial_smoothing):

            print("-- lambda_spatial_smoothing provided as file")

            # load file for dates
            d = np.genfromtxt(model.lambda_spatial_smoothing, dtype=str, usecols=(0, 1))
            # join the two columns
            l = []
            for i in np.arange(d.shape[0]):
                l.append(("%s %s" % (d[i, 0], d[i, 1])))
            d = np.array(l)
            # parser
            parse = lambda x: pandas.datetime.strptime(x, '%Y-%m-%d %H:%M:%S')
            # get np_datetime
            np_datetime = np.array(list(map(parse, d)))
            # convert to seconds
            np_date_s = np.array(list(map(int, [x.total_seconds() for x in (np_datetime - ref_date_time)])),
                                 dtype=np.int64)

            # load file for sigma
            np_lambda_spatial_smoothing_s = np.genfromtxt(model.lambda_spatial_smoothing, usecols=(2))


        else:
            # not a file
            print("!!!ERROR lambda_spatial smoothing option error. %s is not a file" % model.lambda_spatial_smoothing)
            sys.exit()

        # make interpolation for each model time step
        np_mid_model_date_s = (model.np_model_date_s[:-1] + np.diff(model.np_model_date_s) / 2.).astype(int)
        f = interpolate.interp1d(np_date_s, np_lambda_spatial_smoothing_s)
        model.lambda_spatial_smoothing = f(np_mid_model_date_s)  # use interpolation function returned by `interp1d`

        print("-- np_lambda_spatial_smoothing_s range: %.2f - %.2f" % (np.min(model.lambda_spatial_smoothing), np.max(model.lambda_spatial_smoothing)))

    model.lambda_temporal_smoothing = float( model.lambda_temporal_smoothing )
    model.lambda_damping = float( model.lambda_damping )


    ###########################################################################
    # SPATIAL
    ###########################################################################

    print("-- Spatial constraints")
    if isinstance(model.lambda_spatial_smoothing,float):
        if model.lambda_spatial_smoothing>0:
            model = pyeq.lib.regularization.add_spatial_constraint_average_operator(model)
    if isinstance(model.lambda_spatial_smoothing,np.ndarray):
        model = pyeq.lib.regularization.add_spatial_constraint_average_operator( model )

    ###########################################################################
    # TEMPORAL
    ###########################################################################

    print("-- Temporal constraints")
    if model.lambda_temporal_smoothing > 0:
        model = pyeq.lib.regularization.add_temporal_constraint_average_operator( model )

    ###########################################################################
    # DAMPING
    ###########################################################################

    print("-- Damping constraints")
    if model.lambda_damping > 0:
        model = pyeq.lib.regularization.add_damping_constraint( model )

    return model