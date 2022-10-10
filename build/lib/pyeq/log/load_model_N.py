def load_model_N( model ):
    """
    load model as pck and npy file for the normal system
    """

    # import
    from datetime import datetime
    import pickle
    import numpy as np
    import os

    import pyeq.message.message as MESSAGE
    import pyeq.message.verbose_message as VERBOSE
    import pyeq.message.error as ERROR
    import pyeq.message.warning as WARNING
    import pyeq.message.debug_message as DEBUG

    # load model pck file

    # loads the pck
    MESSAGE("Loading %s (%.2f Gb) " % ( model.mpck, os.path.getsize( model.mpck ) /1024 / 1024 / 1024 ) )
    with open( model.mpck, "rb") as f:
        model = pickle.load( f )
    f.close()
    VERBOSE("model object loaded.")

    # load npy files
    VERBOSE("Loading Observation Normal Matrix and Vector :\n- %s\n- %s " % ( model.N_fn, model.Nd_fn ))
    try:
        model.N = np.load(model.N_fn)
        model.Nd = np.load(model.Nd_fn)
    except:
        ERROR("Loading Observation Normal Matrix and Vector :\n- %s\n- %s " % (model.N_fn, model.Nd_fn),exit=1)

    return model
