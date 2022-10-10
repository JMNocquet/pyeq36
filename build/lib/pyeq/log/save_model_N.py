def save_model_N( model):
    """
    Save model as pck and npy file for the normal system
    """

    # import
    from datetime import datetime
    import pickle
    import numpy as np

    import pyeq.message.message as MESSAGE
    import pyeq.message.verbose_message as VERBOSE
    import pyeq.message.error as ERROR
    import pyeq.message.warning as WARNING
    import pyeq.message.debug_message as DEBUG

    # N & Nd file names
    str_time = s1 = datetime.now().strftime("%Y-%m-%d-%H-%M-%S")

    N_fn = ("N_%s" % (str_time))
    Nd_fn = ("Nd_%s" % (str_time))

    # save npy files
    MESSAGE("saving Observation Normal Matrix as %s (%.01lf  Gb)" % ( N_fn,  model.N.nbytes / 1024 / 1024 / 1024 ))
    try:
        np.save(N_fn,model.N)
    except:
        ERROR(("Could not write npy file %s " % N_fn), exit=True)

    MESSAGE("saving Observation Normal Vector as %s (%.01lf  Gb)" % ( Nd_fn,  model.Nd.nbytes / 1024 / 1024 / 1024 ))
    try:
        np.save(Nd_fn,model.Nd)
    except:
        ERROR(("Could not write npy file %s " % Nd_fn), exit=True)

    # save model


    MESSAGE("Deleting Normal System from model")
    try:
        N_save = np.copy(model.N)
        delattr(model, 'N')
        Nd_save = np.copy(model.Nd)
        delattr(model, 'Nd')
    except:
        ERROR("deleting normal system attributes")

    model_fn = ("model_%s.mpck" % (str_time))

    model.N_fn = N_fn+'.npy'
    model.Nd_fn = Nd_fn+'.npy'

    MESSAGE("saving model pickle as %s " % ( model_fn ))

    ofile = open(model_fn, 'wb')
    pickle.dump(    model, ofile)
    ofile.close()

    # re-populate model
    model.N = N_save
    model.Nd = Nd_save

    return model
