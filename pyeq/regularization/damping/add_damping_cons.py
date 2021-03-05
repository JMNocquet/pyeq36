def add_damping_cons(model):
    """
    Add damping constraint.
    """



    # NOTES

    ###############################################################################
    # IMPORT
    ###############################################################################

    import numpy as np
    import pyeq.log

    import pyeq.message.message as MESSAGE
    import pyeq.message.verbose_message as VERBOSE
    import pyeq.message.error as ERROR
    import pyeq.message.warning as WARNING
    import pyeq.message.debug_message as DEBUG

    ###########################################################################
    # INITIALIZE SOME PARAMETERS
    ###########################################################################

    if not model.np_sigma.all():
        MESSAGE("No damping constraints. user provided sigma is %.3f" % float(model.sigma))
    else:
        MESSAGE("Adding damping constraints to normal system. Type is %s" % model.sigma_type )
        np.fill_diagonal(model.N, model.N.diagonal() + 1./model.np_sigma**2 )

    return model
