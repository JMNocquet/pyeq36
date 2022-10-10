def run_from_pck( model ):
    """
    Runs pyeq from a model.pck file
    """
    
    ###########################################################################
    # import
    ###########################################################################

    import pyeq.log
    import pyeq.log

    import pyeq.message.message as MESSAGE
    import pyeq.message.verbose_message as VERBOSE
    import pyeq.message.error as ERROR
    import pyeq.message.warning as WARNING
    import pyeq.message.debug_message as DEBUG

    print("###############################################################################")
    MESSAGE("RUN FROM MODEL.PCK",level=1)
    print("###############################################################################")
    
    model_pck = pyeq.log.load_model_N(model)
    
    ###########################################################################
    # NORMAL OBSERVATION SYSTEM IN MODE.PCK 
    # WILL REDO REGULARIZATION
    # we need to update all information related to the regularization
    ###########################################################################
    
    #l_regularization_parameters = ['dc','sigma','cm_type','cm_norm','tau','lambda_spatial_smoothing']
    #print('-- updating regularization parameters for model_Nobs.pck')
    #for param in l_regularization_parameters:
    #    print("-- %s = %s" % ( param , str(model.__dict__[ param ]) ))
    #    model_pck.__dict__[ param ] = model.__dict__[ param ]
    #    model_pck.__dict__[ 'no_opt' ] = model.__dict__[ 'no_opt' ]

    # populates with the right attributes

    for attr in model.__dict__:
        try:
            model_pck.__dict__[attr] = model.__dict__[attr]
        except:
            pass

    model_pck.save = 'False'

    ###########################################################################
    # MAKE THE MODEL CURRENT 
    ###########################################################################
    #model = model_pck
    
    return model_pck