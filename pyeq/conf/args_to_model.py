def args_to_model( model , args ):
    """
    Fills model information from args
    """

    ###########################################################################
    # IMPORT
    ###########################################################################
    
    # import
    import pyeq.message.message as MESSAGE
    import pyeq.message.verbose_message as VERBOSE
    import pyeq.message.debug_message as DEBUG
    import pyeq.message.error as ERROR

    # import pyacs
    import pyeq.conf
    
    ###########################################################################
    # DEFAULTS
    ###########################################################################

    #model = pyeq.conf.default(model)

    if args.conf is not None:
    ###########################################################################
    # CONF_FILE
    ###########################################################################
        model.conf = args.conf
        MESSAGE(("reading conf file %s" % model.conf ))
        model = pyeq.conf.read_conf(model)

    ###########################################################################
    # ARGS
    ###########################################################################

    for option in dir( args ):
        if (option[0] != '_') and (option != 'conf'):
            option_value = getattr( args, option )
            if option_value is not None:
                MESSAGE("conf changed from command line option -%s = %s " % (option, option_value ) )
                model.__dict__[option] = option_value
    
    ###################################################################
    # WARNING
    ###################################################################
    
    model.warning = '# Warning information\n'

    ###################################################################
    # LAST CHECK
    ###################################################################

    model = pyeq.conf.check_model_conf(model)

    MESSAGE("regularization parameters after check_model_conf: %s %s %s %s %s " % ( model.sigma, model.lambda_spatial_smoothing,
                model.lambda_spatial_smoothing, model.lambda_final_spatial_smoothing, model.lambda_stf ))


    return model 