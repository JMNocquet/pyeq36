def args_to_model( model , args ):
    """
    Fills model information from args
    """

    ###########################################################################
    # IMPORT
    ###########################################################################
    
    # import
    import numpy as np

    # import pyacs
    import pyeq.lib.conf
    
    ###########################################################################
    # DEFAULTS
    ###########################################################################

    model = pyeq.lib.conf.default( model )

    if args.conf is not None:
    ###########################################################################
    # CONF_FILE
    ###########################################################################
        model.conf = args.conf
        print("-- reading conf file %s" % model.conf )
        model = pyeq.lib.conf.read_conf( model )

    ###########################################################################
    # ARGS
    ###########################################################################

    for option in dir( args ):
        if (option[0] != '_') and (option != 'conf'):
            option_value = getattr( args, option )
            if option_value is not None:
                print("-- conf changed from command line option -%s = %s " % (option, option_value ) )
                model.__dict__[option] = option_value
    
    ###################################################################
    # WARNING
    ###################################################################
    
    model.warning = '# Warning information\n'

    ###################################################################
    # LAST CHECK
    ###################################################################

    model = pyeq.lib.conf.check_model_conf( model )

    return model 