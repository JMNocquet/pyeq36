def print_conf( model ):
    """
    print the conf file used for controlling the inversion
    """
    ###########################################################################
    # IMPORT
    ###########################################################################

    import pyeq.message.message as MESSAGE
    import pyeq.message.verbose_message as VERBOSE
    import pyeq.message.error as ERROR
    import pyeq.message.warning as WARNING
    import pyeq.message.debug_message as DEBUG

    import time

    ###########################################################################
    # ADD INFORMATION
    ###########################################################################

    hash_line = ("###########################################################################\n")

    new_header_conf = hash_line
    new_header_conf += ("# conf file by user %s\n" % model.username )
    new_header_conf += ("# %s://%s\n" % (model.hostname,model.conf) )
    new_header_conf += ("# run from %s:%s\n" % (model.hostname, model.wdir))
    str_st = time.strftime("%Y/%m/%d  %H:%M:%S", time.localtime(model.start_time))
    new_header_conf += ("# on %s\n" % (str_st))
    new_header_conf += hash_line

    model.conf_content = new_header_conf + model.conf_content


    ###########################################################################
    # POPULATES CONF DIRECTORY
    ###########################################################################

    try:
        fc = open(model.odir+'/info/conf.dat','w')
        fc.write(model.conf_content)
        fc.close()
    except:
        ERROR("Could not write %s" % model.odir+'/info/conf.dat')

