"""
Reads a configuration file for pyeq
"""

###############################################################################
def read_conf( model ):
###############################################################################
    """
    Reads a configuration file for pyeq.
    
    """

    import pyeq.message.message as MESSAGE
    import pyeq.message.verbose_message as VERBOSE
    import pyeq.message.debug_message as DEBUG
    import pyeq.message.error as ERROR

    #TODO: the conf file is not actually updated with parameters passed through the command line

    # log conf file
    try:
        lcf = open(model.conf,mode='r')
        model.conf_content = lcf.read()
        lcf.close()
    except:
        ERROR(("could not read %s" % model.conf), exit=True)

    # actually reads conf file
    try:
        cf = open( model.conf , 'r')
    except:
        ERROR(("could not read %s" % model.conf ),exit=True )

    H = {}

    try:
        for line in cf:
            ###############################################################################
            # not to be read
            ###############################################################################
            
            # comment
            if line[0] == '#':continue
            
            # blank line
            if len( line.strip() ) ==0:continue
            
            # incomplete line
            lline = line.split('=')
            if len( lline ) < 2:
                ERROR(("reading file: %s line |%s|" % ( model.conf , line ) ),exit=True)

            key = lline[0].strip()
            # get the value
            # get the position of the first occurrence of = 
            try: 
                idx = line.index('=') 
            except:
                ERROR(("= (equal) character missing for line: %s " % line ),exit=True)

            value = line[idx+1:].strip()
            DEBUG(("reading key in conf : %s" % key))
            if hasattr(model,key) and isinstance(model.__dict__[ key ],str) and (key in ['lexclude_gps','lunweight_gps']):
                model.__dict__[ key ] = model.__dict__[ key ] + ' ' + value
            else:
                model.__dict__[ key ] = value
    except:
        ERROR(("reading %s" % model.conf ),exit=True )
    
    
    return model
        