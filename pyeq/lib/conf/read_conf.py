"""
Reads a configuration file for pyeq
"""

###############################################################################
def read_conf( model ):
###############################################################################
    """
    Reads a configuration file for pyeq.
    
    """
    
    try:
        cf = open( model.conf , 'r')
    except:
        print("!!!ERROR: could not read %s" % ( conf_file ) )
        raise FileError
    
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
                print("!!!ERROR reading file: %s line |%s|" % ( fname , line ) )
                import sys
                sys.exit()

            key = lline[0].strip()
            # get the value
            # get the position of the first occurrence of = 
            try: 
                idx = line.index('=') 
            except ValueError as e: 
                print ("No = character ain line. Error: %s " % line ) 
                import sys
                sys.exit()

            value = line[idx+1:].strip()
            
            model.__dict__[ key ] = value
    except:
        print("!!!ERROR reading %s" % ( model.conf ) )
        raise FileError
    
    
    return model
        