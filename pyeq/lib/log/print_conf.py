def print_conf( model ):
    
    ###########################################################################
    # IMPORT
    ###########################################################################
    
    import shutil
    import sys
    from colors import red
    
    ###########################################################################
    # POULATES CONF DIRECTORY
    ###########################################################################
    
    # conf file
    if model.conf is not None:
        try:
            shutil.copyfile( model.conf , model.odir+'/conf/conf.dat' )
        except:
            print( red("[PYEQ ERROR] Could not copy conf file %s into %s " % ( model.conf , model.odir+'/conf/conf.dat' ) ))
            sys.exit()
            
    # dates
    
    
    
    # sigma