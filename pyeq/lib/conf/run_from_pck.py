def run_from_pck( model ):
    """
    Runs pyeq from a model.pck file
    """
    
    # import
    import os
    import pickle
    import pyeq.lib.log
    import pyeq.lib.regularization
    import sys
    from time import time
    import numpy as np
    

    print("###############################################################################")
    print("RUN FROM MODEL.PCK")
    print("###############################################################################")
    
    # loads the pck
    print("-- Loading %s (%.1f Gb) " % ( model.pck, os.path.getsize( model.pck ) /1024 / 1024 / 1024 ) )
    with open( model.pck, "rb") as f:
        model_pck = pickle.load( f )
    f.close()
    print("-- model object loaded.")

    
    ###########################################################################
    # NORMAL OBSERVATION SYSTEM IN MODE.PCK 
    # WILL REDO REGULARIZATION
    # we need to update all information related to the regularization
    ###########################################################################
    
    l_regularization_parameters = ['dc','sigma','cm_type','cm_norm','tau']
    print('-- updating regularization parameters for model_Nobs.pck')
    for param in l_regularization_parameters:
        print("-- %s = %s" % ( param , str(model.__dict__[ param ]) ))
        model_pck.__dict__[ param ] = model.__dict__[ param ]
        model_pck.__dict__[ 'no_opt' ] = model.__dict__[ 'no_opt' ]
    
    ###########################################################################
    # MAKE THE MODEL CURRENT 
    ###########################################################################
    model = model_pck
    
    return model