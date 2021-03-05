def print_model_tensors_shape( model ):

    import textwrap

    terminal_width = 80

    print("--- model.green    shape: " , model.green.shape )
    print("--- model.obs      shape: " , model.obs.shape )
    print("--- model.dm       shape: " , model.dm.shape )
    print("--- model.name_obs      : %d" % model.name_obs.shape[0] )
    print("%s" % ( textwrap.fill(" ".join( model.name_obs) , width=terminal_width ) ))
