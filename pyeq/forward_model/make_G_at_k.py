def make_G_at_k( model , model_step , W_Sdk , gamma_matrix ):
    """
    Builds all the Green matrix required for observation k
    """
    
    # import
    
    import numpy as np
    
    G = {}
    
    # print
    #print( 'model.G\n' , model.G )
    #print( 'W_Sdk\n' , W_Sdk )
    #print( 'gamma_matrix\n' , gamma_matrix )
    #print( ( W_Sdk.T * gamma_matrix[ : , model_step ]).T.flatten() )
    
    # loop
    
    for row in np.flip( np.arange( model_step +1 ) ):
        #print('row ' , row , gamma_matrix[ : , row ] )
        G[row] = ( model.G.T *  ( W_Sdk.T * gamma_matrix[ : , row ]).T.flatten() ).T
                            
    return G 
            
             
    