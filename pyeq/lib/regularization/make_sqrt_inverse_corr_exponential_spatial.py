def make_sqrt_inverse_corr_exponential_spatial( DM , dc , verbose=False ):
    """
    Given a distance matrix NDM, computes the associated inverse square root correlation matrix for assuming an exponential kernel, 
    through Laplacian like differenciating operator.
    
    :param DM: symmetric distance matrix
    :param dc: critical length for correlation
    
    :return: the inverse square root correlation matrix.
    :note: the Laplacian operator used here is crude and will work with quasi-equilateral triangles meshes.
    """
    
    # import
    import numpy as np
    
    
    # Here is a simplication. Since, we do not have a connection matrix, it is assummed that
    # adjacent cells can be selected using a distance criterion dsel only
        

    dsel = np.median(np.sort(DM,axis=0)[2])*1.2
    if verbose:
        print("-- distance used to select adjacent patches is: %10.3lf km " % dsel )
    
    # set all non-adjacent triangles to 0
    L = np.where( DM > dsel , 0 , DM ) 

    EXP = np.exp(-L/dc)
    N_S = np.zeros( ( DM.shape[0], DM.shape[0] ) )
    
    # loop for non-diagonal elements
    for i in np.arange( DM.shape[0] ):
        for j in np.arange( i+1 , DM.shape[0] ):
        # non diagonal are directly the weights with a minus sign
            if L[i,j] != 0:
                N_S[ i , j ] = N_S[ i , j ] - EXP[ i , j ]
    
    # diagonal
    
    # ndim = 2 check that it works for ndim = 3
    # ndim = 2
    
    # SOL 1 - possible
    #    for i in np.arange( L.shape[0] ):
    #        N_S[ i,i ] = - (np.sum(N_S[i,:]) + np.sum( N_S[:,i]) )
    #        N_S[ i,i ] = ( N_S[ i,i ] + 1 ) /np.sqrt( ndim * np.pi )

    # SOL 2
    #N = np.outer( np.diagonal(N_S)*0.+ np.sqrt(np.mean(np.diagonal(N_S))) , np.diagonal(N_S)*0.+np.sqrt(np.mean(np.diagonal(N_S))) )
    #N_S = N_S / np.sqrt( N )                           

    # SOL 3 - intuition from from Tarantola - Square Roots of Exponential Covariance 2008
    # put the diagonal
    for i in np.arange( DM.shape[0] ):
        N_S[ i,i ] = - (np.sum(N_S[i,:]) + np.sum( N_S[:,i]) ) +1.
    # renormalize
    np.fill_diagonal( N_S , np.diagonal(N_S)/N_S[0,0] )
    for i in np.arange( 1 , DM.shape[0] ):
        N_S[ i,i ] = N_S[ i,i ] / np.sqrt( 1 - np.exp(-dsel/dc)**2 )
    
    return N_S