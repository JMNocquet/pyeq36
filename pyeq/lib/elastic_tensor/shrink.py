
###############################################################################
def shrink(G , lfaults=None, lobs=None, lcomponent=None, lrake=None):
###############################################################################
    """
    
    From a 4-dimension elastic tensor G, creates a 2-dimension matrix G' so that np.dot(G',slip)=d
    where slip is the slip component vector and d the displacement vector
    
    The input G tensor is assumed to be organized as follows:
    
    G(i,j,k,l) is the prediction for dislocation j at site i component k for rake l
    k=0,1,2 = east, north, up
    l=0,1 : rake_principal & rake_conjugate

    The return matrix (2D-numpy array) is organized so that:
    
    the slip vector is partitioned into a main rake component and optionally a conjugate rake component m=[ m_mr , m_cr ]
    the observation displacement vector is organized by components d=[de,dn,du] where de,dn,du are the vectors of E,N,U displacements ( the length of de is the length of lobs )

    :param lfaults   : list of fault indexes to be kept
    :param lobs      : list of observation indexes sites to be kept
    :param lcomponent: list of component indexes to be kept
    :param lrake     : list of rake component indexes to be kept
    
    """

    ###########################################################################
    # import 
    ###########################################################################
    
    import numpy as np
    
    ###########################################################################
    # check input argument
    ###########################################################################
    
    # Dimension of input G
    if G.ndim != 4:
        print("!!! ERROR: Input array MUST have 4 dimensions")
        raise ValueError()
    
    # lfaults
    
    if lfaults is None:
        np_fault = np.arange( G.shape[1] )
    else:
        np_fault = np.unique( lfaults )

    # lobs
    
    if lobs is None:
        np_obs = np.arange( G.shape[0] )
    else:
        np_obs = np.unique( lobs )

    # lcomponent

    if lcomponent is None:
        np_component = np.arange( G.shape[2] )
    else:
        np_component = np.unique( lcomponent )
    
    # lrake

    if lrake is None:
        np_rake = np.arange( G.shape[3] )
    else:
        np_rake = np.unique( lrake )

    ###########################################################################
    # remove unwanted component
    ###########################################################################

    G4 = G[ np_obs][:, np_fault][:,:,np_component][:,:,:,np_rake]


    ###########################################################################
    # convert to numpy 2D array
    ###########################################################################

    # lengths
    nf = np_fault.size
    no = np_obs.size
    nc = np_component.size
    nr = np_rake.size

    # slow approach
    #
    #G = np.zeros( ( no*nc , nf*nr ) )
    #
    #for nobs in np.arange(no):
    #    for nfault in np.arange(nf):
    #        for ncomponent in np.arange(nc):
    #            for nrake in np.arange(nr):
    #                G[ ncomponent*no + nobs , nrake*nf + nfault ] = G4[ nobs, nfault, ncomponent, nrake ]
    
    # fast approach

    G = np.moveaxis(G4,[0,1,2,3],[1,3,0,2]).reshape( no*nc , nf*nr ) 

    return( G ) 
    

    


