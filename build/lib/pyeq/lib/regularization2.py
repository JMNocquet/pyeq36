"""
Handles regularization constraints
"""


###############################################################################
def add_Cminv_to_ATPA_temporal_smoothing(ATPA,CM,w,delta_time_in_days,sigma_tau,tau):
###############################################################################
    """
    Adds simultaneous spatial and temporal smoothing to ATPA
    
    ATPA: 2D numpy square array of shape (nstep x nfaults, nstep x nfaults)
    Cm_inv: Cm^{-1} the inverse of slip covariance matrix 2D numpy square array of shape (nfaults, nfaults)
    w : a scalar or 1D array of factors applied to Cm_inv
    """
    
    
    import numpy as np
    import pyacs.lib.glinalg

    def time_to_2D_diff_time_in_days( delta_time_in_days ):
        
        delta_day_2D = np.zeros( ( delta_time_in_days.shape[0] , delta_time_in_days.shape[0]  ) )

        print( delta_day_2D.shape )

        for i in np.arange( delta_day_2D.shape[0] ):
            for j in np.arange( i , delta_day_2D.shape[0] ):
                delta_day_2D[i,j] = delta_time_in_days[j] - delta_time_in_days[i]
                delta_day_2D[j,i] = delta_day_2D[i,j] 

        return delta_day_2D
    
    delta_day_2D = time_to_2D_diff_time_in_days( delta_time_in_days )


    nstep = int( ATPA.shape[0] / CM.shape[0] )
    nfaults = CM.shape[0]

    BIG_CM = np.copy( ATPA ) * 0.

    if not isinstance(w,np.ndarray):
        w =  (np.zeros(nstep) + 1. ) * w


    for i in np.arange(nstep): 
        for j in np.arange( i , nstep): 
    
            print(("-- adding spatial & time regularization constraints at step # %04d / %04d over %04d time steps" %(i+1, j+1 , nstep) ))
            STEP_CM = CM * 1/w[i]**2 * np.exp( -delta_day_2D[i,j] / tau )
            BIG_CM[i*nfaults:(i+1)*nfaults,j*nfaults:(j+1)*nfaults] = STEP_CM
            BIG_CM[j*nfaults:(j+1)*nfaults , i*nfaults:(i+1)*nfaults ] = STEP_CM
    
    print('-- Inverting model matrix ')
    ATPA = ATPA + pyacs.lib.glinalg.cov_to_invcov(BIG_CM)
    
    return(ATPA)


###############################################################################
def add_Cminv_to_ATPA(ATPA,Cm_inv,w):
###############################################################################
    """
    Adds Cm^{-1} to ATPA
    
    ATPA: 2D numpy square array of shape (nstep x nfaults, nstep x nfaults)
    Cm_inv: Cm^{-1} the inverse of slip covariance matrix 2D numpy square array of shape (nfaults, nfaults)
    w : a scalar or 1D array of factors applied to Cm_inv
    """
    
    
    import numpy as np

    nstep = int( ATPA.shape[0] / Cm_inv.shape[0] )
    nfaults = Cm_inv.shape[0]

    if not isinstance(w,np.ndarray):
        w =  (np.zeros(nstep) + 1. ) * w


    for i in np.arange(nstep): 
    
        print(("-- adding spatial regularization constraints at step # %04d  over %04d time steps" %(i+1, nstep) ))
        
        ATPA[i*nfaults:(i+1)*nfaults,i*nfaults:(i+1)*nfaults] += Cm_inv * w[i]
        
    return(ATPA)

###############################################################################
def renorm_w_geometry(SGEOMETRY, Dm, Cm_type , dc):
###############################################################################
    """
    calculates a renormalization factor according to the discretization
    
    :param SGEOMETRY: rec array of pyea geometry
    :param Dm: Distance matrix as 2D numpy array
    :param Cm_type: 'exponential' or 'm_exponential'
    :param dc: correlation distance
    """

    import numpy as np
    
    # gets dc0 = minimum characteristic length for discretization
    dc0=np.sqrt(np.min(SGEOMETRY.rdis_area))
    
    # exponential case w= dc0**2 * n / np.sum(np.exp(-Dm/dc))
    if Cm_type=='exponential':
        w= dc0**2 * SGEOMETRY.rdis_area.shape[0]**2 / np.sum(np.exp(-Dm/dc))
    
    # Radiguet case w= dc**2 * n / np.sum(np.exp(-Dm/dc))
    if Cm_type=='m_exponential':
        w= dc**2 * SGEOMETRY.rdis_area.shape[0] / np.sum(np.exp(-Dm/dc))

    return(w)

