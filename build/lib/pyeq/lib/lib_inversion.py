"""
# VARIOUS INVERSION SCHEMES TO RETRIEVE EQ SLIP DISTRIBUTION FROM STATIC OFFSETS
# AUTHOR: J.-M. NOCQUET GEOAZUR-IRD-CNRS-OCA, Univ. Nice, FRANCE
# DATE  : February 2013
"""


#############################################################################################################
# MATRIX AND VECTORS
#############################################################################################################

###############################################################################
def sp(v,W):
###############################################################################
    """
    scalar product with weight.
    
    :param v: vector as 1D numpy array
    :param w: weight vector as 1D numpy array
    
    :return : V.T W V 
    """
    import numpy as np
    
    V=vector_2_array(v)
    return(np.dot(np.dot(V.T,W),V))

###############################################################################
def Sm(G,m,d,Cd,m0,Cm):
###############################################################################
    
    import numpy as np
    
    Wd=np.linalg.inv(Cd)
    LL=np.linalg.cholesky(Cm)
    SQRT_Wm=np.linalg.inv(LL)
    Wm=np.dot(SQRT_Wm.T,SQRT_Wm)
    
    Gm_d=np.dot(G,m)-d

    return(sp(Gm_d,Wd)+sp(m-m0,Wm))

###############################################################################
def Sm_min(G,d,Cd,m0,Cm):
###############################################################################

    import numpy as np
    
    A=np.dot(np.dot(G,Cm),G.T)+Cd
    W=np.linalg.inv(A)
    
    Gm0_d=np.dot(G,m0)-d
    
    return(sp(Gm0_d,W))


###############################################################################
def Corr_to_Cov(corr,sigma_m):
###############################################################################

    """
    Correlation to covariance matrix
    sigma_m is a vector with sqrt(diag(Cov))
    or a scalar
    """
    import numpy as np

    D=np.diag(sigma_m)
    #print 'D ',D
    TEMP=np.dot(D,corr)
    corr=np.dot(TEMP,D)
    return(corr)

###############################################################################
def Cov_to_Corr(Cov):
###############################################################################
    """
    Covariance to correlation transformation
    """
    import numpy as np

    sigma_m=np.sqrt(np.diag(Cov))
    #print 'sigma_m'
    #print sigma_m
    inv_sigma_m=1./sigma_m
    #print inv_sigma_m
    inv_sigma_m=np.diag(inv_sigma_m)
    TEMP=np.dot(inv_sigma_m,Cov)
    #print 'TEMP ',TEMP.shape
    corr=np.dot(TEMP,inv_sigma_m)
    sigma_m=np.sqrt(np.diag(Cov)).reshape(Cov.shape[0],1)
    return(corr,sigma_m)

###############################################################################
def vector_2_array(M):
###############################################################################
    """
    if M is a 1D vector, converts it to a 2D array of single column
    """

    if len(M.shape)==1:return(M.reshape(M.size,1))
    else:return(M)

###############################################################################
def array_2_vector(M):
###############################################################################
    """If M is a 2D array with a single column, return the column as a vector, return M else"""

    if len(M.shape) != 2 or M.shape[1]!=1:return(M) 
    return(M[:,0])

###############################################################################
def numpy_array_2_numpy_recarray(A,names):
###############################################################################
    """
    Converts a numpy array to a numpy recarray
    names is the names of each field
    """
    import numpy as np

    return(np.rec.array(np.core.records.array(list(tuple(A[:,:].T)),dtype={'names':names,'formats':list(map(np.dtype,A[0,:]))})))

###############################################################################
def numpy_recarray_2_numpy_array(A):
###############################################################################
    """
    Converts a structured array (with homogeneous dtype) to a np.array
    """
    import numpy as np
    return(np.array(A.view(A.dtype[0]).reshape(A.shape + (-1,))))
    
#############################################################################################################
# Manipulation pyeq_model objects
#############################################################################################################

#############################################################################################################
def rectangular_dislocations_from_geometry(GEOMETRY):
#############################################################################################################
    """
    Extract the rectangular dislocations from a pyeq_model GEOMETRY numpy array
    """
    import numpy as np
    DISLOCATIONS=np.zeros((GEOMETRY.shape[0],11))
    DISLOCATIONS[:,0]=np.arange(GEOMETRY.shape[0])
    DISLOCATIONS[:,1]=GEOMETRY[:,0]
    DISLOCATIONS[:,2]=GEOMETRY[:,1]
    DISLOCATIONS[:,3]=GEOMETRY[:,2]
    DISLOCATIONS[:,4]=GEOMETRY[:,7]
    DISLOCATIONS[:,5]=GEOMETRY[:,8]
    DISLOCATIONS[:,6]=GEOMETRY[:,3]
    DISLOCATIONS[:,7]=GEOMETRY[:,4]
    DISLOCATIONS[:,8]=GEOMETRY[:,23]
    DISLOCATIONS[:,9]=GEOMETRY[:,9]
    DISLOCATIONS[:,10]=GEOMETRY[:,10]

    return(DISLOCATIONS)

#############################################################################################################
def centroid_from_geometry(GEOMETRY):
#############################################################################################################
    """
    Extract centroid from a pyeq_model GEOMETRY numpy array
    """
    import numpy as np

    SOURCES=np.zeros((GEOMETRY.shape[0],3))
    SOURCES[:,0]=np.arange(GEOMETRY.shape[0])
    SOURCES[:,1]=GEOMETRY[:,11]
    SOURCES[:,2]=GEOMETRY[:,12]

    return(SOURCES)

#############################################################################################################
def make_isotropic_exponential_Cm(Dm,dc,sigma):
#############################################################################################################

    """
    Creates an isotropic exponential matrix Cm_ij=np.exp(-Dm_ij/dc)*sigma**2

    :param Dm: distance matrix between centroid of subfaults, that is Dm[i,j]=distance in km between subfaults i and j
    :param dc: dc, critical ditance in km
    :param sigma: sigma
    """
    import numpy as np

    return(np.exp(-Dm/dc)*sigma**2)


#############################################################################################################
def normalize_linear_system(G,d,m0,Cd,Cm,rake_constraint, insar=False, verbose=True):
#############################################################################################################
    """
    Normalize the linear system using the approach described in Nocquet (2018)
    
    :param G: model matrix
    :param d: observation vector
    :param Dm: distance matrix between centroid of subfaults
    :param Cd: observation covariance matrix
    :param sigma: constraint on a priori model
    :param rake_constraint (0 fixed rake)
    """
#     if debug:
#         print '=>  Saving Matrices'
#         Cm=corrm*args.sigma**2
#         np.save('G_normalized.npy',G)
#         np.save('d_normalized.npy',d)
#         np.save('Cm.npy',Cm)
#         np.save('Cd.npy',Cd)
#         np.save('lower.npy',np.zeros((G.shape[1])))
#         np.save('upper.npy',MAX_SLIP)
    if verbose:
        print('-- Normalizing the linear system')
    
    if verbose:
        print('  -- Inverting data covariance Cd matrix and taking a root-square')

    import numpy as np

    Pd=np.linalg.inv(Cd)
    SQRT_Wd=np.sqrt(Pd)

    print('  -- Dealing with the model covariance matrix Cm')

    if rake_constraint != 0.0:
        
        print('-- Variable rake case')
        print('-- Calculating Cm for the conjugate rake')
        
        if not insar:
            Cm_rake2= rake_constraint**2 * Cm / Cm[0,0]
        else:
            Cm_rake2= rake_constraint**2 * Cm[:-3,:-3] / Cm[0,0]
            
            print('-- Calculating Cm-1/2 for the conjugate rake')
    
        LL_CONJUGATE=np.linalg.cholesky(Cm_rake2)
        SQRT_Wm_CONJUGATE=np.linalg.inv(LL_CONJUGATE)

        if insar:
            print('-- InSAR case')
            LL_insar = np.linalg.cholesky(Cm[-3:,-3:])
            SQRT_Wm_insar = np.linalg.inv(LL_insar) 

            LL=np.linalg.cholesky(Cm[:-3,:-3])
            SQRT_Wm=np.linalg.inv(LL)
        
            SQRT_Wm_ALL=np.zeros((2*SQRT_Wm.shape[0]+3,2*SQRT_Wm.shape[0]+3))

            SQRT_Wm_ALL[0:SQRT_Wm.shape[0],0:SQRT_Wm.shape[0]]=SQRT_Wm
            SQRT_Wm_ALL[SQRT_Wm.shape[0]:2*SQRT_Wm.shape[0],SQRT_Wm.shape[0]:2*SQRT_Wm.shape[0]]=SQRT_Wm_CONJUGATE
            SQRT_Wm_ALL[-3:,-3:] = SQRT_Wm_insar

            SQRT_Wm = SQRT_Wm_ALL
            
            
        else:
            
            LL=np.linalg.cholesky(Cm)
            SQRT_Wm=np.linalg.inv(LL)

            SQRT_Wm_ALL=np.zeros((2*SQRT_Wm.shape[0],2*SQRT_Wm.shape[0]))
            SQRT_Wm_ALL[0:SQRT_Wm.shape[0],0:SQRT_Wm.shape[0]]=SQRT_Wm
            SQRT_Wm_ALL[SQRT_Wm.shape[0]:2*SQRT_Wm.shape[0],SQRT_Wm.shape[0]:2*SQRT_Wm.shape[0]]=SQRT_Wm_CONJUGATE
            
            SQRT_Wm = SQRT_Wm_ALL
        
    else:
        # fixed rake case
        
        if not insar:
            LL=np.linalg.cholesky(Cm)
            SQRT_Wm=np.linalg.inv(LL)
        else:
            LL=np.linalg.cholesky(Cm[:-3,:-3])
            SQRT_Wm=np.linalg.inv(LL)
            LL_insar=np.linalg.cholesky(Cm[-3:,-3:])
            SQRT_Wm_insar = np.linalg.inv(LL_insar) 

            SQRT_Wm_ALL=np.zeros((SQRT_Wm.shape[0]+3,SQRT_Wm.shape[0]+3))
            SQRT_Wm_ALL[0:SQRT_Wm.shape[0],0:SQRT_Wm.shape[0]]=SQRT_Wm
            SQRT_Wm_ALL[-3:,-3:] = SQRT_Wm_insar
            
            SQRT_Wm = SQRT_Wm_ALL
        

    if verbose:
        print('    -- Size of Cm-1 ',SQRT_Wm.shape)
    if verbose:
        print(' -- Modifying the linear system')

    Ad=np.dot(SQRT_Wd,G)
    Am=SQRT_Wm

    A=np.vstack((Ad,Am))
    
    Bd=np.dot(SQRT_Wd,d).flatten()
    Bm=np.dot(SQRT_Wm,m0).flatten()
    
    B=np.hstack((Bd,Bm))

    return(A,B,SQRT_Wm,SQRT_Wd)


###############################################################################
def make_inversion_linear_bayesian(G,d,Cd,Cm,scalar_m0,m0,verbose=True):
###############################################################################
    """
    Simple linear inversion, with model covariance matrix Cm, data covariance matrix Cd, and a priori model m0
    
    m = m0 + ( Cm * Gt * ( (G * Cm * Gt) + Cd )-1 ) * (d - G * m0 )
    
    """

    import numpy as np
    
    CmGt=np.dot(Cm,G.T)
    inv=np.linalg.inv((np.dot(G,CmGt)+Cd))
    TEMP=np.dot(CmGt,inv)

    if scalar_m0 !=0 and isinstance(m0,np.ndarray):
        print("  -- Using alpha=",scalar_m0,"* m0 as apriori model")
        m0=scalar_m0 * m0
        SLIP=m0+np.dot(TEMP,(d-np.dot(G,m0)))
    else:
        print("  -- Using m0=0 as apriori model")
        m0=scalar_m0 * m0
        SLIP=np.dot(TEMP,d)
        Cov_SLIP=Cm - np.dot(np.dot(TEMP,G),Cm)
    return(SLIP)

###############################################################################
def sym_inv_chol(A,debug=False):
###############################################################################
    """
    Computes the inverse of a symmetric positive definite matrix, using Cholesky factorization
    Uses LAPACK DPOTRF (Cholesky factorization) and DPOTRI routines
    """
    import numpy as np
    from scipy.linalg import lapack
    

    (L,info)=lapack.dpotrf(A)
    if debug:print('L' , L)
    (INV,info)=lapack.dpotri(L)
    if debug:print('INV ',INV)
    INV_T=np.copy(INV.T)
    if debug:print('INV_T ',INV_T)
    np.fill_diagonal(INV_T,0.0)
    if debug:print('INV_T ',INV_T)
    if debug:print('INV ',INV)
    return(INV+INV_T)

###############################################################################
def sol_unbounded_gaussian(m0,G,Cm,Cd,d,return_precision_matrix=False,verbose=True,debug=False):
###############################################################################
    """
    Calculates the solution mm and cov(mm) and optionally Hmm of the unbounded inverse problem, using Tarantola & Valette (1982) solution
    """

    import numpy as np
    
    CmGt=np.dot(Cm,G.T)
    if debug:print('CmGt ',CmGt)
    A=np.dot(G,CmGt)+Cd
    if debug:print('A ',A)

    print('numpy.linalg.cond(A, p=None) ', np.linalg.cond(A))

    
    inv_A=sym_inv_chol(A)

    if debug:
        print('A-1 ',inv_A)
        print('A-1 using Linalg ',np.linalg.inv(A))
    
    P=np.dot(CmGt,inv_A)
    
    mm=m0+np.dot(P,(d-np.dot(G,m0)))
    cov_mm=Cm-np.dot(P,CmGt.T)
    
    if debug:
        print('cov_mm method1')
        print(cov_mm)
        print('cov_mm method2')

    Wd=sym_inv_chol(Cd)
    Wm=sym_inv_chol(Cm)
    WCmm=np.dot(np.dot(G.T,Wd),G)+Wm
    cov_mm=sym_inv_chol(WCmm)
    
    if debug:
        print(cov_mm)

    if return_precision_matrix:
        Hmm=np.dot(np.dot(G.T,sym_inv_chol(Cd)),G) + sym_inv_chol(Cm) 
        return(mm,cov_mm,Hmm)
    else:
        return(mm,cov_mm)



###############################################################################
def renormalization_bounded_inversion(A,SLIP_NO_BOUNDS,SLIP_BOUNDS):
###############################################################################
    """
    Calculates the renormalization factor for a bounded domain
    """

    from scipy.stats import mvn
    import numpy as np
    low = SLIP_BOUNDS[:,0]
    print(low)
    upp = SLIP_BOUNDS[:,1]
    print(upp)
    mu = SLIP_NO_BOUNDS
    print('mu')
    print(mu)
    print(A.shape)
    AtA=np.dot(A.T,A)
    print(AtA.shape) 
    print("-- Calculating the a posteriori covariance matrix")
    S=np.linalg.inv(AtA)
    print(S.shape)
    print(low.shape)
    print(upp.shape)
    print(mu.shape)
    p,i = mvn.mvnun(low,upp,mu,S)
    print("-- Integral over the bounded domain ", p,np.sqrt(-2*np.log(p)))
    return(np.sqrt(-2*np.log(p)))


