"""
Wrapper for NNLS & BVLS inversion routines
Added CVXOPT inversion 15/11/2020
"""

###############################################################################
# CVXOPT
###############################################################################

def pyeq_cvxopt( model ):
    """
    performs an inversion using cvxopt
    """

    from time import time
    import pyeq.lib.inversion
    import pyeq.lib.regularization

    t0 = time()

    SLIP = pyeq.lib.inversion.inversion_cvxopt( model.N,
                     model.Nd,
                     model.DLS,
                     model.DLT,
                     model.D,
                     pyeq.lib.regularization.normalize_lambda(model)[0],
                     pyeq.lib.regularization.normalize_lambda(model)[1],
                     pyeq.lib.regularization.normalize_lambda(model)[2],
                     model.norm_N,
                     model.norm_DLS,
                     model.norm_DLT,
                     model.norm_D,
                     model.m0)

    dtime = time() - t0
    print(("-- cvxopt duration (s): %.1lf" % dtime))

    return SLIP , dtime

###############################################################################
# NNLS
###############################################################################

def pyeq_nnls(ATPA, ATPB, nnls, verbose=False):
    """
    solves the normalized linear system ATPA m = ATPB using non-negative least-squares
    This is a wrapper to various solvers. See notes.
    
    :param ATPA: m x m normal matrix as a numpy 2D array
    :param ATPB: 1D numpy array of normalized right-hand-side vector
    :param nnls: string, routine used for nnls. See note
    :param verbose: verbose mode
    :return X, time: solution vector as a 1D numpy array of lenght m, time of inversion
    
    :note: available methods:    
    scipy_nnls : the classical Lawson and Hansen (1974) nnls. Slow for large system
    nnls_predotted: from https://github.com/nnls - Alexfields
    lsqnonneg: python transcription of Matlab lsqnonneg, usually very slow 
    nnlsm_block_pivot: the preferred method on 02/03/2018
    ref: J. Kim and H. Park, Fast nonnegative matrix factorization: An active-set-like method and comparisons,
    SIAM Journal on Scientific Computing, vol. 33, no. 6, pp. 3261-3281, 2011.
    nnlsm_activeset: J. Kim and H. Park, Fast nonnegative matrix factorization: An active-set-like method and comparisons,
    SIAM Journal on Scientific Computing, vol. 33, no. 6, pp. 3261-3281, 2011.
    scipy_lsq_linear: scipy linear least-squares routine
    scipy_lbfgsb: scipy lbfgs routine)

    """
    # import
    import numpy as np
    import datetime
    from time import time

    import pyeq.message.message as MESSAGE
    import pyeq.message.verbose_message as VERBOSE
    import pyeq.message.error as ERROR
    import pyeq.message.warning as WARNING
    import pyeq.message.debug_message as DEBUG


    # SLIP
    SLIP = None

    # NNLS ALGORITHM
    SCIPY_NNLS = False
    SCIPY_BVLS = False
    SCIPY_TRF = False
    NNLS_PREDOTTED = False
    LSQNONNEG = False
    NNLSM_BLOCK_PIVOT = False    
    NNLSM_ACTIVESET = False
    SCIPY_LSQ_LINEAR = False
    SCIPY_LBFGSB = False
    SCIPY_CHOL = False

    MESSAGE("nnls algorithm: %s" % nnls )

    # NNLS ALGORITHM
    if nnls == 'scipy_nnls':
        SCIPY_NNLS = True
    if nnls == 'scipy_bvls':
        SCIPY_BVLS = True
    if nnls == 'scipy_trf':
        SCIPY_TRF = True
    if nnls == 'nnls_predotted':
        NNLS_PREDOTTED = True
    if nnls == 'lsqnonneg':
        LSQNONNEG = True
    if nnls == 'nnlsm_block_pivot':
        NNLSM_BLOCK_PIVOT = True    
    if nnls == 'nnlsm_activeset':
        NNLSM_ACTIVESET = True
    if nnls == 'scipy_lsq_linear':
        SCIPY_LSQ_LINEAR = True
    if nnls == 'scipy_lbfgsb':
        SCIPY_LBFGSB = True
    if nnls == 'scipy_chol':
        SCIPY_CHOL = True


    if SCIPY_NNLS:
        import scipy.optimize
        MESSAGE("Now doing the inversion using scipy nnls "+datetime.datetime.now().strftime("%Y-%m-%d %H:%M"))
        t0 = time()
        SLIP,w = scipy.optimize.nnls(ATPA, ATPB)
        dtime=time()-t0
        print(("-- scipy_nnls duration (s): %.1lf" % dtime))

    if SCIPY_BVLS:
        import scipy.optimize
        MESSAGE("Now doing the inversion using scipy bvls "+datetime.datetime.now().strftime("%Y-%m-%d %H:%M"))
        t0 = time()
        bounds = ( np.zeros( ATPA.shape[0] ) , np.ones( ATPA.shape[0] ) * np.inf )

        verbose = 0
        RES= scipy.optimize.lsq_linear(ATPA, ATPB, bounds=bounds, method='bvls', tol=1e-10, lsq_solver=None, lsmr_tol=None, max_iter=None, verbose= verbose )
        dtime=time()-t0
        MESSAGE(("scipy_bvls duration (s): %.1lf" % dtime))
        SLIP = RES.x
    
    if NNLS_PREDOTTED:
        MESSAGE("Now doing the inversion using nnls_predotted "+datetime.datetime.now().strftime("%Y-%m-%d %H:%M"))
        import pyeq.optimization.nnls.nnls_predotted
        t0 = time()
        SLIP= pyeq.optimization.nnls.nnls_predotted.nnls_predotted(ATPA, ATPB)
        dtime=time()-t0
        MESSAGE(("nnls_predotted duration (s): %.1lf" % dtime))
    
    if LSQNONNEG:
        MESSAGE("Now doing the inversion using lsqnonneg "+datetime.datetime.now().strftime("%Y-%m-%d %H:%M"))
        import pyeq.optimization.nnls.lsqnonneg
        t0 = time()
        (SLIP, resnorm, residual)= pyeq.optimization.nnls.lsqnonneg.lsqnonneg(ATPA, ATPB)
        dtime=time()-t0
        MESSAGE(("lsqnonneg duration (s): %.1lf" % dtime))
    
    if NNLSM_BLOCK_PIVOT:
        MESSAGE(("Now doing the inversion using nnlsm_blockpivot "+datetime.datetime.now().strftime("%Y-%m-%d %H:%M")))
        import pyeq.optimization.nnls.nnlsm
        t0 = time()
        SLIP,_ = pyeq.optimization.nnls.nnlsm.nnlsm_blockpivot(ATPA, ATPB.reshape(-1, 1), is_input_prod=True, init=None)
        SLIP = SLIP.flatten()
        dtime=time()-t0
        MESSAGE(("nnlsm_blockpivot duration (s): %.1lf" % dtime))
        
        
    if NNLSM_ACTIVESET:
        MESSAGE(("Now doing the inversion using nnlsm_blockpivot "+datetime.datetime.now().strftime("%Y-%m-%d %H:%M")))
        import pyeq.optimization.nnls.nnlsm
        t0 = time()
        SLIP,_ = pyeq.optimization.nnls.nnlsm.nnlsm_activeset(ATPA, ATPB.reshape(-1, 1), is_input_prod=True, init=None)
        SLIP = SLIP.flatten()
        dtime=time()-t0
        MESSAGE(("nnlsm_activeset duration (s): %.1lf" % dtime))
    
    if SCIPY_LSQ_LINEAR:
        MESSAGE(("Now doing the inversion using scipy_lsq_linear "+datetime.datetime.now().strftime("%Y-%m-%d %H:%M")))
        from scipy.optimize import lsq_linear
        import numpy as np
        t0 = time()
        RES = lsq_linear(ATPA, ATPB, bounds=(0, +np.inf), method='trf', tol=1e-15, lsq_solver='lsmr', lsmr_tol=None, max_iter=None, verbose=0)
        SLIP = RES.x
        SLIP = SLIP.flatten()
        dtime=time()-t0
        MESSAGE(("scipy_lsq_linear duration (s): %.1lf" % dtime))
    
    if SCIPY_LBFGSB:
        MESSAGE("Now doing the inversion using scipy_lbfgsb "+datetime.datetime.now().strftime("%Y-%m-%d %H:%M"))

        import numpy as np
                
        def fun(x,ATA,ATB):
            R = ATB - np.dot(ATA,x)
            print(np.dot(R.T,R))
            return np.dot(R.T,R)
    
        from scipy.optimize import minimize
        
        BNDS=[]
        for i in np.arange(ATPA.shape[0]):
            BNDS.append( (0.,None) )
        
        BNDS=tuple(BNDS)
        from scipy import linalg
    
        t0 = time()
        
        x0 = linalg.solve(ATPA, ATPB, sym_pos = True)
        
        RES = minimize(fun, x0, args=(ATPA,ATPB), method='L-BFGS-B', \
                                      jac=None, bounds=BNDS, tol=None, callback=None, \
                                      options={'disp': None, 'maxls': 20, 'iprint': -1, 'gtol': 1e-04, 'eps': 1e-08, \
                                               'maxiter': 15000, 'ftol': 2.220446049250313e-04, 'maxcor': 10, 'maxfun': 1500000})
        
        dtime=time()-t0
        SLIP = RES.x
        SLIP = SLIP.flatten()
        
        MESSAGE(("scipy_lbfgsb duration (s): %.1lf" % dtime))

    if SCIPY_CHOL:
        MESSAGE(("Now doing the inversion using scipy.linalg.solve "+datetime.datetime.now().strftime("%Y-%m-%d %H:%M")))
        from scipy.linalg import solve
        t0 = time()
        SLIP = solve(ATPA, ATPB, sym_pos=True, lower=False, overwrite_a=True, overwrite_b=False, debug=False, check_finite=False)
        dtime=time()-t0
        MESSAGE(("scipy_lsq_linear duration (s): %.1lf" % dtime))

    return(SLIP, dtime)


###############################################################################
# BVLS
###############################################################################

###############################################################################
def pyeq_bvls(A,B, bounds_slip,verbose=True):
###############################################################################

    """
       Stark & Parker algorithm
       Bounds are assumed to be between 0 and max_slip*m0, where max_slip is a vector of dim(m)
       This used to be a wrapper to use bvls.so f2py compiled from the original bvls.f routine. 
       Looks like scipy.optimize.lsq_linear does the bvls problem.
       

       G          :     Green's functions
       d          :     Observation vector
       max_slip   :     vector of max_slip
    """
    
    # import
    from scipy.optimize import lsq_linear    
    from time import time
    
    t0 = time()
    print("-- Inversion using BVLS and bounded values for slip")

    print("  -- Creating bounds vector for slip")
    lower_bounds=[]
    upper_bounds=[]
    for i in range(bounds_slip.shape[0]):
        lower_bounds.append(bounds_slip[i,0])
        upper_bounds.append(bounds_slip[i,1])
    
    print("  -- Running bvls optimization")

    res = lsq_linear(A, B.flatten(), bounds=(lower_bounds, upper_bounds), method='trf', tol=1e-5, lsq_solver=None, lsmr_tol=None, max_iter=None, verbose=0)
    
#    print "-- Results: ", msg, "\n"
    print("-- Elapsed time for inversion : %10.2f s" % (time()-t0))
#    print "-- Number of iterations       : %12d " % (loopa)
    
    
    
    SLIP = res.x
    norm = res.cost
    return(SLIP,norm)
