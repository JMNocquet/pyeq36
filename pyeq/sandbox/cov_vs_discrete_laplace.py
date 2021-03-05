"""
This sandbox package shows tests to check the equivalence of covariance matrix vs discrete laplacian approach
for smoothing a model.

"""

def make_rectangular_grid( x0,y0, nx, ny, stepx, stepy ):
    """
    defines a rectangular grid
    """

    import numpy as np

    x = np.arange(x0, nx*stepx , stepx)
    y = np.arange(y0, ny*stepy , stepy)

    # creates the distance matrix
    from scipy.spatial import distance_matrix
    mesh = np.meshgrid(x, y)
    v = np.vstack((mesh[0].flatten(), mesh[1].flatten())).T
    Dm = distance_matrix(v, v, p=2)

    return Dm

def make_inv_corr( Dm ):

    import numpy as np
    corr =  np.exp(-Dm)
    icorr = np.linalg.inv(corr)

    return icorr

def make_make_inv_sqrt( Dm ):
    import numpy as np
    import scipy.linalg.lapack
    corr =  np.exp(-Dm)

    L,_ = scipy.linalg.lapack.dpotrf( corr )

    return( scipy.linalg.lapack.dpotrf( L )[0] )


def make_discrete_laplace_matrix( Dm , n=5 ):
    """
    n should be 4 or 9
    """

    import numpy as np

    L = np.copy(Dm)

    if n ==  9:
        # laplacian Operator
        # 9 points stencil
        L[L > 1.5] = 0
        L[L == 1] = 0.5
        L[L > 1] = 0.25
        np.fill_diagonal(L, -3)

    if n == 5:
        # 5 points stencil
        L[L > 1.] = 0
        L[L == 1] = 1
        np.fill_diagonal(L, -4)

    return L

def normalize( C ):
    import numpy as np
    return C / np.median( np.diag(C) )

def corr_from_discrete_laplace( L ):
    """
    computes a correlation matrix from laplace operator
    """

    import numpy as np

    L = -L
    np.fill_diagonal(L,np.diagonal(L)+1)
    CL = np.dot(L.T,L)

    return(normalize(CL) )

def dlo_ec_ig( Dm , dc ):
    """
    Computes the discrete Laplace operator from a correlation matrix with
    decreasing exponential kernel on an irregular grid.

    :param Dm: distance matrix
    :param dc: critical distance for the decreasing exponential correlation matrix
    :return DLO: the discrete Laplace operator as a 2D numpy array
    """

    # import
    import numpy as np
    import scipy.linalg

    # inverse correlation matrix
    icorr = make_inv_corr( Dm )
    # compute the inverse
    #import scipy.linalg
    #icoor = scipy.linalg.inv(icorr)

    # compute the eigen values and vectors
    e,v = scipy.linalg.eigh(icorr, eigvals_only=False)

    # check the results
    print(e)
    print(v)

    sqrt_e = np.sqrt(e)
    sqrt_icorr = np.dot( np.dot( v, np.diag(sqrt_e) ) , v.T )

    print('check')
    print(np.allclose(icorr,np.dot(sqrt_icorr,sqrt_icorr)))

    # computes the sum of non-diagonal elements
    diag = np.diagonal( sqrt_icorr )
    S = np.sum(sqrt_icorr,axis=1)

    alpha = S

    #alpha = np.diag(sqrt_icorr) / ( 1+S)
    print("alpha")
    print(alpha)
    best_alpha = np.median(alpha)
    print("best alpha: %.3lf " % (best_alpha))

    DL = np.eye(Dm.shape[0]) - 1./best_alpha * sqrt_icorr

    return DL



def plot_01( CL5, CL9, icorr ):
    import matplotlib.pyplot as plt

    plt.figure(figsize=(10,7))

    # plot 1,1
    plt.subplot(231)
    plt.imshow(CL5,cmap='seismic',vmin=-1.2,vmax=1.2)
    plt.colorbar(fraction=0.046, pad=0.04)
    plt.title('Discrete Laplace n=5 (DL5)')


    # plot 1,2
    plt.subplot(232)
    plt.imshow(CL9,cmap='seismic',vmin=-1.2,vmax=1.2)
    plt.colorbar(fraction=0.046, pad=0.04)
    plt.title('Discrete Laplace n=9 (DL9)')

    # plot 1,3
    plt.subplot(233)
    plt.imshow(icorr,cmap='seismic',vmin=-1.2,vmax=1.2)
    plt.colorbar(fraction=0.046, pad=0.04)
    plt.title('Inverse Exponential')


    # plot 2,1
    plt.subplot(234)
    plt.imshow(CL5-CL9,cmap='seismic',vmin=-1.2,vmax=1.2)
    plt.colorbar(fraction=0.046, pad=0.04)
    plt.title('DL5 - DL9')

    # plot 2,2
    plt.subplot(235)
    plt.imshow(CL5-icorr,cmap='seismic',vmin=-1.2,vmax=1.2)
    plt.colorbar(fraction=0.046, pad=0.04)
    plt.title('DL5 - Inverse Exponential')

    # plot 2,2
    plt.subplot(236)
    plt.imshow(CL9-icorr,cmap='seismic',vmin=-1.2,vmax=1.2)
    plt.colorbar(fraction=0.046, pad=0.04)
    plt.title('DL9 - Inverse Exponential')

    plt.tight_layout()
    plt.show()

def make_example_01():

    x0 = 0
    nx = 10
    stepx = 1

    y0 = 0
    ny = 5
    stepy =1

    # make distance matrix
    Dm = make_rectangular_grid( x0,y0, nx, ny, stepx, stepy )
    # inverse correlation matrix
    icorr = normalize( make_inv_corr( Dm ) )
    # discrete laplace
    DL5 = make_discrete_laplace_matrix(Dm, n=5)
    DL9 = make_discrete_laplace_matrix(Dm, n=9)
    # corr_from_discrete_laplace
    CL5 = corr_from_discrete_laplace(DL5)
    CL9 = corr_from_discrete_laplace(DL9)
    plot_01(CL5, CL9, icorr)

def make_example_02():

    import numpy as np

    x0 = 0
    nx = 20
    stepx = 1

    y0 = 0
    ny = 10
    stepy =1

    # make distance matrix
    Dm = make_rectangular_grid( x0,y0, nx, ny, stepx, stepy )
    # inverse correlation matrix
    icorr = make_inv_corr( Dm )
    # compute the choleski decomposition
    import scipy.linalg
    L = scipy.linalg.cholesky(icorr)
    # compute the laplacian
    DL = dlo_ec_ig( Dm , 1 )
    DL5 = make_discrete_laplace_matrix(Dm, n=5)
    DL9 = make_discrete_laplace_matrix(Dm, n=9)

    plot_01( DL5,DL9,DL)


    nDL5 = normalize(np.dot(DL5.T,DL5))
    nDL9 = normalize(np.dot(DL9.T,DL9))
    nDL  = normalize( np.dot(DL.T,DL) )

    print( np.max(nDL9-nDL))
    print( np.min(nDL9-nDL))

    plot_01( nDL5,nDL9,nDL)


