"""
Python 3.6 wrapper to edcmp
"""


###############################################################################
def disp_tilt_strain_stress_from_edcmp(array_faults, array_gps, verbose=False):
###############################################################################
    """Calculates the 3D displacement at an array of points according to fault parameters

    :param array_faults : 2D numpy array of faults parameters [slip,xf,yf,depth,length,width,strike,dip,rake]
    :param array_gps: 2D numpy array with x y

    :note: xf & yf are E, N, U (pyacs convention). The program change it to edcmp convention (x=N,y=E and z positive downward)
    """

    # import
    import numpy as np

# input edcmp file
    fedcmp_input = 'input_edcmp.dat'
    fs = open(fedcmp_input, 'w')

    # comment intro
    intro_edcmp = ''
    intro_edcmp += '#===============================================================================' + "\n"
    intro_edcmp += '# This is the input file of FORTRAN77 program "edcomp" for calculating' + "\n"
    intro_edcmp += '# earthquakes static deformations (3 displacement components, 6 strain/stress' + "\n"
    intro_edcmp += '# components and 2 vertical tilt components) based on the dislocation theory.' + "\n"
    intro_edcmp += '# The earth model used is either a homogeneous or multi-layered, isotropic and' + "\n"
    intro_edcmp += '# elastic half-space. The earthquake source is represented by an arbitrary number' + "\n"
    intro_edcmp += '# of rectangular dislocation planes.' + "\n"
    intro_edcmp += '#' + "\n"
    intro_edcmp += '# Note that the Cartesian coordinate system is used and the seismological' + "\n"
    intro_edcmp += '# convention is adopted, that is, x is northward, y is eastward, and z is' + "\n"
    intro_edcmp += '# downward.' + "\n"
    intro_edcmp += '#' + "\n"
    intro_edcmp += '# First implemented in Potsdam, Feb, 1999' + "\n"
    intro_edcmp += '# Last modified: Potsdam, Nov, 2001' + "\n"
    intro_edcmp += '#' + "\n"
    intro_edcmp += '# by' + "\n"
    intro_edcmp += '# Rongjiang Wang, Frank Roth, & Francisco Lorenzo' + "\n"
    intro_edcmp += '# GeoForschungsZetrum Potsdam, Telegrafenberg, 14473 Potsdam, Germany' + "\n"
    intro_edcmp += '#' + "\n"
    intro_edcmp += '# For questions and suggestions please send e-mails to wang@gfz-potsdam.de' + "\n"
    intro_edcmp += '#===============================================================================' + "\n"

    # comment obs
    obs_edcmp = ''
    obs_edcmp += '# OBSERVATION ARRAY' + "\n"
    obs_edcmp += '# =================' + "\n"
    obs_edcmp += '# 1. switch for irregular positions (0) or a 1D profile (1)' + "\n"
    obs_edcmp += '#    or a rectangular 2D observation array (2): ixyr' + "\n"
    obs_edcmp += '#' + "\n"
    obs_edcmp += '#    IF (1 for irregular observation positions) THEN' + "\n"
    obs_edcmp += '#    ' + "\n"
    obs_edcmp += '# 2. number of positions: nr' + "\n"
    obs_edcmp += '# 3. coordinates of the observations: (xr(i),yr(i)),i=1,nr' + "\n"
    obs_edcmp += '#' + "\n"
    obs_edcmp += '#    ELSE IF (the switch = 1 for a 1D profile) THEN' + "\n"
    obs_edcmp += '#' + "\n"
    obs_edcmp += '# 2. number of position samples of the profile: nr' + "\n"
    obs_edcmp += '# 3. the start and end positions: (xr1,yr1), (xr2,yr2)' + "\n"
    obs_edcmp += '#' + "\n"
    obs_edcmp += '#    ELSE IF (2 for rectanglular 2D observation array) THEN' + "\n"
    obs_edcmp += '#' + "\n"
    obs_edcmp += '# 2. number of xr samples, start and end values [m]: nxr, xr1,xr2' + "\n"
    obs_edcmp += '# 3. number of yr samples, start and end values [m]: nyr, yr1,yr2' + "\n"
    obs_edcmp += '#' + "\n"
    obs_edcmp += '#    Note that the total number of observation positions (nr or nxr*nyr)' + "\n"
    obs_edcmp += '#    should be <= NRECMAX (see edcglobal.h)!' + "\n"
    obs_edcmp += '#===============================================================================' + "\n"
    obs_edcmp += '#' + "\n"

    # comment output
    output_edcmp = ''
    output_edcmp += '#===============================================================================' + "\n"
    output_edcmp += '# OUTPUTS' + "\n"
    output_edcmp += '# =======' + "\n"
    output_edcmp += '# 1. output directory in char format: outdir' + "\n"
    output_edcmp += '# 2. select the desired outputs (1/0 = yes/no)' + "\n"
    output_edcmp += '# 3. the file names in char format for displacement vector, strain tensor,' + "\n"
    output_edcmp += '#    stress tensor, and vertical tilts:' + "\n"
    output_edcmp += '#    dispfile, strainfile, stressfile, tiltfile' + "\n"
    output_edcmp += '#' + "\n"
    output_edcmp += '#    Note that all file or directory names should not be longer than 80' + "\n"
    output_edcmp += '#    characters. Directories must be ended by / (unix) or \ (dos)!' + "\n"
    output_edcmp += '#===============================================================================' + "\n"
    output_edcmp += '  \'./\'' + "\n"
    output_edcmp += '        1               1              1              1' + "\n"

    # comment source
    source_edcmp = ''

    source_edcmp += '#===============================================================================' + "\n"
    source_edcmp += '# RECTANGLAR DISLOCATION SOURCES' + "\n"
    source_edcmp += '# ==============================' + "\n"
    source_edcmp += '# 1. number of the source rectangles: ns (<= NSMAX in edcglobal.h)' + "\n"
    source_edcmp += '# 2. the 6 parameters for the 1. source rectangle:' + "\n"
    source_edcmp += '#    Slip [m],' + "\n"
    source_edcmp += '#    coordinates of the upper reference point for strike (xs, ys, zs) [m],' + "\n"
    source_edcmp += '#    length (strike direction) [m], and width (dip direction) [m],' + "\n"
    source_edcmp += '#    strike [deg], dip [deg], and rake [deg]' + "\n"
    source_edcmp += '# 3. ... for the 2. source ...' + "\n"
    source_edcmp += '# ...' + "\n"
    source_edcmp += '#                   N' + "\n"
    source_edcmp += '#                  /' + "\n"
    source_edcmp += '#                 /| strike' + "\n"
    source_edcmp += '#         Ref:-> @------------------------' + "\n"
    source_edcmp += '#                |\        p .            \ W' + "\n"
    source_edcmp += '#                :-\      i .              \ i' + "\n"
    source_edcmp += '#                |  \    l .                \ d' + "\n"
    source_edcmp += '#                :90 \  S .                  \ t' + "\n"
    source_edcmp += '#                |-dip\  .                    \ h' + "\n"
    source_edcmp += '#                :     \. | rake               \ ' + "\n"
    source_edcmp += '#                Z      -------------------------' + "\n"
    source_edcmp += '#                              L e n g t h' + "\n"
    source_edcmp += '#' + "\n"
    source_edcmp += '#    Note that if one of the parameters length and width = 0, then a line source' + "\n"
    source_edcmp += '#    will be considered and the displocation parameter Slip has the unit m^2 if' + "\n"
    source_edcmp += '#    both length and width = 0, then a point source will be considered and the' + "\n"
    source_edcmp += '#    Slip has the unit m^3.' + "\n"
    source_edcmp += '#===============================================================================' + "\n"

    # comment fault
    fault_edcmp = ''
    fault_edcmp += '#         coord. origin: (40.739N, 30.05E)' + "\n"
    fault_edcmp += '#-------------------------------------------------------------------------------' + "\n"
    fault_edcmp += '# no  Slip   xs        ys       zs        length    width   strike   dip  rake' + "\n"
    fault_edcmp += '#-------------------------------------------------------------------------------' + "\n"
    fault_edcmp += '#   1   1.00 0.0d+0  50.0d+03  0.0d+03   64.0d+03  15.0d+03  45.0   45.0  -90.0' + "\n"

    # comment earth_model
    earth_model_edcmp = ''
    earth_model_edcmp += '#===============================================================================' + "\n"
    earth_model_edcmp += '# If the earth model used is a layered half-space, then the numerical Green\'s' + "\n"
    earth_model_edcmp += '# function approach is applied. The Green\'s functions should have been prepared' + "\n"
    earth_model_edcmp += '# with the program "edgrn" before the program "edcmp" is started. In this case,' + "\n"
    earth_model_edcmp += '# the following input data give the addresses where the Green\'s functions have' + "\n"
    earth_model_edcmp += '# been stored and the grid side to be used for the automatic discretization' + "\n"
    earth_model_edcmp += '# of the finite rectangular sources.' + "\n"
    earth_model_edcmp += '#' + "\n"
    earth_model_edcmp += '# If the earth model used is a homogeneous half-space, then the analytical' + "\n"
    earth_model_edcmp += '# method of Okada (1992) is applied. In this case, the Green\'s functions are' + "\n"
    earth_model_edcmp += '# not needed, and the following input data give the shear modulus and the' + "\n"
    earth_model_edcmp += '# Poisson ratio of the model.' + "\n"
    earth_model_edcmp += '#===============================================================================' + "\n"
    earth_model_edcmp += '# CHOICE OF EARTH MODEL' + "\n"
    earth_model_edcmp += '# =====================' + "\n"
    earth_model_edcmp += '# 1. switch for layered (1) or homogeneous (0) model' + "\n"
    earth_model_edcmp += '#' + "\n"
    earth_model_edcmp += '#    IF (layered model) THEN' + "\n"
    earth_model_edcmp += '#' + "\n"
    earth_model_edcmp += '# 2. directory of the Green\'s functions and the three files for the' + "\n"
    earth_model_edcmp += '#    fundamental Green\'s functions: grndir, grnfiles(3)' + "\n"
    earth_model_edcmp += '#' + "\n"
    earth_model_edcmp += '#    Note that all file or directory names should not be longer than 80' + "\n"
    earth_model_edcmp += '#    characters. Directories must be ended by / (unix) or \ (dos)!' + "\n"
    earth_model_edcmp += '#' + "\n"
    earth_model_edcmp += '#    ELSE (homogeneous model) THEN' + "\n"
    earth_model_edcmp += '#' + "\n"
    earth_model_edcmp += '# 2. the observation depth, the two Lame constants parameters of the homogeneous' + "\n"
    earth_model_edcmp += '#    model: zrec [m], lambda [Pa], mu [Pa]' + "\n"
    earth_model_edcmp += '#===============================================================================' + "\n"
    earth_model_edcmp += '#  1' + "\n"
    earth_model_edcmp += '#  \'./grnfcts/\'  \'izmhs.ss\'  \'izmhs.ds\'  \'izmhs.cl\'' + "\n"
    earth_model_edcmp += '  0' + "\n"
#    earth_model_edcmp += '  0.00d+00  0.28758E+11  0.29353E+11' + "\n"
    earth_model_edcmp += '  0.00d+00  0.28758E+11  0.29353E+11' + "\n"
    earth_model_edcmp += '#================================end of input===================================' + "\n"

    # begin filling file

    fs.write(intro_edcmp)
    fs.write(obs_edcmp)
    fs.write("  0\n")

    # observation file
    # edcmp convention X = north, Y = East
    fs.write("  %d\n" % array_gps.shape[0])
    line = 0
    # change on July 2019
    # xf0=yf0=0 and gps locations are array_gps[i,1] - yf0 ,array_gps[i,0] -xf0


    xf0 = array_faults[0, 1]
    yf0 = array_faults[0, 2]

    if verbose:
        print("-- reference coordinates: xf0 yf0 %lf %lf " % (xf0, yf0))

    for i in range(array_gps.shape[0] - 1):
        #fs.write("  (%5.4lfd+03,%5.4lfd+03)," % (array_gps[i, 1] - yf0, array_gps[i, 0] - xf0))
        fs.write("  (%5.4lf,%5.4lf)," % (array_gps[i, 1] - yf0, array_gps[i, 0] - xf0))
        line += 1
        if (line == 3):
            line = 0;
            fs.write("\n")
    i = array_gps.shape[0] - 1
    #fs.write("  (%5.4lfd+03,%5.4lfd+03)" % (array_gps[i, 1] - yf0, array_gps[i, 0] - xf0))
    fs.write("  (%5.4lf,%5.4lf)" % (array_gps[i, 1] - yf0, array_gps[i, 0] - xf0))
    fs.write("\n")
    fs.write(output_edcmp)
    fs.write("  \'edcmp.disp\'    \'edcmp.strn\'   \'edcmp.strss\'  \'edcmp.tilt\'" + "\n")
    fs.write(source_edcmp)
    # number of faults
    nf = array_faults.shape[0]
    fs.write(" %d\n" % nf)
    fs.write(fault_edcmp)

    for idxf in np.arange(nf):
        [slip, xf, yf, depth, length, width, strike, dip, rake] = array_faults[idxf]

        # change on July 2019
        # xf=yf=0 and gps locations are array_gps[i,1] - yf ,array_gps[i,0] -xf
        #        fs.write("  %d   %5.2lf %5.1lfd+03 %5.1lfd+03 %5.1lfd+03 %5.1lfd+03 %5.1lfd+03 %5.1lf %5.1lf %5.1lf\n" % \
        #                 (idxf+1,slip,yf,xf,depth,length,width,strike,dip,rake))
        # change 14/12/2020 - changed sign of depth to account for PYACS to EDCMP convention
#        fs.write("  %d   %5.2lf %5.1lfd+03 %5.1lfd+03 %5.1lfd+03 %5.1lfd+03 %5.1lfd+03 %5.1lf %5.1lf %5.1lf\n" % \
#                 (idxf + 1, slip, yf - yf0, xf - xf0, -depth, length, width, strike, dip, rake))
        fs.write("  %d   %5.2lf %5.1lf %5.1lf %5.1lf %5.1lf %5.1lf %5.1lf %5.1lf %5.1lf\n" % \
                 (idxf + 1, slip, yf - yf0, xf - xf0, -depth, length, width, strike, dip, rake))

    # earth model
    fs.write(earth_model_edcmp)
    fs.close()

    # now runs edcmp
    fs = open('linput_edcmp', 'w')
    fs.write("%s\n" % fedcmp_input)
    fs.close()

    from pyacs.lib import syscmd
    syscmd.getstatusoutput('edcmp < linput_edcmp')

    # reads disp results
    R = np.genfromtxt('edcmp.disp', dtype='float', comments='#')
    #### if edcmp.disp only only one record a 1-D array (vector) is returned to R instead of a 2-D array
    if (R.ndim == 1): R.resize(1, 5)
    DISP = np.zeros(R.shape)
    DISP[:, 0] = R[:, 1] * 1.E-3
    DISP[:, 1] = R[:, 0] * 1.E-3
    DISP[:, 2] = R[:, 3]
    DISP[:, 3] = R[:, 2]
    DISP[:, 4] = -R[:, 4]

    # read strain results
    # edcmp convention: #   X_m         Y_m         Exx         Eyy         Ezz         Exy         Eyz         Ezx with X=north, Y=East, Z=downward
    # pyacs convention: #   X_m         Y_m         Exx         Eyy         Ezz         Exy         Exz         Eyz with X=East, Y=North, Z=upward
    R = np.genfromtxt('edcmp.strn', dtype='float', comments='#')
    #### if edcmp.disp only only one record a 1-D array (vector) is returned to R instead of a 2-D array
    if (R.ndim == 1): R.resize(1, 8)
    STRAIN = np.zeros(R.shape)
    STRAIN[:, 0] = R[:, 1] * 1.E-3 # y->x
    STRAIN[:, 1] = R[:, 0] * 1.E-3 # x->y
    STRAIN[:, 2] = R[:, 3] # Exx # yy -> xx
    STRAIN[:, 3] = R[:, 2] # Eyy # xx-> yy
    STRAIN[:, 4] = R[:, 4] # Ezz # zz-> zz
    STRAIN[:, 5] = R[:, 5] # Exy # yx -> xy
    STRAIN[:, 6] = R[:, 6] # Exz # yz -> xz
    STRAIN[:, 7] = R[:, 7] # Eyz # xz -> yz


    # read stress results
    # edcmp convention: #   X_m         Y_m         Exx         Eyy         Ezz         Exy         Eyz         Ezx with X=north, Y=East, Z=downward
    # pyacs convention: #   X_m         Y_m         Exx         Eyy         Ezz         Exy         Exz         Eyz with X=East, Y=North, Z=upward
    R = np.genfromtxt('edcmp.strss', dtype='float', comments='#')
    #### if edcmp.disp only only one record a 1-D array (vector) is returned to R instead of a 2-D array
    if (R.ndim == 1): R.resize(1, 8)
    STRESS = np.zeros(R.shape)
    STRESS[:, 0] = R[:, 1] * 1.E-3  # y->x
    STRESS[:, 1] = R[:, 0] * 1.E-3  # x->y
    STRESS[:, 2] = R[:, 3]  # Exx # yy -> xx
    STRESS[:, 3] = R[:, 2]  # Eyy # xx-> yy
    STRESS[:, 4] = R[:, 4]  # Ezz # zz-> zz
    STRESS[:, 5] = R[:, 5]  # Exy # yx -> xy
    STRESS[:, 6] = R[:, 6]  # Exz # yz -> xz
    STRESS[:, 7] = R[:, 7]  # Eyz # xz -> yz


    # read tilt results
    R = np.genfromtxt('edcmp.tilt', dtype='float', comments='#')
    #### if edcmp.disp only only one record a 1-D array (vector) is returned to R instead of a 2-D array
    if (R.ndim == 1): R.resize(1, 4)
    TILT = np.zeros(R.shape)
    TILT[:, 0] = R[:, 1] * 1.E-3  # y->x
    TILT[:, 1] = R[:, 0] * 1.E-3  # x->y
    TILT[:, 2] = R[:, 3]  # Exx # yy -> xx
    TILT[:, 3] = R[:, 2]  # Eyy # xx-> yy

    return ( DISP, STRAIN, STRESS, TILT )
