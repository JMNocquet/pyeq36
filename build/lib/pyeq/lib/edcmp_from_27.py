"""
A few functions to help creating edcmp input and read results
"""
###############################################################################
# IMPORT
###############################################################################

import numpy as np


###############################################################################
def get_edcmp_comments():
    ###############################################################################

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

    fault_edcmp = ''
    fault_edcmp += '#         coord. origin: (40.739N, 30.05E)' + "\n"
    fault_edcmp += '#-------------------------------------------------------------------------------' + "\n"
    fault_edcmp += '# no  Slip   xs        ys       zs        length    width   strike   dip  rake' + "\n"
    fault_edcmp += '#-------------------------------------------------------------------------------' + "\n"
    fault_edcmp += '#   1   1.00 0.0d+0  50.0d+03  0.0d+03   64.0d+03  15.0d+03  45.0   45.0  -90.0' + "\n"

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

    end_edcmp = "#================================end of input===================================\n"

    return (intro_edcmp, obs_edcmp, output_edcmp, source_edcmp, fault_edcmp, earth_model_edcmp, end_edcmp)


###############################################################################
def make_edcmp_input_file(obs_x, obs_y, obs_z, fslip, fx, fy, fdepth, flength, fwidth, fstrike, fdip, frake, \
                          DISP=True, STRAIN=False, STRESS=False, TILT=False, elambda=0.32000E+11, emu=0.32000E+11):
    ###############################################################################

    # All input parameters are assumed to be 1-D array
    # obs_x, obs_y, obs_z, fx, fy, fdepth are assumed to be in local cartesian coordinates
    # with units km

    # check that obs_z is constant, otherwise it is a problem
    if isinstance(obs_z, np.ndarray):
        test_obs_z = np.zeros((obs_z.shape[0])) + obs_z[0]
        if not np.array_equal(obs_z, test_obs_z):
            import sys
            print
            '!!! ERROR in edcmp.make_edcmp_input_file: observation MUST be provided at the same constant depth'
            sys.exit()

    obs_x = np.array([obs_x]).flatten()
    obs_y = np.array([obs_y]).flatten()
    obs_z = np.array([obs_z]).flatten()

    fedcmp_input = 'input_edcmp.dat'
    fs = open(fedcmp_input, 'w')
    (intro_edcmp, obs_edcmp, output_edcmp, source_edcmp, fault_edcmp, earth_model_edcmp,
     end_edcmp) = get_edcmp_comments()
    # ---------------------------
    fs.write(intro_edcmp)
    # ---------------------------
    fs.write(obs_edcmp)
    # ---------------------------
    fs.write("  0\n")  # 0 indicates irregular positions option
    # ---------------------------
    fs.write("  %d\n" % obs_x.shape[0])  # number of observation points

    flag_full_line = 0
    for i in np.arange(obs_x.shape[0] - 1):
        fs.write("  (%5.4lfd+03,%5.4lfd+03)," % (obs_x[i], obs_y[i]))
        flag_full_line += 1
        if (flag_full_line == 3):
            flag_full_line = 0;
            fs.write("\n")
    i = obs_x.shape[0] - 1
    fs.write("  (%5.4lfd+03,%5.4lfd+03)" % (obs_x[i], obs_y[i]))
    fs.write("\n")
    # ---------------------------
    fs.write(output_edcmp)
    output_options = ("        %d               %d              %d              %d\n" % (
    int(DISP), int(STRAIN), int(STRESS), int(TILT)))

    fs.write(output_options)
    fs.write("  \'edcmp.disp\'    \'edcmp.strn\'   \'edcmp.strss\'  \'edcmp.tilt\'" + "\n")
    # ---------------------------
    fs.write(source_edcmp)
    fs.write(" %d\n" % fslip.shape[0])
    # ---------------------------
    fs.write(fault_edcmp)
    for i in np.arange(fslip.shape[0]):
        fs.write("  %d   %5.2lf %5.1lfd+03 %5.1lfd+03 %5.1lfd+03 %5.1lfd+03 %5.1lfd+03 %5.1lf %5.1lf %5.1lf\n" % \
                 (i + 1, fslip[i], fx[i], fy[i], fdepth[i], flength[i], fwidth[i], fstrike[i], fdip[i], frake[i]))
    # ---------------------------
    fs.write(earth_model_edcmp)
    # ---------------------------
    #    fs.write("  %8.3lfd+03  0.28758E+11  0.29353E+11\n" % rdepth)
    fs.write("  %8.3lfd+03  %13.5e %13.5e \n" % (obs_z[0], elambda, emu))
    fs.write(end_edcmp)

    fs.close()


###############################################################################
def to_1Dnparray(A):
    ###############################################################################
    """
    Returns a 1-D numpy array whatever is passed as argument
    """
    return (np.array(A.reshape(-1)))


###############################################################################
def run():
    ###############################################################################
    from pyacs.lib import syscmd
    syscmd.getstatusoutput('echo input_edcmp.dat > linput_edcmp')
    syscmd.getstatusoutput('edcmp < linput_edcmp')

#    import commands
#
#    commands.getstatusoutput('echo input_edcmp.dat > linput_edcmp')
#    commands.getstatusoutput('edcmp < linput_edcmp')


###############################################################################
def get_stress_from_edcmp_output():
    ###############################################################################
    """
    Converts edcmp stress file to numpy.array
    Format:
    #   X_m         Y_m         Sxx_Pa      Syy_Pa      Szz_Pa      Sxy_Pa      Syz_Pa      Szx_Pa
    """
    R = np.genfromtxt('edcmp.strss', dtype='float', comments='#')

    return (R)


###############################################################################
def get_disp_from_edcmp_output():
    ###############################################################################
    """
    Converts edcmp disp file to numpy.array
    Format:
    #   X_m         Y_m        Ux_m        Uy_m        Uz_m
    """
    U = np.genfromtxt('edcmp.disp', dtype='float', comments='#')

    return (U)


###############################################################################
def get_stress_tensor(R, i):
    ###############################################################################
    """
    Makes a 3x3 stres tensor from a np.array result from get_stress_from_edcmp_output
    """

    R = R.reshape(-1, 8)
    STRESS_TENSOR = np.zeros((3, 3))

    STRESS_TENSOR[0, 0] = R[i, 2]
    STRESS_TENSOR[1, 1] = R[i, 3]
    STRESS_TENSOR[2, 2] = R[i, 4]

    STRESS_TENSOR[0, 1] = R[i, 5]
    STRESS_TENSOR[1, 0] = R[i, 5]

    STRESS_TENSOR[1, 2] = R[i, 6]
    STRESS_TENSOR[2, 1] = R[i, 6]

    STRESS_TENSOR[0, 2] = R[i, 7]
    STRESS_TENSOR[2, 0] = R[i, 7]

    return (STRESS_TENSOR)