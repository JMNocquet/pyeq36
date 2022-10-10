#!/usr/local/geodesy/anaconda38/bin/python

###################################################################
# SCRIPT    : pyeq_coulomb.py
# AUTHOR    :
# DATE      : January 2017 - refactor to python3 June 2020
# INPUT     :
# OUTPUT    :
# NOTE      :
###################################################################

###################################################################
# MODULES IMPORT
###################################################################
import sys, os
from pyeq.lib import dislocation as DL

one_degree = 111.1

from pyeq.lib import Coulomb

import argparse
from argparse import RawTextHelpFormatter
import numpy as np

from time import time

t0 = time()

np.set_printoptions(4)

###################################################################
# PARSE ARGUMENT LINE
###################################################################

prog_info = "pyeq_coulomb.py : calculates coulomb stress\n"

prog_epilog = "J.-M. Nocquet (IRD-Geoazur & IGEPN-Ecuador) - January 2017 - refactoring June 2020"

# ----------------------------------------------------------
# INPUT
# ----------------------------------------------------------

parser = argparse.ArgumentParser(description=prog_info, epilog=prog_epilog, formatter_class=RawTextHelpFormatter)

# slip
parser.add_argument('--model', action='store', dest='model', type=str, default=None, help='slip model pyacs format')
parser.add_argument('--disloc', action='append', dest='ldisloc', type=str, default=[], \
                    help='dislocation parameters /fault_long/fault_lat/fault_depth/fault_length/fault_width/fault_strike/fault_dip/fault_rake/fault_slip')

# receivers faults
parser.add_argument('--receiver_fault', action='append', dest='lrcef', type=str, default=[], \
                    help='receiver fault (Coulomb stress will be calculated at its center) /fault_long/fault_lat/fault_depth/fault_strike/fault_dip/fault_rake')

parser.add_argument('--rcf_file', action='store', dest='rcf_file', type=str, default='', \
                    help='receiver faults file long,lat,depth,strike,dip,rake')

parser.add_argument('--rcf_geometry_file', action='store', dest='rcf_geometry_file', type=str, default='', \
                    help='receiver faults file as pyacs geometry file')

parser.add_argument('--rake', action='store', dest='rake', type=float, default=None, \
                    help='constant rake. Overwrite any previous values if provided')

parser.add_argument('-mu', action='store', dest='mu', type=float, required=True, \
                    help='friction coefficient')

parser.add_argument('--coor_type', action='store', dest='coor_type', type=str, default='geo', \
                    help='coordinates type : geo=geographic  or xyz=local cartesian')

# ----------------------------------------------------------
# general options
# ----------------------------------------------------------

parser.add_argument('--verbose', '-v', action='count', default=0, help='verbose mode')
parser.add_argument('--debug', action='count', default=0, help='debug mode')

parser.add_argument('-e', action='store', type=str, required=True, dest='experiment', help='experiment name')

args = parser.parse_args()

if (len(sys.argv) < 2): parser.print_help();sys.exit()

#############################################################################################################
# INITIALIZATION
#############################################################################################################
import os

# verbose
if args.verbose > 0:

    "-- Verbose mode"
    verbose = True
else:
    verbose = False
# debug
if args.debug > 0:

    "-- Debug mode"
    verbose = True
    debug = True
else:
    debug = False

#############################################################################################################
# DECIPHER SOURCE INPUT DATA
#############################################################################################################

NP_SOURCE = None

np_fault_source = None
np_record_source = None

# DISLOC OPTION
for disloc in args.ldisloc:
    if disloc[0] == '/': disloc = disloc[1:]
    [fault_long, fault_lat, fault_depth, fault_length, fault_width, fault_strike, fault_dip, fault_rake,
     fault_slip] = map(float, disloc.split('/'))

    # creates a dislocation object
    rdis_area = fault_length * fault_width
    max_slip = 1.
    disloc = DL.Dislocation(1, fault_long, fault_lat, fault_depth / one_degree, fault_strike, fault_dip,
                            fault_length / one_degree, fault_width / one_degree, rdis_area, fault_rake, max_slip)
    # get the corners
    (X1, X2, X3, X4) = disloc.corners()

    # append information

    SOURCE = np.array([[fault_long, fault_lat, fault_depth, fault_length, fault_width, fault_strike, fault_dip,
                        fault_rake, fault_slip]])
    fault_source = np.array([[X1[0], X1[1]], [X2[0], X2[1]], [X3[0], X3[1]], [X4[0], X4[1]], [X1[0], X1[1]]])
    record_source = np.array([[fault_depth, fault_slip, fault_rake]])

    if NP_SOURCE == None:
        NP_SOURCE = SOURCE
        np_fault_source = fault_source
        np_record_source = record_source
    else:
        NP_SOURCE = np.vstack([NP_SOURCE, SOURCE])
        np_fault_source = np.vstack([np_fault_source, fault_source])
        np_record_source = np.vstack([np_record_source, record_source])

if args.ldisloc != []:
    print('-- ', NP_SOURCE.shape[0], ' source dislocation read from option --disloc')

# MODEL OPTION

if args.model:
    try:
        SOL_PYEQ_MODEL = np.genfromtxt(args.model, comments='#')
    except:
        print('!!! Could not read ', args.model)
        import sys;
        sys.exit()

    max_slip = np.max(SOL_PYEQ_MODEL[:, -1])

    NP_SOURCE = np.zeros((SOL_PYEQ_MODEL.shape[0], 9))
    # long, lat, depth, length, width
    NP_SOURCE[:, :5] = SOL_PYEQ_MODEL[:, :5]
    # strike
    NP_SOURCE[:, 5] = SOL_PYEQ_MODEL[:, 7]
    # dip
    NP_SOURCE[:, 6] = SOL_PYEQ_MODEL[:, 8]
    # rake
    NP_SOURCE[:, 7] = SOL_PYEQ_MODEL[:, -5]
    # slip
    NP_SOURCE[:, 8] = SOL_PYEQ_MODEL[:, -1] / 1.E3

if args.model != '': print
'-- ', NP_SOURCE.shape[0], ' source dislocation read from ', args.model

print
"--------------------------------------------------------------------"

np.set_printoptions(precision=3, threshold=100000, linewidth=150)
if verbose:
    print("-------------------------------------------------------------------")
    print(" SOURCE DISLOCATION")
    print(" fault_long,fault_lat,fault_depth,fault_length,fault_width,fault_strike,fault_dip,fault_rake,fault_slip")
    print("-------------------------------------------------------------------")
    print(NP_SOURCE)
    print("-------------------------------------------------------------------")

#############################################################################################################
# DECIPHER RECEIVER INPUT DATA
#############################################################################################################

NP_RECEIVER = None  # long,lat,depth,strike,dip,rake

##### MANUAL ENTRIES /fault_long/fault_lat/fault_depth/fault_strike/fault_dip/fault_rake
for recf in args.lrcef:
    if recf[0] == '/': recf = recf[1:]
    [fault_long, fault_lat, fault_depth, fault_strike, fault_dip, fault_rake] = map(float, recf.split('/'))
    RECEIVER = np.array([[fault_long, fault_lat, fault_depth, fault_strike, fault_dip, fault_rake]])

    if NP_RECEIVER == None:
        NP_RECEIVER = RECEIVER
    else:
        NP_RECEIVER = np.vstack([NP_RECEIVER, RECEIVER])

if args.lrcef != []:
    print('-- ', NP_RECEIVER.shape[0], ' receiver points read from option --receiver_fault')

##### PYEQ_MODEL2 GEOMETRY FORMAT
##### !!! REQUIRES RAKE TO BE PROVIDED

if args.rcf_geometry_file:
    # args.rake is required

    if args.rake == None:
        print('!!! ERROR: rake option is required with option rcf_geometry_file ')
        import sys;

        sys.exit()

    try:
        GEOMETRY_PYEQ_MODEL = np.genfromtxt(args.rcf_geometry_file, comments='#')
    except:
        print('!!! Could not read ', args.rcf_geometry_file)
        import sys;

        sys.exit()

    NP_RECEIVER = np.zeros((GEOMETRY_PYEQ_MODEL.shape[0], 6))
    # fault_long,fault_lat,fault_depth
    NP_RECEIVER[:, :3] = GEOMETRY_PYEQ_MODEL[:, :3]
    # striken dip
    NP_RECEIVER[:, 3:4] = GEOMETRY_PYEQ_MODEL[:, 7:8]
    # rake
    NP_RECEIVER[:, -1] = args.rake

if args.rcf_geometry_file != '':
    print('-- ', NP_RECEIVER.shape[0], ' receiver points read from  geometry file ', args.rcf_geometry_file)

##### FREE FILE
if args.rcf_file:
    # format is assumed to be long,lat,depth,strike,dip,rake
    try:
        NP_RECEIVER = np.genfromtxt(args.rcf_file, comments='#')
    except:
        print('!!! Could not read ', args.rcf_file)
        import sys;

        sys.exit()

if args.rcf_file != '':
    print('-- ', NP_RECEIVER.shape[0], ' receiver points read from file ', args.rcf_file)

if verbose:
    print("-------------------------------------------------------------------")
    print(" RECEIVERS")
    print(" long  lat  depth strike dip rake")
    print("-------------------------------------------------------------------")
    print(NP_RECEIVER)
    print("-------------------------------------------------------------------")

#############################################################################################################
# DO THE ACTUAL CACULATION
#############################################################################################################


NP_COULOMB = np.zeros((NP_RECEIVER.shape[0]))
NP_SHEAR = np.zeros((NP_RECEIVER.shape[0]))
NP_NORMAL = np.zeros((NP_RECEIVER.shape[0]))

NP_REC_LON = NP_RECEIVER[:, 0]
NP_REC_LAT = NP_RECEIVER[:, 1]
# First, make sure that depth is positive downward
NP_REC_DEPTH = np.sqrt(NP_RECEIVER[:, 2] ** 2)
NP_REC_STRIKE = NP_RECEIVER[:, 3]
NP_REC_DIP = NP_RECEIVER[:, 4]
NP_REC_RAKE = NP_RECEIVER[:, 5]

NP_SOURCE_LON = NP_SOURCE[:, 0]
NP_SOURCE_LAT = NP_SOURCE[:, 1]
# First, make sure that depth is positive downward
NP_SOURCE_DEPTH = np.sqrt(NP_SOURCE[:, 2] ** 2)
NP_SOURCE_LENGTH = NP_SOURCE[:, 3]
NP_SOURCE_WIDTH = NP_SOURCE[:, 4]
NP_SOURCE_STRIKE = NP_SOURCE[:, 5]
NP_SOURCE_DIP = NP_SOURCE[:, 6]
NP_SOURCE_RAKE = NP_SOURCE[:, 7]
NP_SOURCE_SLIP = NP_SOURCE[:, 8]

header = "    long        lat      depth     strike        dip       rake     Coulomb (Pa)    Coulomb (bar)      Shear (bar)     Normal (Bar)"

print('-- Running the Coulomb calculation')

for i in np.arange(NP_RECEIVER.shape[0]):
    print('  -- receiver fault #', i, ' / ', NP_RECEIVER.shape[0])
    NP_COULOMB[i], NP_SHEAR[i], NP_NORMAL[i] = \
        Coulomb.coulomb_for_single_receiver_edcmp( \
            NP_REC_LON[i], NP_REC_LAT[i], NP_REC_DEPTH[i], NP_REC_STRIKE[i], NP_REC_DIP[i], NP_REC_RAKE[i], \
            NP_SOURCE_LON, NP_SOURCE_LAT, NP_SOURCE_DEPTH, NP_SOURCE_LENGTH, NP_SOURCE_WIDTH, NP_SOURCE_STRIKE,
            NP_SOURCE_DIP, NP_SOURCE_RAKE, NP_SOURCE_SLIP \
            , mu=args.mu, coor_type=args.coor_type)

print('-- End Coulomb calculation')

#############################################################################################################
# SAVE RESULTS
#############################################################################################################
print('-- Saving results')

# format long,lat,depth, strike,dip,rake,coulomb
COULOMB = np.zeros((NP_RECEIVER.shape[0], 10))

COULOMB[:, :6] = NP_RECEIVER[:, :6]
COULOMB[:, 6] = NP_COULOMB
COULOMB[:, 7] = NP_COULOMB * 1E-5
COULOMB[:, 8] = NP_SHEAR * 1E-5
COULOMB[:, 9] = NP_NORMAL * 1E-5

coulomb_file = shapefile = ("%s_coulomb.dat" % args.experiment)

print("-- Saving results to ", coulomb_file)
header = "    long        lat      depth     strike        dip       rake     Coulomb (Pa)    Coulomb (bar)      Shear (bar)     Normal (Bar)"
np.savetxt(coulomb_file, COULOMB, ("%10.4lf %10.4lf %10.2lf %10.2lf %10.2lf %10.2lf %16.2lf %16.4lf %16.4lf %16.4lf"),
           header=header)

# QGIS RESULTS



