#!/usr/bin/env python


###################################################################
# SCRIPT    : pyeq_magnitude.py
# AUTHOR    : JM NOCQUET
# DATE        : June 2012
# OUTPUT    :
# NOTE        :
#
###################################################################

###################################################################
# MODULES IMPORT
###################################################################

import sys
import argparse

import pyacs.lib.units

###################################################################
# PARSE ARGUMENT LINE
###################################################################

prog_info="Converts Mo and calculated magnitude from slip, width, length"
prog_epilog="J.-M. Nocquet (Geoazur-CNRS) - August 2012"


parser = argparse.ArgumentParser(description=prog_info,epilog=prog_epilog,prefix_chars='-')
parser.add_argument('-m0', action='append', type=float, dest='M0',help='M0 in N.m',default=[])
parser.add_argument('-m', action='append', type=float, dest='magnitude',help='magnitude',default=[])
parser.add_argument('-s', action='store', type=float, dest='slip',help='slip in meters',default=None)
parser.add_argument('-w', action='store', type=float, dest='width',help='Fault width on km',default=None)
parser.add_argument('-l', action='store', type=float, dest='length',help='Fault length on km',default=None)
parser.add_argument('-mu', action='store', type=float, dest='mu',help='rigidity module, default 3.1E10',default=3.1E10)
parser.add_argument('--verbose', '-v', action='count',default=0,help='verbose mode')
parser.add_argument('--debug', action='count',default=0,help='debug mode')

args = parser.parse_args()

if (len(sys.argv)<2):parser.print_help();sys.exit()

# verbose
if args.verbose>0:
    print("-- Verbose mode")
    verbose=True
else:
    verbose=False

###################################################################
# CASE M0
###################################################################

if args.M0 != []:
    sum_m0=0.0
    for M0 in args.M0:
        magnitude=pyacs.lib.units.moment_to_magnitude(M0)
        print( "-- moment %5.2E N.m            =>            magnitude %5.3lf" % (M0,magnitude))

        if len(args.M0)>1:
            sum_m0=sum_m0+M0
            print( "-- cumulative moment %5.2E N.m => cumulative magnitude %5.3lf" % (sum_m0,pyacs.lib.units.moment_to_magnitude(sum_m0)))
        
###################################################################
# CASE magnitude
###################################################################
        
if args.magnitude != []:
    sum_m0=0.0
    for magnitude in args.magnitude:
        m0=pyacs.lib.units.magnitude_to_moment(magnitude)
        print("-- magnitude %5.2lf             =>            moment %5.2E N.m " % (magnitude,m0))

        if len(args.magnitude)>1:
            sum_m0=sum_m0+m0
            print( "-- cumulative moment %5.2E N.m => cumulative magnitude %5.3lf" % (sum_m0,pyacs.lib.units.moment_to_magnitude(sum_m0)))
    
###################################################################
# CASE SLIP, AREA
###################################################################

if args.slip and args.width and args.length: 
    M0=args.mu*1.E6*args.length*args.width*args.slip
    magnitude=pyacs.lib.units.moment_to_magnitude(M0)
    print( "-- moment %5.2E N.m            =>            magnitude %5.3lf" % (M0,magnitude))
