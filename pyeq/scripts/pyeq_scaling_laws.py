#!/usr/bin/env python

###################################################################
# SCRIPT    : pyacs_eq_scaling_laws.py
# AUTHOR    : JM NOCQUET
# DATE      : June 2012
# OUTPUT    : 
# NOTE      :
#         
###################################################################

###################################################################
# MODULES IMPORT
###################################################################

import sys
import argparse
import math


###################################################################
# PARSE ARGUMENT LINE
###################################################################

prog_info="Provides prediction from published Earthquakes scaling laws"
prog_epilog="J.-M. Nocquet (Geoazur-CNRS) - June 2012"

list_eq_law={'WC94': 'Wells & Coppersmith (1994)'}

parser = argparse.ArgumentParser(description=prog_info,epilog=prog_epilog,prefix_chars='-')
parser.add_argument('--law', action='store', type=str, dest='law',help='Earthquake scaling laws; choose among ',choices=list_eq_law.keys(), default='WC94')
parser.add_argument('--magnitude', action='store', type=float, dest='magnitude',help='magnitude',default=6.)
parser.add_argument('--type', action='store', type=str, dest='type',help='type, choose among SS,N,R,ALL',default='ALL')
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

print("-- Earthquake scaling law used : ", args.law)
print("-- Fault type : ", args.type)

print("-- Magnitude : ", args.magnitude)
M=math.pow(10,3./2.*args.magnitude+9.1)
print ("-- Moment Magnitude %5.2E [N.m]" % M)

###################################################################
# WC94
###################################################################
if args.law == 'WC94':
    
    # Maximum Displacement
    # log(MD)= a + b * M
    A={'SS':-7.03,'R':-1.84,'N':-5.90,'ALL':-5.46}
    B={'SS':1.03,'R':0.29,'N':0.89,'ALL':0.82}
    log_MD= A[args.type] + B[args.type]*args.magnitude
    MD=math.pow(10,log_MD)
    print("-- Maximum displacement (m): %5.1lf " % MD)
    
    # Avergae Displacement
    # log(AD)= a + b * M

    A={'SS':-6.32,'R':-0.74,'N':-4.45,'ALL':-4.80}
    B={'SS':0.90,'R':0.08,'N':0.63,'ALL':0.69}
    log_AD= A[args.type] + B[args.type]*args.magnitude
    AD=math.pow(10,log_AD)
    print("-- Average displacement (m): %5.3lf " % AD )
    
    # RLD Subsurface rupture length
    # log(RLD) = a + b * M
    A={'SS':-2.57,'R':-2.42,'N':-1.88,'ALL':-2.44}
    B={'SS':0.62,'R':0.58,'N':0.50,'ALL':0.59}

    log_RLD= A[args.type] + B[args.type]*args.magnitude
    RLD=math.pow(10,log_RLD)
    print("-- Subsurface rupture length (km): %5.1lf " % RLD)
    
    # RW downdip rupture width
    # log(RW) = a + b * M
    A={'SS':-0.76,'R':-1.61,'N':-1.14,'ALL':-1.01}
    B={'SS':0.27,'R':0.41,'N':0.35,'ALL':0.32}

    log_RW= A[args.type] + B[args.type]*args.magnitude
    RW=math.pow(10,log_RW)
    print("-- Downdip rupture width (km): %5.1lf " % RW )
    
    # Check consistency
    
    mu=3.E10 # shear modulus in Pa
    
    M0_law= mu * AD * RLD * RW * 1.E6
    print ("-- Retrieved Moment Magnitude %5.2E [N.m]" % M)
    magnitude_law= 2./3.*(math.log10(M0_law)-9.1)
    print("-- Retrieved moment magnitude %8.1f " % magnitude_law )
    
    
    
    
