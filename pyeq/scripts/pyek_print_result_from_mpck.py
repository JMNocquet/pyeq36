#!/usr/bin/env python

###################################################################
# SCRIPT    : pyek_print_result_from_mpck.py
# AUTHOR    : JM NOCQUET
# DATE      : December 2020
###################################################################

###################################################################
# MODULES IMPORT
###################################################################

import os
import argparse
import pyeq.lib.objects.pyeq_model
import pyeq.log.print_results

###################################################################
# PARSE ARGUMENT LINE
###################################################################

prog_info="Print the results from a kinematic inversion"
prog_epilog="J.-M. Nocquet (Geoazur-UCA-IRD-CNRS-OCA) - December 2020"


parser = argparse.ArgumentParser(description=prog_info,epilog=prog_epilog,prefix_chars='-')
parser.add_argument('-mpck', action='store', type=str, dest='mpck',required=True,help='model pickle')
parser.add_argument('-verbose', '-v', action='count',default=0,help='verbose mode')

args = parser.parse_args()

# verbose
if args.verbose>0:
    print("-- Verbose mode")
    verbose=True
else:
    verbose=False


###################################################################
# LOAD MODEL PCK
###################################################################

if verbose:
    print("-- Loading %s (%.2f Gb) " % ( args.mpck , os.path.getsize( args.mpck ) /1024 / 1024 / 1024 ) )

model = pyeq.lib.objects.pyeq_model.load( args.mpck )

###################################################################
# PRINT RESULTS
###################################################################

print("###############################################################################")
print("PRINTING RESULTS")
print("###############################################################################")

pyeq.log.print_results(model)
