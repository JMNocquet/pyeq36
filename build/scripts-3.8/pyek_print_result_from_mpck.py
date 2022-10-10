#!/usr/local/geodesy/anaconda38/bin/python

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

import logging
import pyeq.message.message as MESSAGE
import pyeq.message.verbose_message as VERBOSE
import pyeq.message.error as ERROR
import pyeq.message.warning as WARNING
import pyeq.message.debug_message as DEBUG

###################################################################
# PARSE ARGUMENT LINE
###################################################################

prog_info="Print the results from a kinematic inversion"
prog_epilog="J.-M. Nocquet (Geoazur-UCA-IRD-CNRS-OCA) - December 2020"


parser = argparse.ArgumentParser(description=prog_info,epilog=prog_epilog,prefix_chars='-')
parser.add_argument('-mpck', action='store', type=str, dest='mpck',required=True,help='model pickle')
parser.add_argument('-odir', action='store', type=str, dest='odir',default=None,help='output directory')
parser.add_argument('-verbose', '-v', action='count',default=0,help='verbose mode')
parser.add_argument('-debug', action='count',default=0,help='debug mode')

args = parser.parse_args()

# verbose & debug
logging.getLogger("my_logger").setLevel(logging.WARNING)

if args.verbose>0:
    verbose=True
    logging.getLogger("my_logger").setLevel(logging.INFO)
else:
    verbose=False
if args.debug>0:
    debug=True
    logging.getLogger("my_logger").setLevel(logging.DEBUG)
else:debug=False


###################################################################
# LOAD MODEL PCK
###################################################################

if verbose:
    print("-- Loading %s (%.2f Gb) " % ( args.mpck , os.path.getsize( args.mpck ) /1024 / 1024 / 1024 ) )

model = pyeq.lib.objects.pyeq_model.load( args.mpck )

###################################################################
# CHANGE NAME & ODIR
###################################################################

if args.odir is not None:
    MESSAGE(("CHANGING MODEL NAME TO %s" % args.odir))
    model.odir = args.odir
    model.name = args.odir


###################################################################
# PRINT RESULTS
###################################################################

MESSAGE(("PRINTING RESULTS TO %s" % model.odir),level=1)

pyeq.log.print_results(model)
