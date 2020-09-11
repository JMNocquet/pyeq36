#!/usr/bin/env python

###################################################################
# SCRIPT    : pyeq_static_model_to_qgis.py
# AUTHOR    : JM NOCQUET
# DATE      : September 2012 - Feruary 2018
###################################################################


import os

def file_base_name(file_name):
    if '.' in file_name:
        separator_index = file_name.index('.')
        base_name = file_name[:separator_index]
        return base_name
    else:
        return file_name

def path_base_name(path):
    file_name = os.path.basename(path)
    return file_base_name(file_name)



###################################################################
# MODULES IMPORT
###################################################################

import sys
import argparse
import pyacs.lib.shapefile

###################################################################
# PARSE ARGUMENT LINE
###################################################################

parser = argparse.ArgumentParser()
parser.add_argument('-dat', action='store', type=str, dest='dat',required=True,help='solution (.dat) output from pyeq_model_*inversion.py')
parser.add_argument('-shapefile', action='store', type=str, dest='shapefile',default=None,help='output shapefile name')
parser.add_argument('--t', action='count', default=0, help='output triangular dislocation (default is rectangular)')
parser.add_argument('--verbose', action='count', default=0, help='verbose mode')

if (len(sys.argv)<2):parser.print_help()
args = parser.parse_args()

if args.shapefile is None:
    shp = path_base_name( args.dat )

if args.t>0:
    dis_type='tde'
else:
    dis_type='rec'

if args.verbose>0:
    verbose=True
else:
    verbose=False

pyacs.lib.shapefile.static_slip_to_shapefile(args.dat, shp, dis_type, verbose)