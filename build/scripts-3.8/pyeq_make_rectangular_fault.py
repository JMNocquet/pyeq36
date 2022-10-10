#!/usr/local/geodesy/anaconda38/bin/python

###################################################################
# SCRIPT    : pyeq_model2_make_geometry_rectangular_fault.py
# AUTHOR    : J.-M. Nocquet IRD Geoazur & IGEPN Ecuador
# DATE      : 25/01/2017
###################################################################

###################################################################
# MODULES IMPORT
###################################################################

import sys
import argparse
from argparse import RawTextHelpFormatter
import numpy as np
from pyacs.lib.gmtpoint import GMT_Point
import pyacs.lib.faultslip


###################################################################
# PARSE ARGUMENT LINE
###################################################################

prog_info="pyeq_model2_make_geometry_rectangular_fault.py\n"
prog_info+="--------------------------------------------------------------------------------------------------------------------------------------\n"
prog_info+="Creates a geometry of rectangular fault(s) for slip inversion\n"

prog_epilog="J.-M. Nocquet (IRD-Geoazur & IGEPN Cuador) - January 2017"

#----------------------------------------------------------
#REQUIRED INPUT
#----------------------------------------------------------

parser = argparse.ArgumentParser(description=prog_info,epilog=prog_epilog,formatter_class=RawTextHelpFormatter)

parser.add_argument('-corner',  action='store', dest='corner',    type=str,required=True,help='/long/lat/depth of top left corner of the rectangular fault dec.deg and km (positive)')
parser.add_argument('-cornerf', action='store', dest='cornerf',   type=str, default='' , help='/long/lat of top right corner of the rectangular fault dec.deg and km (positive)')
parser.add_argument('-length',  action='store', dest='length',    type=float,help='length in km')
parser.add_argument('-width',   action='store', dest='width',     type=float,required=True,help='width in km')
parser.add_argument('-strike',  action='store', dest='strike',    type=float,help='strike in dec.deg')
parser.add_argument('-dip',     action='store' ,dest='dip',       type=float,required=True, help='dip  in dec.deg')
parser.add_argument('-e',       action='store', dest='experiment',type=str,required=True,help='experiment name')

#----------------------------------------------------------
# OPTIONAL INPUT
#----------------------------------------------------------

parser.add_argument('--nstrike', action='store' ,dest='nstrike',  default=1, type=int, help='number of subfaults division along strike')
parser.add_argument('--ndip', action='store' ,dest='ndip',  default=1, type=int, help='number of subfaults division along dip')
parser.add_argument('--qgis', action='count' , help='generate shapefile for visualization in QGIS')

#----------------------------------------------------------
# GENERAL INPT
#----------------------------------------------------------

parser.add_argument('--verbose', '-v', action='count', default=0, help='verbose mode')
parser.add_argument('--debug', action='count', default=0, help='debug mode')

args = parser.parse_args()

if (len(sys.argv)<2):parser.print_help();sys.exit()

# verbose
if args.verbose>0:
    print("-- Verbose mode")
    verbose=True
else:
    verbose=False
# debug
if args.debug>0:
    print("=> Debug mode")
    verbose=True
    debug=True
else:
    debug=False

#############################################################################################################
# INITIALIZATION
#############################################################################################################

fault_long,fault_lat,depth=list(map(float,args.corner.split('/')[1:]))

# case cornerf has been provided instead of strike and length

if args.cornerf !='':
        lon_f,lat_f = list(map(float,args.cornerf.split('/')[1:]))
        # strike
        args.strike   = pyacs.lib.faultslip.geo_to_strike( fault_long, fault_lat ,lon_f, lat_f)
    
        # length
        Mi = GMT_Point(code='Mi' , lon=fault_long , lat=fault_lat)
        Me = GMT_Point(code='Mf' , lon=lon_f , lat=lat_f)

        args.length   = Mi.spherical_distance( Me )
        # in km
        args.length   = args.length / 1.E3




l0=fault_long

if debug:print("-- corner coordinates ",fault_long,fault_lat,depth)

###################################################################
# CREATES DISLOCATION OBJECT AND SUBFAULTS IF REQUIRED
###################################################################

# Dislocations sources files
from pyeq.lib import eq_disloc_3d as dislocation
# one degree in km
one_degree=111.

dislocation=dislocation.Dislocation(index=0, x=fault_long, y=fault_lat, depth=depth, \
                                    strike=args.strike, dip=args.dip, \
                                    length=args.length , width=args.width, area = args.length * args.width )

dislocation.print_info()

if args.nstrike > 1 or args.ndip > 1 :
    if verbose:print("-- Dividing faults into ",args.nstrike, ' x ', args.ndip, ' subfaults')
    ldisloc=dislocation.subfaults(args.nstrike,args.ndip,coor_type='geo')
else:
    ldisloc=[dislocation]

    if verbose:
        for disloc in ldisloc:disloc.print_info()
    
###################################################################
# OPEN OUTPUT FILES
###################################################################

f_sources=args.experiment+'_sources_lp.dat'
f_dislocations=args.experiment+'_dislocations.dat'
f_gmt=args.experiment+'_rectangular_dislocations.gmt'
f_geometry=args.experiment+'_geometry.dat'
f_geometry_npy=args.experiment+'_geometry.npy'

n_dis = len( ldisloc )

print(("-- Geometry will include %04d dislocations" % n_dis))

Hdtype={'names':['    rdis_long','   rdis_lat','rdis_depth',\
                'rdis_length','rdis_width',' rdis_area','ratio_rdis_tdis',\
                '   strike','       dip',\
                'centroid_long','centroid_lat','centroid_depth',\
                'tdis_long1', ' tdis_lat1','tdis_depth1', \
                'tdis_long2',' tdis_lat2','tdis_depth2', \
                'tdis_long3',' tdis_lat3','tdis_depth3', \
                ' tdis_area'],\
         'formats':['f','f','f','f','f','f','f','f','f','f','f','f','f','f','f','f','f','f','f','f','f','f']}

GEOMETRY=np.zeros((n_dis,22))


for index_source in np.arange( n_dis ):

    disloc = ldisloc[index_source]
    
    ratio_rdis_tdis = 0
    
    X=disloc.centroid(coor_type='geo')
    
    centroid_long  = X[0]
    centroid_lat   = X[1]
    centroid_depth = X[2]
    
    tdis_long1 = tdis_lat1 = tdis_depth1 = 0.
    tdis_long2 = tdis_lat2 = tdis_depth2 = 0.
    tdis_long3 = tdis_lat3 = tdis_depth3 = 0.
    tdis_area = 0.
    
    GEOMETRY[index_source,:]=[disloc.x,disloc.y,disloc.depth,\
                              disloc.length,disloc.width,disloc.area, ratio_rdis_tdis,\
                              disloc.strike, disloc.dip,\
                              centroid_long,centroid_lat,centroid_depth,\
                              tdis_long1, tdis_lat1, tdis_depth1, \
                              tdis_long2, tdis_lat2, tdis_depth2, \
                              tdis_long3, tdis_lat3, tdis_depth3, \
                              tdis_area]


print("-- Writing NEW GEOMETRY FORMAT ",f_geometry,f_geometry_npy)

print("-- Sorting geometry array by descending latitudes")

GEOMETRY = GEOMETRY[GEOMETRY[:,10].argsort()[::-1]]

names=['rdis_long','rdis_lat','rdis_depth','rdis_length','rdis_width',\
       'rdis_area','ratio_rdis_tdis','strike','dip',\
       'centroid_long','centroid_lat','centroid_depth',\
       'tdis_long1','tdis_lat1','tdis_depth1',\
       'tdis_long2','tdis_lat2','tdis_depth2',\
       'tdis_long3','tdis_lat3','tdis_depth3','tdis_area']

np.save(f_geometry_npy, GEOMETRY)

GEOMETRY_TXT=np.zeros((GEOMETRY.shape[0],GEOMETRY.shape[1]+1))
GEOMETRY_TXT[:,0]=np.arange(GEOMETRY.shape[0])
GEOMETRY_TXT[:,1:] = GEOMETRY
header_cols=' '.join(Hdtype['names'])
format_header='%04d %10.5lf %10.5lf      %6.2lf \
     %6.2lf     %6.2lf    %6.2lf          %6.2lf\
    %6.2lf %10.2lf \
   %10.5lf   %10.5lf         %6.2lf \
%10.5lf %10.5lf      %6.2lf \
%10.5lf %10.5lf      %6.2lf \
%10.5lf %10.5lf      %6.2lf \
   %6.2lf '

np.savetxt(f_geometry, GEOMETRY_TXT, fmt=format_header, header=header_cols)

