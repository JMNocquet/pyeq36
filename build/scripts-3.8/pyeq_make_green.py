#!/usr/local/geodesy/anaconda38/bin/python

###################################################################
# SCRIPT    : pyeq_make_green.py
# AUTHOR    : 
# DATE      : December 2011 - Update Fall 2017 - December 2020
# OUTPUT    : 
# NOTE      :
#         
###################################################################


###################################################################
# PYTHON MODULES IMPORT
###################################################################

import argparse, sys
import numpy as np
from time import time
from colors import red
from progress.bar import Bar
import logging

###################################################################
# PYACS ADDITIONAL MODULES IMPORT
###################################################################

import pyacs.lib.coordinates
import pyeq.lib.geometry.to_np_array
import pyeq.message.message as MESSAGE
import pyeq.green.make
import pyeq.message.verbose_message as VERBOSE
import pyeq.message.error as ERROR
import pyeq.message.warning as WARNING
import pyeq.message.debug_message as DEBUG
#from pyacs.lib import gmtpoint as GMT_Point
#from pyeq.lib import eq_disloc_3d as Dislocation
#from pyeq.lib import lib_inversion

###################################################################
# PARSE ARGUMENT LINE
###################################################################

prog_info="pyeq_make_green.py : calculates the transfert (Green) functions relating surface deformation observations to slip on a list of individual faults\n"
prog_info+="                    reads either a pyeq geometry dat or npy file\n"
prog_info+="                    format for geometry is:\n"
prog_info+="00: dislocation id (start at 0):\n"
prog_info+="01: dislocation id (start at 0):\n"


prog_epilog="J.-M. Nocquet (Geoazur-CNRS) - May 2012"

parser = argparse.ArgumentParser(description=prog_info,epilog=prog_epilog)
parser.add_argument('-gps_h', action='store', dest='gmth',default=None,help='gmt psvelo files for horizontal displacement/velocity (mm or mm/y)')
parser.add_argument('-gps_u', action='store', dest='gmtu',default=None,help='files for vertical displacement/velocity - positive upwards with 4 columns: lon, lat, v_up, s_v_up (mm or mm/y)')
parser.add_argument('-insar', action='store', dest='insar',default=None,help='files for insar data. Format: Number xind yind east north data err wgt Elos Nlos Ulos (in metres)')
parser.add_argument('-g', action='store', dest='geometry',required=True,help='geometry file including dislocations and sources in npy format; see help for format')
parser.add_argument('-type', action='store',default='tde',help='dislocation elements type: tde (triangles) or rde (rectangles)')
parser.add_argument('-method', action='store',default='nikkhoo',help='method for Green tensor calculation: nikkhoo (default), meade (tde only) or edcmp (rde only)')
parser.add_argument('--verbose', '-v', action='count',default=0,help='verbose mode')
parser.add_argument('--debug', action='count',default=0,help='debug mode')
parser.add_argument('-e', action='store',type=str,required=True, dest='experiment',help='experiment name')

args = parser.parse_args()

if (len(sys.argv)<2):parser.print_help();sys.exit()

###############################################################################
# VERBOSE MODE
###############################################################################

if args.debug>0:
    logging.getLogger("my_logger").setLevel(logging.DEGUG)
    MESSAGE("verbose mode: DEBUG")
    verbose = True
    debug = True
else:
    debug = False
    if args.verbose > 0:
        logging.getLogger("my_logger").setLevel(logging.INFO)
        MESSAGE("verbose mode: VERBOSE")
        verbose = True
    else:
        logging.getLogger("my_logger").setLevel(logging.WARNING)
        MESSAGE("verbose mode: WARNING")
        verbose = False

###############################################################################
# CHECK TYPE AND METHOD COMPATIBILITY
###############################################################################

if ('tde' in args.type.lower()) and ('edcmp' in args.method.lower()):
    ERROR("edcmp method only is for rectangular dislocations",exit=True)

if ('rde' in args.type.lower()) and ('meade' in args.method.lower()):
    ERROR("meade method only is for triangular dislocations",exit=True)

###############################################################################
# READS GEOMETRY FILE AS NPY
###############################################################################

MESSAGE("Loading %s" % args.geometry )

if args.geometry[-3:] == 'dat':
    # reads the dat text file
    GEOMETRY , SGEOMETRY = pyeq.lib.geometry.to_np_array.dat_geometry_to_array_and_recarray( args.geometry , verbose=verbose )

if args.geometry[-3:] == 'npy':
    # reads the npy
    GEOMETRY , SGEOMETRY = pyeq.lib.geometry.to_np_array.npy_geometry_to_array_and_recarray( args.geometry , verbose=verbose )

MESSAGE("Number of subfaults: %d" % GEOMETRY.shape[0])

###############################################################################
# READS GMT FILE FOR HORIZONTAL COMPONENTS
###############################################################################


MESSAGE("Reading horizontal observation coordinates %s" % args.gmth)

OBS      = np.atleast_2d( np.genfromtxt(args.gmth,comments='#',usecols=(list(range(7)))) )
NAME_OBS = np.atleast_1d( np.genfromtxt(args.gmth,comments='#',usecols=(7),dtype=str) )

MESSAGE("Number of horizontal observations: %s" % OBS.shape[0])

# sorting
lindex = np.argsort( NAME_OBS )
OBS = OBS[ lindex ]
NAME_OBS = NAME_OBS[ lindex ]

###############################################################################
# READS UP FILE IF PROVIDED
###############################################################################

OBS_UP=np.zeros((0,0))
NAME_OBS_UP=np.zeros((0))

if args.gmtu != None:

    MESSAGE("Reading vertical observation coordinates %s" % args.gmtu)

    OBS_UP=np.genfromtxt(args.gmtu,comments='#',usecols=(list(range(4))))
    NAME_OBS_UP=np.genfromtxt(args.gmtu,comments='#',usecols=(4),dtype=str)

    MESSAGE("Number of vertical observations: %s" % OBS_UP.shape[0])

###############################################################################
# READS INSAR DATA FILE IF PROVIDED
###############################################################################

OBS_INSAR=np.zeros((0,0))

if args.insar != None:
    
    MESSAGE("Reading InSAR data %s " % args.insar)

    # data format must be index xind yind east north data err wgt Elos Nlos Ulos
    OBS_INSAR = np.genfromtxt(args.insar,comments='#',usecols = ( (3,4,5,6,8,9,10) ))
    OBS_INSAR[:,2] = OBS_INSAR[:,2]*1.E3
    OBS_INSAR[:,3] = OBS_INSAR[:,3]*1.E3
    # INSAR is now lon lat los sigma Elos Nlos Ulos with los in mm

    MESSAGE("Number of InSAR points : %d " % OBS_INSAR.shape[0])

###############################################################################
# CREATES GREEN FUNCTIONS FOR HORIZONTAL COMPONENTS
###############################################################################

n_dislocations=SGEOMETRY.shape[0]
n_gps=OBS.shape[0]
slip=1.0


# GREEN IS A TENSOR OF DIM 4
# GREEN(i,j,k,l) is the prediction for dislocation i at site j component k for rake l
# k=0,1,2 = east, north, up
# l=0,1 : rake_00 & rake_90

GREEN=np.zeros((n_dislocations, n_gps, 3,2))

if OBS.shape[0] > 0:
    MESSAGE("Creating Green tensor for GPS horizontal components")
    if 'meade' in args.method:
        GREEN = pyeq.green.make.meade_tde( GEOMETRY , OBS[:,:2], coor_type='geo' , disp=True, strain=False, stress= False, verbose=verbose )

    if ('tde' in args.type) and ('nikkhoo' in args.method):
        GREEN = pyeq.green.make.nikkhoo_tde( GEOMETRY , OBS[:,:2], coor_type='geo' , disp=True, strain=False, stress= False, verbose=verbose )

    if ('rde' in args.type) and ('nikkhoo' in args.method):
        GREEN = pyeq.green.make.nikkhoo_rde(GEOMETRY, OBS[:, :2], coor_type='geo', disp=True, strain=False,
                                            stress=False, verbose=verbose)

    if ('rde' in args.type) and ('edcmp' in args.method):
        GREEN, GREEN_STRAIN = pyeq.green.make.edcmp_rde( GEOMETRY , OBS[:,:2], coor_type='geo' , disp=True, strain=True, stress= False, verbose=verbose )

###############################################################################
# CREATES OBSERVATION FILE FOR VERTICAL COMPONENT IF PROVIDED
n_gps_up = OBS_UP.shape[0]
GREEN_UP = np.zeros((n_dislocations, n_gps_up, 3, 2))
if OBS_UP.shape[0] > 0:
###############################################################################

    MESSAGE("Creating Green tensor for the vertical component")
    if 'meade' in args.method:
        GREEN_UP = pyeq.green.make.meade_tde( GEOMETRY , OBS_UP[:,:2], coor_type='geo' , disp=True, strain=False, stress= False, verbose=verbose )

    if ('tde' in args.type) and ('nikkhoo' in args.method):
        GREEN_UP = pyeq.green.make.nikkhoo_tde( GEOMETRY , OBS_UP[:,:2], coor_type='geo' , disp=True, strain=False, stress= False, verbose=verbose )

    if ('rde' in args.type) and ('nikkhoo' in args.method):
        GREEN_UP = pyeq.green.make.nikkhoo_rde(GEOMETRY, OBS_UP[:, :2], coor_type='geo', disp=True, strain=False,
                                            stress=False, verbose=verbose)

    if ('rde' in args.type) and ('edcmp' in args.method):
        GREEN_UP, GREEN_STRAIN = pyeq.green.make.edcmp_rde( GEOMETRY , OBS_UP[:,:2], coor_type='geo' , disp=True, strain=True, stress= False, verbose=verbose )

###############################################################################
# CREATES GREEN TENSOR FOR INSAR IF PROVIDED
n_insar = OBS_INSAR.shape[0]
GREEN_INSAR = np.zeros((n_dislocations, n_insar, 3, 2))
if OBS_INSAR.shape[0] > 0:
###############################################################################

    MESSAGE("Creating Green tensor for the InSAR obs")

    if 'meade' in args.method:
        GREEN_INSAR = pyeq.green.make.meade_tde(GEOMETRY, OBS_INSAR[:, :2], coor_type='geo', disp=True, strain=False,
                                             stress=False, verbose=verbose)

    if ('tde' in args.type) and ('nikkhoo' in args.method):
        GREEN_INSAR = pyeq.green.make.nikkhoo_tde(GEOMETRY, OBS_INSAR[:, :2], coor_type='geo', disp=True, strain=False,
                                               stress=False, verbose=verbose)

    if ('rde' in args.type) and ('nikkhoo' in args.method):
        GREEN_INSAR = pyeq.green.make.nikkhoo_rde(GEOMETRY, OBS_INSAR[:, :2], coor_type='geo', disp=True, strain=False,
                                               stress=False, verbose=verbose)

    if ('rde' in args.type) and ('edcmp' in args.method):
        GREEN_INSAR, GREEN_STRAIN = pyeq.green.make.edcmp_rde(GEOMETRY, OBS_INSAR[:, :2], coor_type='geo', disp=True, strain=True,
                                                           stress=False, verbose=verbose)

###############################################################################
# MAKE THE DISTANCE MATRIX Dm
# THE DISTANCE MATRIX HAS SIZE mxm (m = number of model parameters)
# Dm is useful to calculate Cm
###############################################################################

MESSAGE("Now calculating the distance matrices")

nfaults=SGEOMETRY.shape[0]
Dm=np.zeros((nfaults,nfaults))

# test to speed up things - JMN 26/11/2020

(x_ref, y_ref, z_ref) = pyacs.lib.coordinates.geo2xyz(np.radians(GEOMETRY[:,9]), np.radians(GEOMETRY[:,10]), -np.fabs(GEOMETRY[:,11]))
import scipy.spatial

X = np.stack([x_ref,y_ref,z_ref]).T

Dm = scipy.spatial.distance_matrix(X,X) * 1.E-3

Dm_name=args.experiment+'_Dm.npy'

VERBOSE("Saving the relative Distance matrix: %s" % Dm_name)
try:
    np.save(Dm_name,Dm)
except:
    ERROR("Could not save %s " % Dm_name)

###############################################################################
# MAKE THE DISTANCE MATRIX D0 relative to the barycenter
# D0 is useful to calculate sigma
###############################################################################
D0_name=args.experiment+'_D0.npy'

MESSAGE("Creating the absolute Distance matrix %s"  % D0_name)

long0=np.mean(SGEOMETRY.centroid_long)
lat0=np.mean(SGEOMETRY.centroid_lat)
depth0=np.mean(SGEOMETRY.centroid_depth)

(x0,y0,z0)= pyacs.lib.coordinates.geo2xyz(np.radians(long0),np.radians(lat0),-np.fabs(depth0))

D0=np.zeros(nfaults)
for i in range(nfaults):
    DEBUG("fault #%05d over %05d" %(i,nfaults))

    (long_ref,lat_ref,depth_ref)=(SGEOMETRY.centroid_long[i],SGEOMETRY.centroid_lat[i],SGEOMETRY.centroid_depth[i])
    (x_ref,y_ref,z_ref)= pyacs.lib.coordinates.geo2xyz(np.radians(long_ref),np.radians(lat_ref),-np.fabs(depth_ref))

    distance_meters=np.sqrt(np.sum( (x_ref-x0)**2 + (y_ref-y0)**2 + (z_ref-z0)**2 ))
    distance_km=distance_meters / 1000.0

    D0[i]=distance_km

try:
    np.save(D0_name,D0)
except:
    ERROR("Could not save %s " % D0_name)

######################################################################################
# SAVE A BIG STRUCTURE ARRAY
# ORDER IS: GEOMETRY, Dm, ARRAY_SOURCES, DISLOCATIONS, GREEN, GREEN_UP,MAX_SLIP,OBS,NAME_OBS,OBS_UP,NAME_OBS_UP
######################################################################################

NPZ_OUT = ("%s_input.npz" % args.experiment)

MESSAGE("Saving the inversion structured array %s " % NPZ_OUT)

#np.save('GREEN.npy',GREEN)
#np.save('GREEN_UP.npy',GREEN_UP)
#np.save('GREEN_INSAR.npy',GREEN_INSAR)

# PYHTON 3.6
try:
    np.savez(NPZ_OUT, \
            GEOMETRY=GEOMETRY,\
            Dm=Dm,\
            GREEN=GREEN,\
            GREEN_UP=GREEN_UP,\
            OBS=OBS,\
            NAME_OBS=NAME_OBS,\
            OBS_UP=OBS_UP,\
            NAME_OBS_UP=NAME_OBS_UP,\
            GREEN_INSAR=GREEN_INSAR,\
            OBS_INSAR=OBS_INSAR,\
            )
except:
    ERROR("Could not save %s " % NPZ_OUT)


# check NPZ
try:
    npz = np.load(NPZ_OUT)
except:
    ERROR("Something went wrong. Could not load %s " % NPZ_OUT)
