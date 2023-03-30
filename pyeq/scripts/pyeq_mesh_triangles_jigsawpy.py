#!/usr/bin/env python

###################################################################
# SCRIPT    : pyacs_mesh_triangles_jigsawpy.py
# AUTHOR    :
# DATE      : December 2021
# INPUT     : grd_file
#           : experiment name
# OUTPUT    :
# NOTE      :
#
###################################################################


###################################################################
# MODULES IMPORT
###################################################################
import argparse, sys
import numpy as np
import logging
import scipy.interpolate
import netCDF4
import jigsawpy


import os
#import pandas as pd
from descartes import PolygonPatch
import matplotlib.pyplot as plt
import alphashape

from pyacs.lib import coordinates
#from pyacs.lib import icosahedron as Triangular_Mesh_Global
from icecream import ic

import pyeq.message.message as MESSAGE
import pyeq.message.verbose_message as VERBOSE
import pyeq.message.error as ERROR
import pyeq.message.warning as WARNING
import pyeq.message.debug_message as DEBUG

logging.getLogger("my_logger").setLevel(logging.INFO)

###################################################################
# PARSE ARGUMENT LINE
###################################################################

prog_info = "pygps_mesh_triangles_jigsawpy.py: discretizes a curved surface using quasi-equilateral triangles\n"
prog_info += "                     uses jigsawpy (https://github.com/dengwirda/jigsaw-python)\n"
prog_info += "                     First creates a triangular mesh covering the area = grid (intersection) bounds options\n"
prog_info += "                     Then for each center of the triangles, calculates strike, dip and rake (rake can be calculated using a relative Euler pole), and a ratio to be applied for slip scaling.\n"
prog_info += "                     The ratio equals the surface of the considered triangle / surface of the rectangle dislocation actually used\n"
prog_info += "-----------------------------------------------------------\n"
prog_info += "FORMAT FOR OUTPUT DISLOCATION FILE\n"
prog_info += "-----------------------------------------------------------\n"
prog_info += "\%05d \%10.5lf \%10.5lf \%8.1lf \%8.1lf \%8.2lf \%8.2lf \%8.2lf \%8.2lf\n"
prog_info += " index long lat depth strike dip length width rake ratio\n"

prog_epilog = "J.-M. Nocquet (Geoazur-IRD) - December 2021"

parser = argparse.ArgumentParser(description=prog_info, epilog=prog_epilog)
parser.add_argument('-g', action='store', dest='grd', type=str, required=True, help='Netcdf grid, likely slab1 or slab2 or custom grid')
parser.add_argument('-b', action='store', dest='bounds', type=str, required=True,
                    help='Bounds: /min_lon/max_lon/min_lat/max_lat , be careful about -180/+180 and 0/360 grids')
parser.add_argument('-l', action='store', dest='ltedge', required=True, type=float,
                    help='length of triangle edge in km')
parser.add_argument('-d', action='store', dest='depth_range', required=True, type=str,
                    help='Depth range min_depth/max_depth in km')
parser.add_argument('-alpha', action='store', dest='alpha_value', default=0., type=float,
                    help='alpha parameter to simplify the polygon defining the mesh boundary')
parser.add_argument('-gps', action='store', dest='gmth', default=None, type=str,
                    help='gps psvelo gmt file for mesh optimization according to resolution')
parser.add_argument('--verbose', '-v', action='count', default=0, help='verbose mode')
parser.add_argument('--debug', action='count', default=0, help='debug mode')
parser.add_argument('-e', action='store', type=str, required=True, dest='experiment', help='experiment name')

args = parser.parse_args()

if (len(sys.argv) < 2): parser.print_help();sys.exit()


###################################################################
# SET VERBOSE LEVEL
###################################################################
if args.verbose > 0:
    MESSAGE("Verbose mode")
    verbose = True
else:
    verbose = False

# debug
if args.debug > 0:
    MESSAGE("Debug mode")
    verbose = True
    debug = True
else:
    debug = False


#Rt = 6371.0E3

if not(verbose):
    # WARNING LEVEL: ONLY MAJOR STEPS AND WARNINGS WILL BE PRINTED
    logging.getLogger("my_logger").setLevel(logging.WARNING)
else:
    if verbose:
        # VERBOSE MODE
        logging.getLogger("my_logger").setLevel(logging.INFO)
    else:
        # WARNING MODE
        logging.getLogger("my_logger").setLevel(logging.WARNING)

if debug:
    # DEBUG MODE
    logging.getLogger("my_logger").setLevel(logging.DEBUG)

if logging.getLogger("my_logger").getEffectiveLevel() == logging.WARNING:
    MESSAGE("verbose level: WARNING ONLY")

if logging.getLogger("my_logger").getEffectiveLevel() == logging.INFO:
    MESSAGE("verbose level: VERBOSE")

if logging.getLogger("my_logger").getEffectiveLevel() == logging.DEBUG:
    MESSAGE("verbose level: DEBUG")

###################################################################
# CHECK INPUT ARGUMENT
###################################################################


(min_depth, max_depth) = list(map(float, args.depth_range.split('/')))
DEBUG("Bounds (long. in/max lat. min/max): %s" % " ".join(args.bounds.split('/')[1:]))
(lon_min, lon_max, lat_min, lat_max) = list(map(float, args.bounds.split('/')[1:]))
(rlon_min, rlon_max, rlat_min, rlat_max) = list(map(np.radians, (lon_min, lon_max, lat_min, lat_max)))

###################################################################
# READS THE GRID
###################################################################

# 1. creates a regular rectangular grid over using user's provided bounds
x = np.arange(lon_min,lon_max, 0.1)
y = np.arange(lat_min,lat_max, 0.1)
xv, yv = np.meshgrid(x, y, sparse=False, indexing='ij')
a = np.array((xv.flatten(),yv.flatten())).T

#2. reads grid and interpolate depth
try:
    ds = netCDF4.Dataset( args.grd )
except:
    ERROR( ("Could not read grid: %s" % (args.grd)),exit=True )

z = ds.variables['z'][:].T
# change JMN 28/03/2023 - account for 0-360 grid like slab2
longitudes = ds.variables['x'][:]
lidx = np.where( longitudes > 180 )[0]
longitudes[lidx] = longitudes[lidx] - 360.
interp = scipy.interpolate.RegularGridInterpolator(
    tuple((longitudes, ds.variables['y'][:])),
    z.data,
    method='linear',
    bounds_error=False)
############################################################
#interp = scipy.interpolate.RegularGridInterpolator(
#    tuple((ds.variables['x'][:], ds.variables['y'][:])),
#    z.data,
#    method='linear',
#    bounds_error=False)

interpolated = interp( (a[:,0], a[:,1]) )

aa = np.c_[a,interpolated.flatten()]

#3. defines points within the grid and within user's provided bounds
# inside grid
lindex = np.where( ~np.isnan(aa[:,2]) )[0]
aa = aa[lindex,:]
# depth range
lindex = np.where( np.fabs(aa[:,2])<max_depth )[0]
aa = aa[lindex,:]
lindex = np.where( np.fabs(aa[:,2])>min_depth )[0]
aa = aa[lindex,:]
# lon bounds
lindex = np.where( aa[:,0]>lon_min )[0]
aa = aa[lindex,:]
lindex = np.where( aa[:,0]<lon_max )[0]
aa = aa[lindex,:]
# lat bounds
lindex = np.where( aa[:,1]>lat_min )[0]
aa = aa[lindex,:]
lindex = np.where( aa[:,1]<lat_max )[0]
aa = aa[lindex,:]

if aa.shape[0] == 0:
    ERROR("grid bounds and user bounds do not overlap" , exit=True )

if debug:
    import matplotlib.pyplot as plt
    plt.scatter(aa[:, 0], aa[:, 1], c=aa[:, 2], cmap='jet', s=0.1)
    plt.colorbar()
    plt.show()

###################################################################
# use alphashape to get the discretized area enveloppe
###################################################################

points_2d = []
for i in np.arange( aa.shape[0] ):
  points_2d.append( (aa[i,0],aa[i,1]) )

alpha_shape = alphashape.alphashape(points_2d, args.alpha_value)
boundary_xy = np.array(alpha_shape.exterior.coords.xy).T
if debug:
    alpha_shape
    fig, ax = plt.subplots()
    ax.scatter(*zip(*points_2d),s=0.1)
    ax.scatter(boundary_xy[:,0],boundary_xy[:,1],s=3,color='r')
    ax.add_patch(PolygonPatch(alpha_shape, alpha=0.2))
    ax.set_aspect('equal')
    plt.show()

boundary_xy = np.array(alpha_shape.exterior.coords.xy).T

###################################################################
# PREPARE JIGSAWPY
###################################################################

lnode = []
for i in np.arange( boundary_xy.shape[0] ):
  lnode.append( ((boundary_xy[i,0],boundary_xy[i,1]),0) )

ledge = []
for i in np.arange( len(lnode) ):
    ledge.append( ((i,i+1),0) )
ledge[-1] = ((i,0),0)

#------------------------------------ init JIGSAW

opts = jigsawpy.jigsaw_jig_t()

geom = jigsawpy.jigsaw_msh_t()
hmat = jigsawpy.jigsaw_msh_t()
mesh = jigsawpy.jigsaw_msh_t()

#------------------------------------ define JIGSAW geometry

#geom.mshID = "ellipsoid-mesh"
geom.mshID = "euclidean-mesh"
geom.ndims = +2

geom.vert2 = np.array(lnode,dtype=geom.VERT2_t)
geom.edge2 = np.array(ledge,dtype=geom.EDGE2_t)

#------------------------------------ define JIGSAW options

opts.hfun_hmax = args.ltedge / 111. /4.             # push HFUN limits
opts.hfun_hmin = args.ltedge / 111. /8.              # push HFUN limits
opts.mesh_dims = +2                 # 2-dim. simplexes
opts.optm_qlim = +.95
#opts.mesh_top1 = True               # for sharp feat's
#opts.geom_feat = True

#------------------------------------ run JIGSAW
MESSAGE(("running JIGSAW with edge length constraint %.3lf" % args.ltedge))
jigsawpy.lib.jigsaw(opts, geom, mesh)

###################################################################
# SAVE RESULTS
###################################################################

vtkfile = args.experiment+'.vtk'
meshfile = args.experiment+'.msh'
MESSAGE("Saving %s" % vtkfile)
jigsawpy.savevtk(vtkfile, mesh)
MESSAGE("Saving %s" % meshfile)
jigsawpy.savemsh(meshfile, mesh)

###################################################################
# JIGSAW TO PYACS FORMAT
###################################################################
nvert = mesh.vert2.shape[0]
ntria = mesh.tria3.shape[0]
nedge = mesh.edge2.shape[0]
MESSAGE("mesh includes %d vertices and %d triangles" % (nvert,ntria)  )
# vertices
vert = np.zeros(( nvert , 3 ))

for i in np.arange( nvert ):
    vert[i,:2] = mesh.vert2[i][0]
    if debug:ic(i,mesh.vert2[i][0])
# edges
edge = np.zeros(( nedge , 2 ))
for i in np.arange( nedge ):
    edge[i,:2] = mesh.edge2[i][0]

# fill depth from grid
# change JMN 28/03/2023 to accommodate 0-360 grid
interp = scipy.interpolate.RegularGridInterpolator(
    tuple((longitudes, ds.variables['y'][:])),
    z.data,
    method='linear',
    bounds_error=True)
#interp = scipy.interpolate.RegularGridInterpolator(
#    tuple((ds.variables['x'][:], ds.variables['y'][:])),
#    z.data,
#    method='linear',
#    bounds_error=True)
vert[:,2] = interp( (vert[:,0], vert[:,1]) )

# check vertices on the mesh edges
if np.isnan( vert[:,2]).any():
    MESSAGE("Handling vertices slighlty out of the grid")
    wedge = np.copy(edge)
    ledge = []
    lvertice = []
    ledge.append(wedge[0])
    wedge = np.delete(wedge, 0, axis=0)
    next_vertex = ledge[0][-1]
    init_vertex = ledge[0][0]
    lvertice.append(next_vertex)

    while next_vertex != init_vertex:

        index = np.where(wedge == next_vertex)[0]
        ledge.append(wedge[index][0])
        if wedge[index][-1][0] != next_vertex:
            next_vertex = wedge[index][-1][0]
        else:
            next_vertex = wedge[index][-1][1]
        lvertice.append(next_vertex)
        wedge = np.delete(wedge, index, axis=0)

    MESSAGE("Mesh boundary includes %d vertices" % (len(lvertice)))

    np_vert_edge = np.array( lvertice ,dtype=int )
    lindex_nan = np.where( np.isnan(vert[np_vert_edge,2]) )[0]
    MESSAGE("There are %d vertices with depth problems" % (len(lindex_nan)))
    lindex_not_nan = np.where( ~np.isnan(vert[np_vert_edge,2]) )[0]

    if debug:ic(np.arange(np_vert_edge.shape[0])[lindex_nan])
    if debug:ic(np.arange(np_vert_edge.shape[0])[lindex_not_nan])
    if debug:ic(vert[np_vert_edge[lindex_not_nan]])

    vert[np_vert_edge[lindex_nan],2] = np.interp(np.arange(np_vert_edge.shape[0])[lindex_nan],\
                                                 np.arange(np_vert_edge.shape[0])[lindex_not_nan],\
                                                 vert[np_vert_edge[lindex_not_nan],2])
    if debug:ic(vert[np_vert_edge])

if debug:
    import matplotlib.pyplot as plt
    plt.scatter(vert[np_vert_edge, 0], vert[np_vert_edge, 1], s=5)
    plt.colorbar()
    plt.show()


# check nan
if np.isnan( vert[:,2]).any():
    lindex= np.where( np.isnan(vert[:,2]) )[0]
    if debug:ic( vert[lindex,:] )
    ERROR("Some vertices had NaN depth." , exit=True)

# triangles
faces = np.zeros((ntria,3),dtype=int)
for i in np.arange( ntria ):
    faces[i,:] = mesh.tria3[i][0]
MESSAGE("Computing area, strike and dip")

####################################################################
# GEOMETRY
####################################################################

geometry = np.zeros(( ntria, 22 ))

# set depth positive
vert[:,2] = np.fabs(vert[:,2])

for i in np.arange( ntria ):
    # get vertice
    [idx_v1, idx_v2, idx_v3] = mesh.tria3[i][0]
    AA = np.array(coordinates.geo2xyz(vert[idx_v1][0], vert[idx_v1][1], -vert[idx_v1][2]*1.E3,unit='dec_deg'))
    BB = np.array(coordinates.geo2xyz(vert[idx_v2][0], vert[idx_v2][1], -vert[idx_v2][2]*1.E3,unit='dec_deg'))
    CC = np.array(coordinates.geo2xyz(vert[idx_v3][0], vert[idx_v3][1], -vert[idx_v3][2]*1.E3,unit='dec_deg'))

    if debug:
        ic(i)
        ic(vert[idx_v1])
        ic(vert[idx_v2])
        ic(vert[idx_v3])
        ic(AA)
        ic(BB)
        ic(CC)

    # Normal vector

    N = np.cross((BB - AA), (CC - AA))
    length = np.sqrt(N[0] ** 2 + N[1] ** 2 + N[2] ** 2)
    n = N / length

    # area in km^2 S=0.5 * norm (B-A x C-A)
    area = .5 * length / 1.E6
    if debug:ic(area)
    # barycenter of the triangle
    M = (AA + BB + CC) / 3

    # rotate to local frame in ENU
    (l, phi, he) = coordinates.xyz2geo(M[0], M[1], M[2])
    if debug:ic(np.degrees([l,phi]))
    R = coordinates.mat_rot_general_to_local(l, phi)
    ENU = np.dot(R, n)

    # we want the normal vector to be upward
    if ENU[2] < 0: ENU = -ENU

    # dip
    dip = np.degrees(np.arctan((np.sqrt(ENU[0] ** 2 + ENU[1] ** 2)) / ENU[2]))
    if debug:ic(dip)
    # usv the unit strike vector
    usv = np.array([-ENU[1], ENU[0], 0.]) / np.sqrt(ENU[0] ** 2 + ENU[1] ** 2)
    # get strike
    strike = np.degrees(np.arctan2(usv[0], usv[1]))
    if debug:ic(strike)
    if strike > 90.0 or strike < -90.0:
        if verbose: ic('ENU strike dip lp ', ENU, strike, dip, np.degrees(l), np.degrees(phi))

    # 0,1,2:rdis_long,rdis_lat,rdis_depth
    # 3,4:rdis_length,rdis_width
    # 5:rdis_area
    # 6:ratio_rdis_tdis
    # 7:strike
    # 8:dip
    # 9,10,11:centroid_long,centroid_lat,centroid_depth
    # 12,13,14:tdis_long1,tdis_lat1,tdis_depth1,
    # 15,16,17:tdis_long2,tdis_lat2,tdis_depth2
    # 18,19,20: tdis_long3,tdis_lat3,tdis_depth3
    # 21:tdis_area

    geometry[i,7] = strike
    geometry[i,8] = dip
    geometry[i,9:11] = np.degrees([l,phi])
    geometry[i,11] = (vert[idx_v1,2]+vert[idx_v2,2]+vert[idx_v3,2])/3.
    geometry[i,12:15] = vert[idx_v1]
    geometry[i,15:18] = vert[idx_v2]
    geometry[i,18:21] = vert[idx_v3]
    geometry[i,21] = area

# print median edge length
median_area = np.median(geometry[:,21])
min_area = np.min(geometry[:,21])
max_area = np.max(geometry[:,21])

median_edge = np.sqrt(median_area*4/np.sqrt(3))
min_edge = np.sqrt(min_area*4/np.sqrt(3))
max_edge = np.sqrt(max_area*4/np.sqrt(3))

MESSAGE("median min and max edge length in km: %.1lf %.1lf %.1lf" % (median_edge,min_edge,max_edge) )

f_geometry=args.experiment+'_geometry.dat'
f_geometry_npy=args.experiment+'_geometry.npy'

MESSAGE("writing geometry files: %s %s" % (f_geometry, f_geometry_npy))

#print("-- Sorting geometry array by descending latitudes")
#GEOMETRY = geometry[geometry[:, 10].argsort()[::-1]]
GEOMETRY = geometry

names = ['rdis_long', 'rdis_lat', 'rdis_depth', 'rdis_length', 'rdis_width', \
         'rdis_area', 'ratio_rdis_tdis', 'strike', 'dip', \
         'centroid_long', 'centroid_lat', 'centroid_depth', \
         'tdis_long1', 'tdis_lat1', 'tdis_depth1', \
         'tdis_long2', 'tdis_lat2', 'tdis_depth2', \
         'tdis_long3', 'tdis_lat3', 'tdis_depth3', 'tdis_area']

np.save(f_geometry_npy, GEOMETRY)

Hdtype={'names':['    rdis_long','   rdis_lat','rdis_depth',\
                'rdis_length','rdis_width',' rdis_area','ratio_rdis_tdis',\
                '   strike','       dip',\
                'centroid_long','centroid_lat','centroid_depth',\
                'tdis_long1', ' tdis_lat1','tdis_depth1', \
                'tdis_long2',' tdis_lat2','tdis_depth2', \
                'tdis_long3',' tdis_lat3','tdis_depth3', \
                ' tdis_area'],\
         'formats':['f','f','f','f','f','f','f','f','f','f','f','f','f','f','f','f','f','f','f','f','f','f']}

GEOMETRY_TXT = np.zeros((GEOMETRY.shape[0], GEOMETRY.shape[1] + 1))
GEOMETRY_TXT[:, 0] = np.arange(GEOMETRY.shape[0])
GEOMETRY_TXT[:, 1:] = GEOMETRY
header_cols = ' '.join(Hdtype['names'])
format_header = '%04d %10.5lf %10.5lf      %6.2lf \
     %6.2lf     %6.2lf     %6.2lf          %6.2lf\
    %6.2lf %10.2lf \
   %10.5lf   %10.5lf         %6.2lf \
%10.5lf %10.5lf      %6.2lf \
%10.5lf %10.5lf      %6.2lf \
%10.5lf %10.5lf      %6.2lf \
   %6.2lf '

np.savetxt(f_geometry, GEOMETRY_TXT, fmt=format_header, header=header_cols)

#fsources.close()
#fdislocations.close()
#frectangular_gmt.close()

# print gmt and shapefiles

from pyeq.log.geometry2shp_gmt import geometry2shp_gmt

geometry2shp_gmt(GEOMETRY, 'tde', out_shp=args.experiment + '_tde', out_gmt=args.experiment + '_tde.gmt',
                 verbose=verbose)
# geometry2shp_gmt(GEOMETRY, 'rde', out_shp=args.experiment + '_rde', out_gmt=args.experiment + '_rde.gmt',
#                  verbose=verbose)
#
# # test save verts and faces
#
# # np.savetxt('faces.dat',np.array(faces))
# # np.savetxt('verts.dat',np.array(verts))
#
#

if args.gmth is None: sys.exit()

####################################################################
# COMPUTE GREEN FUNCTION AND RESOLUTION MATRIX
####################################################################

import pyeq.green.make

###############################################################################
# READS GMT FILE FOR HORIZONTAL COMPONENTS
###############################################################################


MESSAGE("Reading horizontal observation coordinates %s" % args.gmth)

OBS      = np.atleast_2d( np.genfromtxt(args.gmth,comments='#',usecols=(list(range(7)))) )
NAME_OBS = np.atleast_1d( np.genfromtxt(args.gmth,comments='#',usecols=(7),dtype=str) )
# sorting
lindex = np.argsort( NAME_OBS )
OBS = OBS[ lindex ]
NAME_OBS = NAME_OBS[ lindex ]
MESSAGE("Number of horizontal observations: %s" % OBS.shape[0])

n_dislocations=GEOMETRY.shape[0]
n_gps=OBS.shape[0]
slip=1.0

GREEN=np.zeros((n_dislocations, n_gps, 3,2))

MESSAGE("Creating Green tensor")
GREEN = pyeq.green.make.nikkhoo_tde( GEOMETRY , OBS[:,:2], coor_type='geo' , disp=True, strain=False, stress= False, verbose=verbose )

MESSAGE("Creating resolution map")

resolution = np.sqrt( np.sum( (GREEN*1.E3)**2 , axis=(1,2,3) ))
MESSAGE("Best  resolved triangle %.1lf mm" % (np.max( resolution )) )
MESSAGE("Worst resolved triangle %.1lf mm" % (np.min( resolution )) )
MESSAGE("Ratio %.1lf and log %.1lf " % ( np.max( resolution) / np.min( resolution) , np.log10(np.max( resolution) / np.min( resolution))  ))

ic(resolution.shape)
np.savetxt('spatial_resolution.dat', np.c_[geometry[:, 9:11], resolution],
           fmt="%10.5lf %10.5lf %10.3E")
GEOMETRY[:,-1] = resolution
geometry2shp_gmt(GEOMETRY, 'tde', out_shp=args.experiment + '_resolution', out_gmt=args.experiment + '_resolution.gmt',
                 verbose=verbose)

###############################################################################
# MAKES IMPROVED MESH
###############################################################################

# ------------------------------------ compute HFUN over GEOM
# resolution to length constraint
lcons = np.sqrt(1./resolution)
# normalize
lcons = lcons / np.min(lcons) * args.ltedge / 111.

ic(np.min(lcons))
ic(np.max(lcons))

ic(np.min(lcons)*111.)
ic(np.max(lcons)*111.)

np_sr = np.c_[geometry[:, 9:11], lcons]

from scipy.interpolate import griddata

# define grid.
xi = np.arange(np.min(np_sr[:, 0]), np.max(np_sr[:, 0]), median_edge / 111. /2. )
yi = np.arange(np.min(np_sr[:, 1]), np.max(np_sr[:, 1]), median_edge / 111. /2. )
zi = griddata((np_sr[:, 0], np_sr[:, 1]), np_sr[:, 2], (xi[None, :], yi[:, None]), method='linear',
              fill_value=0.5)

ic(xi.shape)
ic(yi.shape)
ic(zi.shape)

hmat.mshID = "euclidean-grid"
hmat.ndims = +2
hmat.xgrid = np.array(
    xi, dtype=hmat.REALS_t)
hmat.ygrid = np.array(
    yi, dtype=hmat.REALS_t)
hmat.value = np.array(
    zi, dtype=hmat.REALS_t)

# ------------------------------------ build mesh via JIGSAW!

#print("Call libJIGSAW: case 0c.")

opts.hfun_scal = "absolute"
opts.hfun_hmax = float("inf")  # null HFUN limits
#opts.hfun_hmin = float(+0.05)

opts.mesh_dims = +2  # 2-dim. simplexes

opts.optm_qlim = +.95

#opts.mesh_top1 = True  # for sharp feat's
#opts.geom_feat = True

MESSAGE("Running jigsawpy.lib.jigsaw")
jigsawpy.lib.jigsaw(opts, geom, mesh,None, hmat)


jigsawpy.savevtk('test_ireso.vtk', mesh)
jigsawpy.savemsh("test_ireso.msh", mesh)

###################################################################
# JIGSAW TO PYACS FORMAT
###################################################################
nvert = mesh.vert2.shape[0]
ntria = mesh.tria3.shape[0]
nedge = mesh.edge2.shape[0]
MESSAGE("mesh includes %d vertices and %d triangles" % (nvert,ntria)  )
# vertices
vert = np.zeros(( nvert , 3 ))

for i in np.arange( nvert ):
    vert[i,:2] = mesh.vert2[i][0]
    if debug:ic(i,mesh.vert2[i][0])
# edges
edge = np.zeros(( nedge , 2 ))
for i in np.arange( nedge ):
    edge[i,:2] = mesh.edge2[i][0]

# fill depth from grid
interp = scipy.interpolate.RegularGridInterpolator(
    tuple((ds.variables['x'][:], ds.variables['y'][:])),
    z.data,
    method='linear',
    bounds_error=True)
vert[:,2] = interp( (vert[:,0], vert[:,1]) )

# check vertices on the mesh edges
if np.isnan( vert[:,2]).any():
    MESSAGE("Handling vertices slighlty out of the grid")
    wedge = np.copy(edge)
    ledge = []
    lvertice = []
    ledge.append(wedge[0])
    wedge = np.delete(wedge, 0, axis=0)
    next_vertex = ledge[0][-1]
    init_vertex = ledge[0][0]
    lvertice.append(next_vertex)

    while next_vertex != init_vertex:

        index = np.where(wedge == next_vertex)[0]
        ledge.append(wedge[index][0])
        if wedge[index][-1][0] != next_vertex:
            next_vertex = wedge[index][-1][0]
        else:
            next_vertex = wedge[index][-1][1]
        lvertice.append(next_vertex)
        wedge = np.delete(wedge, index, axis=0)

    MESSAGE("Mesh boundary includes %d vertices" % (len(lvertice)))

    np_vert_edge = np.array( lvertice ,dtype=int )
    lindex_nan = np.where( np.isnan(vert[np_vert_edge,2]) )[0]
    MESSAGE("There are %d vertices with depth problems" % (len(lindex_nan)))
    lindex_not_nan = np.where( ~np.isnan(vert[np_vert_edge,2]) )[0]

    if debug:ic(np.arange(np_vert_edge.shape[0])[lindex_nan])
    if debug:ic(np.arange(np_vert_edge.shape[0])[lindex_not_nan])
    if debug:ic(vert[np_vert_edge[lindex_not_nan]])

    vert[np_vert_edge[lindex_nan],2] = np.interp(np.arange(np_vert_edge.shape[0])[lindex_nan],\
                                                 np.arange(np_vert_edge.shape[0])[lindex_not_nan],\
                                                 vert[np_vert_edge[lindex_not_nan],2])
    if debug:ic(vert[np_vert_edge])

if debug:
    import matplotlib.pyplot as plt
    plt.scatter(vert[np_vert_edge, 0], vert[np_vert_edge, 1], s=5)
    plt.colorbar()
    plt.show()


# check nan
if np.isnan( vert[:,2]).any():
    lindex= np.where( np.isnan(vert[:,2]) )[0]
    if debug:ic( vert[lindex,:] )
    ERROR("Some vertices had NaN depth." , exit=True)

# triangles
faces = np.zeros((ntria,3),dtype=int)
for i in np.arange( ntria ):
    faces[i,:] = mesh.tria3[i][0]
MESSAGE("Computing area, strike and dip")

####################################################################
# GEOMETRY
####################################################################

geometry = np.zeros(( ntria, 22 ))

# set depth positive
vert[:,2] = np.fabs(vert[:,2])

for i in np.arange( ntria ):
    # get vertice
    [idx_v1, idx_v2, idx_v3] = mesh.tria3[i][0]
    AA = np.array(coordinates.geo2xyz(vert[idx_v1][0], vert[idx_v1][1], -vert[idx_v1][2]*1.E3,unit='dec_deg'))
    BB = np.array(coordinates.geo2xyz(vert[idx_v2][0], vert[idx_v2][1], -vert[idx_v2][2]*1.E3,unit='dec_deg'))
    CC = np.array(coordinates.geo2xyz(vert[idx_v3][0], vert[idx_v3][1], -vert[idx_v3][2]*1.E3,unit='dec_deg'))

    if debug:
        ic(i)
        ic(vert[idx_v1])
        ic(vert[idx_v2])
        ic(vert[idx_v3])
        ic(AA)
        ic(BB)
        ic(CC)

    # Normal vector

    N = np.cross((BB - AA), (CC - AA))
    length = np.sqrt(N[0] ** 2 + N[1] ** 2 + N[2] ** 2)
    n = N / length

    # area in km^2 S=0.5 * norm (B-A x C-A)
    area = .5 * length / 1.E6
    if debug:ic(area)
    # barycenter of the triangle
    M = (AA + BB + CC) / 3

    # rotate to local frame in ENU
    (l, phi, he) = coordinates.xyz2geo(M[0], M[1], M[2])
    if debug:ic(np.degrees([l,phi]))
    R = coordinates.mat_rot_general_to_local(l, phi)
    ENU = np.dot(R, n)

    # we want the normal vector to be upward
    if ENU[2] < 0: ENU = -ENU

    # dip
    dip = np.degrees(np.arctan((np.sqrt(ENU[0] ** 2 + ENU[1] ** 2)) / ENU[2]))
    if debug:ic(dip)
    # usv the unit strike vector
    usv = np.array([-ENU[1], ENU[0], 0.]) / np.sqrt(ENU[0] ** 2 + ENU[1] ** 2)
    # get strike
    strike = np.degrees(np.arctan2(usv[0], usv[1]))
    if debug:ic(strike)
    if strike > 90.0 or strike < -90.0:
        if verbose: ic('ENU strike dip lp ', ENU, strike, dip, np.degrees(l), np.degrees(phi))

    # 0,1,2:rdis_long,rdis_lat,rdis_depth
    # 3,4:rdis_length,rdis_width
    # 5:rdis_area
    # 6:ratio_rdis_tdis
    # 7:strike
    # 8:dip
    # 9,10,11:centroid_long,centroid_lat,centroid_depth
    # 12,13,14:tdis_long1,tdis_lat1,tdis_depth1,
    # 15,16,17:tdis_long2,tdis_lat2,tdis_depth2
    # 18,19,20: tdis_long3,tdis_lat3,tdis_depth3
    # 21:tdis_area

    geometry[i,7] = strike
    geometry[i,8] = dip
    geometry[i,9:11] = np.degrees([l,phi])
    geometry[i,11] = (vert[idx_v1,2]+vert[idx_v2,2]+vert[idx_v3,2])/3.
    geometry[i,12:15] = vert[idx_v1]
    geometry[i,15:18] = vert[idx_v2]
    geometry[i,18:21] = vert[idx_v3]
    geometry[i,21] = area

# print median edge length
median_area = np.median(geometry[:,21])
min_area = np.min(geometry[:,21])
max_area = np.max(geometry[:,21])

median_edge = np.sqrt(median_area*4/np.sqrt(3))
min_edge = np.sqrt(min_area*4/np.sqrt(3))
max_edge = np.sqrt(max_area*4/np.sqrt(3))

MESSAGE("median min and max edge length in km: %.1lf %.1lf %.1lf" % (median_edge,min_edge,max_edge) )

f_geometry=args.experiment+'ireso_geometry.dat'
f_geometry_npy=args.experiment+'ireso_geometry.npy'

MESSAGE("writing geometry files: %s %s" % (f_geometry, f_geometry_npy))

#print("-- Sorting geometry array by descending latitudes")
#GEOMETRY = geometry[geometry[:, 10].argsort()[::-1]]
GEOMETRY = geometry

names = ['rdis_long', 'rdis_lat', 'rdis_depth', 'rdis_length', 'rdis_width', \
         'rdis_area', 'ratio_rdis_tdis', 'strike', 'dip', \
         'centroid_long', 'centroid_lat', 'centroid_depth', \
         'tdis_long1', 'tdis_lat1', 'tdis_depth1', \
         'tdis_long2', 'tdis_lat2', 'tdis_depth2', \
         'tdis_long3', 'tdis_lat3', 'tdis_depth3', 'tdis_area']

np.save(f_geometry_npy, GEOMETRY)

Hdtype={'names':['    rdis_long','   rdis_lat','rdis_depth',\
                'rdis_length','rdis_width',' rdis_area','ratio_rdis_tdis',\
                '   strike','       dip',\
                'centroid_long','centroid_lat','centroid_depth',\
                'tdis_long1', ' tdis_lat1','tdis_depth1', \
                'tdis_long2',' tdis_lat2','tdis_depth2', \
                'tdis_long3',' tdis_lat3','tdis_depth3', \
                ' tdis_area'],\
         'formats':['f','f','f','f','f','f','f','f','f','f','f','f','f','f','f','f','f','f','f','f','f','f']}

GEOMETRY_TXT = np.zeros((GEOMETRY.shape[0], GEOMETRY.shape[1] + 1))
GEOMETRY_TXT[:, 0] = np.arange(GEOMETRY.shape[0])
GEOMETRY_TXT[:, 1:] = GEOMETRY
header_cols = ' '.join(Hdtype['names'])
format_header = '%04d %10.5lf %10.5lf      %6.2lf \
     %6.2lf     %6.2lf     %6.2lf          %6.2lf\
    %6.2lf %10.2lf \
   %10.5lf   %10.5lf         %6.2lf \
%10.5lf %10.5lf      %6.2lf \
%10.5lf %10.5lf      %6.2lf \
%10.5lf %10.5lf      %6.2lf \
   %6.2lf '

np.savetxt(f_geometry, GEOMETRY_TXT, fmt=format_header, header=header_cols)

#fsources.close()
#fdislocations.close()
#frectangular_gmt.close()

# print gmt and shapefiles

from pyeq.log.geometry2shp_gmt import geometry2shp_gmt

geometry2shp_gmt(GEOMETRY, 'tde', out_shp=args.experiment + 'ireso_tde', out_gmt=args.experiment + '_tde.gmt',
                 verbose=verbose)

n_dislocations=GEOMETRY.shape[0]
n_gps=OBS.shape[0]
slip=1.0

GREEN=np.zeros((n_dislocations, n_gps, 3,2))

MESSAGE("Creating Green tensor")
GREEN = pyeq.green.make.nikkhoo_tde( GEOMETRY , OBS[:,:2], coor_type='geo' , disp=True, strain=False, stress= False, verbose=verbose )

MESSAGE("Creating resolution map")

resolution = np.sqrt( np.sum( (GREEN*1.E3)**2 , axis=(1,2,3) ))
MESSAGE("Best  resolved triangle %.1lf mm" % (np.max( resolution )) )
MESSAGE("Worst resolved triangle %.1lf mm" % (np.min( resolution )) )
MESSAGE("Ratio %.1lf and log %.1lf " % ( np.max( resolution) / np.min( resolution) , np.log10(np.max( resolution) / np.min( resolution))  ))

ic(resolution.shape)
np.savetxt('ireso.dat', np.c_[geometry[:, 9:11], resolution],
           fmt="%10.5lf %10.5lf %10.3E")
GEOMETRY[:,-1] = resolution
geometry2shp_gmt(GEOMETRY, 'tde', out_shp=args.experiment + '_ireso_resolution', out_gmt=args.experiment + '_resolution.gmt',
                 verbose=verbose)
