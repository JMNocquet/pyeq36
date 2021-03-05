#!/usr/bin/env python

###################################################################
# SCRIPT    : pyacs_parametrize_curve_surface_triangles.py
# AUTHOR    : 
# DATE      : December 2011
# INPUT     :  grd_file
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

from pyacs.lib import coordinates
from pyacs.lib import icosahedron as Triangular_Mesh_Global

###################################################################
# PARSE ARGUMENT LINE
###################################################################

prog_info="pygps_parametrize_curve_surface_triangles.py:discretizes a curved surface using equilateral triangles and produces a file to be read by pyeq_model_make_green \n"
prog_info+="                     First creates a triangular mesh covering the area = grid (intersection) bounds options\n"
prog_info+="                     Then for each center of the triangles, calculates strike, dip and rake (rake can be calculated using a relative Euler pole), and a ratio to be applied for slip scaling.\n"
prog_info+="                     The ratio equals the surface of the considered triangle / surface of the rectangle dislocation actually used\n";
prog_info+="-----------------------------------------------------------\n"
prog_info+="FORMAT FOR OUTPUT DISLOCATION FILE\n"
prog_info+="-----------------------------------------------------------\n"
prog_info+="\%05d \%10.5lf \%10.5lf \%8.1lf \%8.1lf \%8.2lf \%8.2lf \%8.2lf \%8.2lf\n" 
prog_info+=" index long lat depth strike dip length width rake ratio\n"

prog_epilog="J.-M. Nocquet (Geoazur-CNRS) - October 2012"

parser = argparse.ArgumentParser(description=prog_info,epilog=prog_epilog)
parser.add_argument('-g', action='store', dest='grd',type=str,required=True,help='Netcdf grid')
parser.add_argument('-b', action='store', dest='bounds',type=str,required=True,help='Bounds: /min_lon/max_lon/min_lat/max_lat')
parser.add_argument('-n', action='store', dest='n_subdivision',required=True,type=int,help='number of subdivision of icosahedron; n=7 -> 57.66 km, n=8 28.83 km')
parser.add_argument('-d', action='store', dest='depth_range',required=True,type=str,help='Depth range min_depth/max_depth in km')
# obsolete now values are calcultated automatically # JMN 29/05/2016
#parser.add_argument('-i', action='store', dest='step',required=True,type=float,help='Length (and width) of the dislocations in km; should be < -d')
parser.add_argument('--verbose', '-v', action='count',default=0,help='verbose mode')
parser.add_argument('--debug', action='count',default=0,help='debug mode')
#parser.add_argument('--netcdf', action='count',default=0,help='uses python netCDF4 module rather than GMT')
parser.add_argument('-e', action='store',type=str,required=True, dest='experiment',help='experiment name')

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
    print("-- Debug mode")
    verbose=True
    debug=True
else:
    debug=False
    
# netcdf
#if args.netcdf>0:
#    print("-- Using python netCDF4 module")
#    netcdf=True
#else:
#    netcdf=False

Rt=6371.0E3

    
#    return (netcdfgrid,grid_steps,experiment, bounds, depth_range,dwidth)
#(netcdfgrid,grid_steps,experiment, bounds, depth_range,dwidth)=get_arg()
#print "(netcdfgrid,grid_steps,experiment)",(netcdfgrid,grid_steps,experiment)

#dx=args.step

(min_depth,max_depth)=list(map(float,args.depth_range.split('/')))
if debug:print("-- Bounds (long. in/max lat. min/max): "," ".join(args.bounds.split('/')[1:]))
(lon_min,lon_max,lat_min,lat_max)=list(map(float,args.bounds.split('/')[1:]))
(rlon_min,rlon_max,rlat_min,rlat_max)=list(map(np.radians,(lon_min,lon_max,lat_min,lat_max)))


###################################################################
# SUBROUTINES
###################################################################

###################################################################
def geo_center_face(face,verts):
###################################################################
    
    """returns long. lat. in degrees of a face barycenter"""

    Rt=6371.0E3
    x_center=0
    y_center=0
    z_center=0
    for j in [0,1,2]:
        (x,y,z)=verts[face[j]]
        x=x*Rt
        y=y*Rt
        z=z*Rt
        x_center=x_center+x
        y_center=y_center+y
        z_center=z_center+z
    x_center=x_center/3.
    y_center=y_center/3.
    z_center=z_center/3.
    (l,phi,he)=coordinates.xyz2geo(x_center,y_center,z_center)
    return(np.degrees(l),np.degrees(phi),he)


###################################################################
def write_triangles(faces,verts,name):
###################################################################

    ftriangles=open(name,'w')

    Rt=6371.0E3
    for i in range(len(faces)):
        ftriangles.write(">\n")
        face=faces[i]
        for j in [0,1,2]:
            (x,y,z)=verts[face[j]]
            x=x*Rt
            y=y*Rt
            z=z*Rt
            (l,phi,_he)=coordinates.xyz2geo(x,y,z)
            #print ("%10.5lf %10.5lf \n" % (math.degrees(l),math.degrees(phi)))
            ftriangles.write("%10.5lf %10.5lf \n" % (np.degrees(l),np.degrees(phi)))
    ftriangles.close()
    

###################################################################
def get_depth(longitude,latitude,grid):
###################################################################
    ftmp_xy=open('tmp.xy','w')
    ftmp_xy.write("%10.5lf %10.5lf\n" % (longitude,latitude))
    ftmp_xy.close()
    import subprocess
    cmd="gmt grdtrack "+'tmp.xy'+ " -G"+grid+ " > tmpgrid.dat"
    #print "- running ", cmd
    subprocess.getstatusoutput(cmd)
    
    fs=open('tmpgrid.dat','r')
    lline = fs.readlines()
    if len(lline) <1:str_depth='NaN'   
    else:str_depth=lline[0].split()[-1]
    return(str_depth)

# ###################################################################
# def get_depth_netcdf(longitude,latitude,grid):
# ###################################################################
#     from mpl_toolkits import basemap
#     lat = grid.variables['y'][:]
#     lon = grid.variables['x'][:]
#     z = grid.variables['z'][:].data
#     return(basemap.interp(z, lon, lat, np.array([longitude]),np.array([latitude]), order=1))

###################################################################
def face_area_strike_dip(face,verts,grid):
###################################################################
    
    Rt=6371.0E3
    # first get XYZ coordinates of face's vertices on the Sphere
    A=np.array(verts[face[0]])
    B=np.array(verts[face[1]])
    C=np.array(verts[face[2]])
    A=A*Rt
    B=B*Rt
    C=C*Rt
    # Get depth
    (l,phi,he)=coordinates.xyz2geo(A[0],A[1],A[2])
    longitude=np.degrees(l)
    latitude=np.degrees(phi)
    depth=get_depth(longitude,latitude,args.grd)
    if debug:print(("    - vertex 1: %10.5lf %10.5lf %s" %(longitude,latitude,depth)))
    if depth=='NaN':return('NaN','NaN','NaN')
    #print longitude,latitude,depth
    rdepth=-np.sqrt(float(depth)**2)*1.E3 # depth in m & must be negative
    
    AA=np.array(coordinates.geo2xyz(l,phi,rdepth))

    (l,phi,he)=coordinates.xyz2geo(B[0],B[1],B[2])
    longitude=np.degrees(l)
    latitude=np.degrees(phi)
    depth=get_depth(longitude,latitude,args.grd)
    if debug:print(("    - vertex 2: %10.5lf %10.5lf %s" %(longitude,latitude,depth)))
    if depth=='NaN':return('NaN','NaN','NaN')
#    if debug:print 'depth ',depth 
    rdepth=-np.sqrt(float(depth)**2)*1.E3 # depth in m & must be negative
    #print longitude,latitude,depth
    BB=np.array(coordinates.geo2xyz(l,phi,rdepth))

    (l,phi,he)=coordinates.xyz2geo(C[0],C[1],C[2])
    longitude=np.degrees(l)
    latitude=np.degrees(phi)
    depth=get_depth(longitude,latitude,args.grd)
    if debug:print(("    - vertex 3: %10.5lf %10.5lf %s" %(longitude,latitude,depth)))
    if depth=='NaN':return('NaN','NaN','NaN')
#    if debug:print 'depth ',depth 
    rdepth=-np.sqrt(float(depth)**2)*1.E3 # depth in m & must be negative
    #print longitude,latitude,depth
    CC=np.array(coordinates.geo2xyz(l,phi,rdepth))
    
    # Normal vector
    
    N=np.cross((BB-AA),(CC-AA))
    length=np.sqrt(N[0]**2+N[1]**2+N[2]**2)
    n=N/length
    
    # area in km^2 S=0.5 * norm (B-A x C-A) 
    area=.5*length/1.E6
    # barycenter of the triangle
    M=(AA+BB+CC)/3
#    lnup= M/np.sqrt(M[0]**2+M[1]**2+M[2]**2)
#    dip=math.degrees(math.acos(np.dot(n,lnup)))

    # rotate to local frame in ENU    
    (l,phi,he)=coordinates.xyz2geo(M[0],M[1],M[2])
    R=coordinates.mat_rot_general_to_local(l,phi)
    ENU=np.dot(R,n)

    # we want the normal vector to be upward
    
    if ENU[2]<0:ENU=-ENU
    
    # dip

    dip=np.degrees(np.arctan( (np.sqrt(ENU[0]**2+ENU[1]**2)) / ENU[2]   ))
    
    # usv the unit strike vector

    usv=np.array([-ENU[1],ENU[0],0.])/np.sqrt(ENU[0]**2+ENU[1]**2)

    # get strike
    
    strike=np.degrees(np.arctan2(usv[0],usv[1]))
    
    if strike>90.0 or strike < -90.0: 
        if verbose:print('ENU strike dip lp ',ENU, strike, dip, np.degrees(l),np.degrees(phi))
        
    # For a vector normal to a plane, V=(East,North,Up)
    # strike = atan2(East,North)
    #strike=math.degrees(math.atan2(ENU[0],ENU[1]))
#     print '***** strike ',strike
#     if strike>180.:strike=strike-360.
#  
#     if strike>90.0:
#         print '!!!! strike > 90 ',strike
#         strike=strike-180.0
#         print '!!! new strike ',strike
#     if strike<-90.0:
#         print '!!!! strike < -90 ',strike
#         strike=strike+180.0
#         print '!!! new strike ',strike
#         
#     if dip>90.:
#         print '!!!! print dip > 90'
#         strike=strike+180.
#         if strike>180.:strike=strike-360.
#         dip=180.-dip
#         print '!!! new dip ',dip
    
    return(area,strike,dip)



###################################################################
# NOW CREATING THE GLOBAL MESH
# Add zone='global' for a mesh covering the whole Earth
###################################################################

(verts,faces)=Triangular_Mesh_Global.mesh_regional(num_subdivisions=args.n_subdivision, bounds=args.bounds)

###################################################################
# NOW SELECTING FACES ON THE GRID
###################################################################

print("-- Now selecting faces on the grid...")

lfaces=[]

print("-- Number of faces to be tested: ",len(faces))

#if netcdf:
#    pass
#     import netCDF4
#     grid = netCDF4.Dataset(args.grd)

for face in faces:
    (longitude,latitude,he)=geo_center_face(face,verts)
    if debug: print(("- Testing long. lat %10.4lf %10.4lf" % (longitude,latitude)))
#    if netcdf:
#        str_depth=get_depth_netcdf(longitude,latitude,grid)
#    else:

    str_depth=get_depth(longitude,latitude,args.grd)
    if str_depth == 'NaN' or str_depth == 'nan':
        if debug : print("         -Outside grid ; rejected")
        continue
    
    if debug : print("         -Inside grid ; kept")
    
    depth=float(str_depth)
    if np.fabs(depth)<max_depth and np.fabs(depth)>min_depth:lfaces.append(face)

faces=lfaces
name_triangles=args.experiment+'_triangles.dat'
write_triangles(faces,verts,name_triangles)

###################################################################
# CALCULATES FAULTS PARAMETERS FOR EACH SOURCE POINT
###################################################################

    
# select triangles with all vertices on grid
# writes dislocation, source file

# Dislocations sources files
from pyeq.lib.eq_disloc_3d import Dislocation as Dislocation
# one degree in km
one_degree=111.

f_sources=args.experiment+'_sources_lp.dat'
f_dislocations=args.experiment+'_dislocations.dat'
f_gmt=args.experiment+'_rectangular_dislocations.gmt'
f_geometry=args.experiment+'_geometry.dat'
f_geometry_npy=args.experiment+'_geometry.npy'

print("-- Now selecting triangles on grid, calculating area, strike and dip, and creating dislocation file ", f_dislocations)
print("-- Number of faces to be tested ",len(faces))

fsources=open(f_sources,'w')
fdislocations=open(f_dislocations,'w')
frectangular_gmt=open(f_gmt,'w')
fgeometry=open(f_geometry,'w')

fdislocations.write("#index longitude   latitude    depth   strike      dip   length    width         area    rake  max_slip\n")

index_source=0
lfaces=[]

# get n_dislocations
n_dis=0
larea=[]

lindex_face=[]


for i in range(len(faces)):
    face=faces[i]
    (area,strike,dip)=face_area_strike_dip(face,verts,args.grd)
    if area!='NaN':
        n_dis+=1
        larea.append(area)
        lindex_face.append(i)
        
median_area=np.median(larea)
print(("    -- median_area: %10.2f km**2" % median_area))
dislocation_length_km=np.sqrt(median_area)
print(("    -- length/width of dislocation: %10.2f km" % dislocation_length_km))

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
rdis_length=dislocation_length_km
rdis_width=dislocation_length_km
rdis_area=rdis_length*rdis_width

# added to get npy for meade block code

lidx_vert = []

for i in sorted(lindex_face):
    face=faces[i]
    print("tde #%4d: %04d %04d %04d " % (i,face[0],face[1],face[2]) )
    lidx_vert.append( face[0] )
    lidx_vert.append( face[1] )
    lidx_vert.append( face[2] )

# make the list uniq
np_idx_vert = np.array( list(set( lidx_vert ) ))
print( np_idx_vert )
# make the numpy array for save
np_v = np.zeros((len(lindex_face),3),dtype=int)
np_c = np.zeros( ( np_idx_vert.shape[0],3 ) )

# fills np_v
# idx_v=0
# for i in sorted(lindex_face):
#     face=faces[i]
#     
#     np_v[idx_v,0] = np.where( np_idx_vert == face[0] )[0] + 1
#     np_v[idx_v,1] = np.where( np_idx_vert == face[1] )[0] + 1
#     np_v[idx_v,2] = np.where( np_idx_vert == face[2] )[0] + 1
# 
#     idx_v = idx_v + 1
# np.savetxt('meade_v.dat',np_v,fmt="%05d  %05d  %05d")
# 
# # fills np_c
# 
# for i in np.arange( np_idx_vert.shape[0] ):
#     (x,y,z)=verts[ np_idx_vert[i] ]
#     x=x*Rt
#     y=y*Rt
#     z=z*Rt
#     (l,phi,he)=coordinates.xyz2geo(x,y,z)
#     longitude=np.degrees(l)
#     latitude=np.degrees(phi)
#     depth= - np.sqrt( float(get_depth(longitude,latitude,args.grd))**2 )
# 
#     np_c[i,:] = [longitude,latitude,depth]
# 
# np.savetxt('meade_c.dat',np_c,fmt="%10.5lf  %10.5lf  %10.5lf")



index_dis=0
for i in sorted(lindex_face):
    if verbose:print(("-- Face #%d over %d" % (i,n_dis)))
    face=faces[i]
    
    (area,strike,dip)=face_area_strike_dip(face,verts,args.grd)

    # Added for weird marine geophysics geometry
    if strike > 90.0: 
        print('!!!! strike dip',strike,dip)
        strike=strike-90.0
        print('!!!! strike corrected to ',strike)
    if strike < -90.0:
        print('!!!! strike dip',strike,dip)
        strike=strike+90.0 
        print('!!!! corrected to ',strike,dip)


    (x,y,z)=verts[face[0]]
    x=x*Rt
    y=y*Rt
    z=z*Rt
    (l,phi,he)=coordinates.xyz2geo(x,y,z)
    longitude=np.degrees(l)
    latitude=np.degrees(phi)
    depth=float(get_depth(longitude,latitude,args.grd))
    (tdis_long1,tdis_lat1,tdis_depth1)=(longitude,latitude,depth)

    (x,y,z)=verts[face[1]]
    x=x*Rt
    y=y*Rt
    z=z*Rt
    (l,phi,he)=coordinates.xyz2geo(x,y,z)
    longitude=np.degrees(l)
    latitude=np.degrees(phi)
    depth=float(get_depth(longitude,latitude,args.grd))
    (tdis_long2,tdis_lat2,tdis_depth2)=(longitude,latitude,depth)
    
    (x,y,z)=verts[face[2]]
    x=x*Rt
    y=y*Rt
    z=z*Rt
    (l,phi,he)=coordinates.xyz2geo(x,y,z)
    longitude=np.degrees(l)
    latitude=np.degrees(phi)
    depth=float(get_depth(longitude,latitude,args.grd))
    (tdis_long3,tdis_lat3,tdis_depth3)=(longitude,latitude,depth)
    
    area_tdis=area
    
    if debug:print(("OK: area,strike,dip: %10.2lf  %5.2lf %5.2lf " % (area,strike,dip)))
    lfaces.append(face)

    # WRITES SOURCE POINT LOCATION
    (lon,lat,he)=geo_center_face(face,verts)
    depth=float(get_depth(lon,lat,args.grd))
    
    centroid_long=lon
    centroid_lat=lat
    centroid_depth=depth
    
    
    X_SOURCE=np.array([lon,lat,depth])
    fsources.write("%05d %10.5lf %10.5lf\n" %(index_source,lon,lat))
    if debug:print('long lat depth ',lon,lat,depth)
    # WRITES DISLOCATIONS PARAMETERS
        # First, calculates origin of small rectangular dislocation having lon, lat as center
    dislocation=Dislocation(index=None, x=lon, y=lat, depth=depth, strike=strike, dip=dip, \
                                        length=rdis_length/one_degree, width=rdis_width/one_degree)
    (X1,X2,X3,X4)=dislocation.corners()
    
    CENTER_DISLOCATION=(X1+X2+X3+X4)/4
    DELTA=CENTER_DISLOCATION-X_SOURCE
    ORIGIN_DISLOCATION=X_SOURCE-DELTA
    
    (X1,X2,X3,X4)=dislocation.corners()
    X=X1-(X3-X1)/2.

    (od_x,od_y,od_z)=list(ORIGIN_DISLOCATION)
    fdislocations.write("%05d %10.5lf %10.5lf %8.1lf %8.1lf %8.2lf %8.2lf %8.2lf %12.2lf\n" \
                        %(index_source,od_x,od_y,od_z,strike,dip,rdis_length,rdis_width, area))
    if debug:print(("source #%05d %10.5lf %10.5lf %8.1lf %8.1lf %8.2lf %8.2lf %8.2lf" %(index_source,X1[0],X1[1],X1[2],strike,dip,rdis_length,rdis_width)))
    dislocation=Dislocation(index=index_source, x=od_x, y=od_y, depth=od_z, strike=strike, dip=dip, \
                                        length=rdis_length, width=rdis_width)
    
    rdis_long =od_x
    rdis_lat  =od_y
    rdis_depth=od_z
    
    # writing gmt psxy rectangular dislocation geometry file 
    frectangular_gmt.write("#%05d\n" % index_source)
    frectangular_gmt.write(">\n")
    frectangular_gmt.write("%10.5lf %10.5lf \n" % ((X1-DELTA)[0],(X1-DELTA)[1]))
    frectangular_gmt.write("%10.5lf %10.5lf \n" % ((X2-DELTA)[0],(X2-DELTA)[1]))
    frectangular_gmt.write("%10.5lf %10.5lf \n" % ((X3-DELTA)[0],(X3-DELTA)[1]))
    frectangular_gmt.write("%10.5lf %10.5lf \n" % ((X4-DELTA)[0],(X4-DELTA)[1]))
    frectangular_gmt.write(">\n")

    # GEOMETRY
    
    tdis_area=area
    ratio_rdis_tdis=tdis_area / rdis_area
    
    GEOMETRY[index_source,:]=[rdis_long,rdis_lat,rdis_depth,\
                              rdis_length,rdis_width,rdis_area,ratio_rdis_tdis,\
                              strike,dip,\
                              centroid_long,centroid_lat,centroid_depth,\
                              tdis_long1, tdis_lat1, tdis_depth1, \
                              tdis_long2, tdis_lat2, tdis_depth2, \
                              tdis_long3, tdis_lat3, tdis_depth3, \
                              tdis_area]

    index_source=index_source+1
        
# print results


faces=lfaces
name_triangles=args.experiment+'_triangles.dat'
write_triangles(faces,verts,name_triangles)

print("-- files for gmt plots ", f_sources, f_dislocations, f_gmt)

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


fsources.close()
fdislocations.close()
frectangular_gmt.close()

# print gmt and shapefiles

from pyeq.log.geometry2shp_gmt import geometry2shp_gmt

geometry2shp_gmt( GEOMETRY, 'tde' , out_shp=args.experiment + '_tde', out_gmt= args.experiment + '_tde.gmt' , verbose=verbose) 
geometry2shp_gmt( GEOMETRY, 'rde' , out_shp=args.experiment + '_rde', out_gmt= args.experiment + '_rde.gmt' , verbose=verbose) 

# test save verts and faces

#np.savetxt('faces.dat',np.array(faces))
#np.savetxt('verts.dat',np.array(verts))


