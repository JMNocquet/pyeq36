#!/usr/bin/env python

###################################################################
# SCRIPT    : pyeq_make_green.py
# AUTHOR    : 
# DATE      : December 2011 - Updated Fall 2017
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

###################################################################
# PYACS ADDITIONAL MODULES IMPORT
###################################################################

import pyacs.lib.coordinates
from pyacs.lib import gmtpoint as GMT_Point
from pyeq.lib import eq_disloc_3d as Dislocation
from pyeq.lib import lib_inversion

###################################################################
# PARSE ARGUMENT LINE
###################################################################

prog_info="pyeq_make_green.py : calculates the transfert (Green) functions relating surface deformation observations to slip on a list of individual faults\n"
prog_info+="                     this version only handles a fixed rake (although individual faults can have different fixed rake)\n"
prog_info+="                     rake can be provided either as a value or a relative motion given by an Euler pole (with option --pole)\n"
prog_info+="                     with the --pole option be careful of the sign of the rake calculated (double check and if not OK , change the sign of the angular velocity)\n"
prog_info+="                     the --pole option will also indicate the maximum velocity (for constrained inversion, useful for interseismic)\n"
prog_info+="                     check options of the inversion program\n"
prog_info+="                     Output: matrices in numpy npy format - name starts with experiment name\n";
prog_info+="-----------------------------------------------------------\n"
prog_info+="FORMAT FOR DISLOCATION FILE\n"
prog_info+="-----------------------------------------------------------\n"
prog_info+="\%05d \%10.5lf \%10.5lf \%8.1lf \%8.1lf \%8.2lf \%8.2lf \%8.2lf\n" 
prog_info+=" index long lat depth strike dip length width\n"

prog_epilog="J.-M. Nocquet (Geoazur-CNRS) - May 2012"

parser = argparse.ArgumentParser(description=prog_info,epilog=prog_epilog)
parser.add_argument('-gps_h', action='store', dest='gmth',default=None,help='gmt psvelo files for horizontal displacement/velocity (mm or mm/y)')
parser.add_argument('-gps_u', action='store', dest='gmtu',default=None,help='files for vertical displacement/velocity - positive upwards with 4 columns: lon, lat, v_up, s_v_up (mm or mm/y)')
parser.add_argument('-insar', action='store', dest='insar',default=None,help='files for insar data. Format: Number xind yind east north data err wgt Elos Nlos Ulos (in metres)')
parser.add_argument('-g', action='store', dest='geometry',required=True,help='geometry file including dislocations and sources in npy format; see help for format')
parser.add_argument('--tde', action='count',default=0,help='triangular dislocations computed using Meade (2007)')
parser.add_argument('--verbose', '-v', action='count',default=0,help='verbose mode')
parser.add_argument('-e', action='store',type=str,required=True, dest='experiment',help='experiment name')

args = parser.parse_args()

if (len(sys.argv)<2):parser.print_help();sys.exit()

#if not (args.pole or args.rake):
#    print "=> ERROR : either --pole or --rake should be provided"
#    sys.exit()
    

# verbose
if args.verbose>0:
    print("-- Verbose mode")
    verbose=True
else:
    verbose=False

# TDE
if args.tde>0:
    from pyeq.lib.meade_tde import tde
    Poisson_Ratio=0.25
    TDE=True
else:
    TDE=False

###############################################################################
# H_FAULTS INITALIZATION
###############################################################################

H_fault_lp={}
H_fault_xy={}

###############################################################################
# READS GEOMETRY FILE AS NPY
###############################################################################

import pyeq.lib.geometry.to_np_array


if args.geometry[-3:] == 'dat':
    # reads the dat text file
    GEOMETRY , SGEOMETRY = pyeq.lib.geometry.to_np_array.dat_geometry_to_array_and_recarray( args.geometry , verbose=verbose )

if args.geometry[-3:] == 'npy':
    # reads the npy
    GEOMETRY , SGEOMETRY = pyeq.lib.geometry.to_np_array.npy_geometry_to_array_and_recarray( args.geometry , verbose=verbose )



print("-- Building dislocations from ",args.geometry)

# TRIANGULAR DISLOCATONS
if TDE:
    X=np.zeros(3)
    Y=np.zeros(3)
    Z=np.zeros(3)
    
    XYZ=np.zeros((3,3))
    
    H_TDE={}
    print("-- Triangular dislocations will be used")

for i in range(SGEOMETRY.shape[0]):

    if verbose:
        print('  -- ',i,' / ',SGEOMETRY.shape[0])
        

    [rdis_long,rdis_lat,rdis_depth,rdis_length,rdis_width,rdis_area,ratio_rdis_tdis,strike,dip,\
     centroid_long,centroid_lat,centroid_depth,\
     tdis_long1,tdis_lat1,tdis_depth1,tdis_long2,tdis_lat2,tdis_depth2,tdis_long3,tdis_lat3,tdis_depth3,\
     tdis_area]=np.array(list(SGEOMETRY[i]))

    # triangular dislocation
    if TDE:
        (X[0],Y[0])=pyacs.lib.coordinates.geo2flat_earth(tdis_long1,tdis_lat1);Z[0]=tdis_depth1
        (X[1],Y[1])=pyacs.lib.coordinates.geo2flat_earth(tdis_long2,tdis_lat2);Z[1]=tdis_depth2
        (X[2],Y[2])=pyacs.lib.coordinates.geo2flat_earth(tdis_long3,tdis_lat3);Z[2]=tdis_depth3
    
    
        XYZ[:,0]=X
        XYZ[:,1]=Y
        XYZ[:,2]=Z
    
        H_TDE[i]=np.copy(XYZ)
        
    # rectangular dislocations
    else:

        depth=np.sqrt(rdis_depth**2)
        index_fault=i
        if (rdis_long > 180.):rdis_long=rdis_long-360.
        
        lon=rdis_long
        lat=rdis_lat
        length=rdis_length
        width=rdis_width
        area=rdis_area
        
        # fake values
        
        rake=0.0
        max_slip=0.0
        
        (x,y) = pyacs.lib.coordinates.geo2flat_earth(lon,lat)
        dislocation_lp=Dislocation.Dislocation(index_fault, lon, lat, depth, strike, dip, length, width,area, rake, max_slip)
        dislocation_xy=Dislocation.Dislocation(index_fault, x, y, depth, strike, dip, length, width,area, rake, max_slip)
        H_fault_xy[index_fault]=dislocation_xy
        H_fault_lp[index_fault]=dislocation_lp
         
###############################################################################
# READS GMT FILE FOR HORIZONTAL COMPONENTS
###############################################################################


print("-- Reading ",args.gmth)


OBS=np.genfromtxt(args.gmth,comments='#',usecols=(list(range(7))))
NAME_OBS=np.genfromtxt(args.gmth,comments='#',usecols=(7),dtype=str)

# sorting
lindex = np.argsort( NAME_OBS )
OBS = OBS[ lindex ]
NAME_OBS = NAME_OBS[ lindex ]



array_gps=np.zeros((OBS.shape[0],2))

if TDE:
    SX=np.zeros(OBS.shape[0])
    SY=np.zeros(OBS.shape[0])
    SZ=np.zeros(OBS.shape[0])

print("-- Number of GPS sites (horizontal components) :", OBS.shape[0])

# 
index_gmt_point=0
H_GMT_Point_xy={}
H_GMT_Point_lp={}

[lon, lat] = [OBS[:,0], OBS[:,1]]
(x,y) = pyacs.lib.coordinates.geo2flat_earth( lon, lat )
array_gps = np.array([x,y]).T

for i in range(OBS.shape[0]):
    [lon, lat, ve, vn, sve, svn, sven]=OBS[i,:]
    code=NAME_OBS[i]
    GPS_Point_lp=GMT_Point.GMT_Point(code=code,lon=lon,lat=lat,he=0., Ve=ve,Vn=vn,Vu=0,SVe=sve,SVn=svn,SVu=0,SVen=sven,Cv_xyz=None, Cv_enu=None, index=index_gmt_point)
    if TDE:
        SX[i]=x[i]
        SY[i]=y[i]
    GPS_Point_xy=GMT_Point.GMT_Point(code=code,lon=x,lat=y,he=0., Ve=ve,Vn=vn,Vu=0,SVe=sve,SVn=svn,SVu=0,SVen=sven,Cv_xyz=None, Cv_enu=None, index=index_gmt_point)
    H_GMT_Point_xy[i]=GPS_Point_xy
    H_GMT_Point_lp[i]=GPS_Point_lp


###############################################################################
# READS UP FILE IF PROVIDED
###############################################################################

OBS_UP=np.zeros((0,0))
NAME_OBS_UP=np.zeros((0))

if args.gmtu != None:
    
    print("-- Reading ",args.gmtu)

    OBS_UP=np.genfromtxt(args.gmtu,comments='#',usecols=(list(range(4))))
    NAME_OBS_UP=np.genfromtxt(args.gmtu,comments='#',usecols=(4),dtype=str)

    array_gps_up=np.zeros((OBS_UP.shape[0],2))

    if TDE:
        SX_U=np.zeros(OBS.shape[0])
        SY_U=np.zeros(OBS.shape[0])
        SZ_U=np.zeros(OBS.shape[0])


    print("-- Number of GPS sites (vertical component) :", OBS_UP.shape[0])
    
    index_up_point=0
    H_GMT_Point_xy_up={}
    H_GMT_Point_lp_up={}
    
    for i in range(OBS_UP.shape[0]):
        [lon, lat, vu, svu]=OBS_UP[i,:]
        code=NAME_OBS_UP[i]
        GPS_Point_lp_up=GMT_Point.GMT_Point(code=code,lon=lon,lat=lat,he=None, Ve=None,Vn=None,Vu=vu,SVe=None,SVn=None,SVu=svu,SVen=None,Cv_xyz=None, Cv_enu=None, index=index_gmt_point)
        (x,y)=pyacs.lib.coordinates.geo2flat_earth(lon,lat)
        array_gps_up[i,:]=[x,y]
        if TDE:
            SX_U[i]=x
            SY_U[i]=y
        GPS_Point_xy_up=GMT_Point.GMT_Point(code=code,lon=x,lat=y,he=0., Ve=ve,Vn=vn,Vu=vu,SVe=None,SVn=None,SVu=svu,SVen=None,Cv_xyz=None, Cv_enu=None, index=index_gmt_point)
        H_GMT_Point_xy_up[index_up_point]=GPS_Point_xy_up
        H_GMT_Point_lp_up[index_up_point]=GPS_Point_lp_up
        index_up_point=index_up_point+1

###############################################################################
# READS INSAR DATA FILE IF PROVIDED
###############################################################################

OBS_INSAR=np.zeros((0,0))

if args.insar != None:
    
    print("-- Reading ",args.insar)

    # data format must be index xind yind east north data err wgt Elos Nlos Ulos
    OBS_INSAR = np.genfromtxt(args.insar,comments='#',usecols = ( (3,4,5,6,8,9,10) ))
    OBS_INSAR[:,2] = OBS_INSAR[:,2]*1.E3
    OBS_INSAR[:,3] = OBS_INSAR[:,3]*1.E3
    # INSAR is now lon lat los sigma Elos Nlos Ulos with los in mm

    array_insar=np.zeros((OBS_INSAR.shape[0],2))

    if TDE:
        SX_INSAR=np.zeros(OBS_INSAR.shape[0])
        SY_INSAR=np.zeros(OBS_INSAR.shape[0])
        SZ_INSAR=np.zeros(OBS_INSAR.shape[0])


    print("-- Number of InSAR points :", OBS_INSAR.shape[0])
    
#    index_insar_point=0
#    H_GMT_Point_xy_insar={}
#    H_GMT_Point_lp_insar={}

    LON_LAT_LOS_SLOS = OBS_INSAR[:,:4]
    (X,Y)=pyacs.lib.coordinates.geo2flat_earth(LON_LAT_LOS_SLOS[:,0],LON_LAT_LOS_SLOS[:,1])
    
    array_insar = np.vstack((X,Y)).T

    if TDE:
        SX_INSAR = X
        SY_INSAR = Y

###############################################################################
# CREATES GREEN FUNCTIONS FOR HORIZONTAL COMPONENTS
###############################################################################

n_dislocations=SGEOMETRY.shape[0]
n_gps=OBS.shape[0]
slip=1.0

print("-- Creating Green's functions matrix for horizontal components")

# GREEN IS A TENSOR OF DIM 4
# GREEN(i,j,k,l) is the prediction for dislocation i at site j component k for rake l
# k=0,1,2 = east, north, up
# l=0,1 : rake_00 & rake_90

GREEN=np.zeros((n_dislocations, n_gps, 3,2))

if TDE:
    # observation points
    SX=array_gps[:,0]*1.E3
    SY=array_gps[:,1]*1.E3
    SZ=array_gps[:,0]*0.0

    for index in range(SGEOMETRY.shape[0]):
        
        if verbose:
            print('  -- ',index,' / ',SGEOMETRY.shape[0])
        
        X=H_TDE[index][:,0]*1.E3
        Y=H_TDE[index][:,1]*1.E3
        Z=H_TDE[index][:,2]*1.E3
        Z=np.sqrt(Z**2)
        
        
        U_rake_00=tde.calc_tri_displacements(SX,SY,SZ, X, Y, Z, Poisson_Ratio, -1.0, 0.0,  0.0)
        U_rake_90=tde.calc_tri_displacements(SX,SY,SZ, X, Y, Z, Poisson_Ratio,  0.0, 0.0, -1.0)
        
        green_rake_00=np.zeros((array_gps.shape[0],5))
        green_rake_90=np.zeros((array_gps.shape[0],5))
        
        green_rake_00[:,0:2]=array_gps
        green_rake_90[:,0:2]=array_gps

        green_rake_00[:,2]=U_rake_00['x']
        green_rake_00[:,3]=U_rake_00['y']
        # Meade tde convention: positive Uz downward - corrected 18/02/2020 by adding -
        green_rake_00[:,4]=-U_rake_00['z']

        green_rake_90[:,2]=U_rake_90['x']
        green_rake_90[:,3]=U_rake_90['y']
        # Meade tde convention: positive Uz downward - corrected 18/02/2020 by adding -
        green_rake_90[:,4]=-U_rake_90['z']

        GREEN[index,:,:,0]=green_rake_00[:,2:5]
        GREEN[index,:,:,1]=green_rake_90[:,2:5]


else:

    slip = 1.0

    for index in H_fault_xy.keys():
        fault_xy=H_fault_xy[index]
        fault_lp=H_fault_lp[index]
        
        
    #    ARRAY_SOURCES[index,:]=SOURCES[index,:]
    
    #     # add 24/03/2016 for the weird geometry associated with the seamount from JY. Collot (JGR, 2017)
    #     print '!!!! JY',index,fault_lp.x,fault_lp.y,fault_lp.strike,rake
    # #         
    #     if fault_lp.strike > 90.0: 
    #         print '!!!! ',index,fault_lp.x,fault_lp.y,fault_lp.strike,rake
    #         rake=rake+180.0 
    #         print '!!!! corrected to ',index,fault_lp.x,fault_lp.y,fault_lp.strike,rake
    #     if fault_lp.strike < -90.0:
    #         print '!!!! ',index,fault_lp.x,fault_lp.y,fault_lp.strike,rake
    #         rake=rake-180.0 
    #         print '!!!! corrected to ',index,fault_lp.x,fault_lp.y,fault_lp.strike,rake
    
            
        green_rake_00=fault_xy.disp_slip_rake(slip,0.0,array_gps)
        green_rake_90=fault_xy.disp_slip_rake(slip,90.0,array_gps)

#        green_rake_00=fault_xy.disp_slip_rake_no_edcmp(slip,0.0,array_gps)
#        green_rake_90=fault_xy.disp_slip_rake_no_edcmp(slip,90.0,array_gps)

        
        green_en=green_rake_00[:,2:4]
        
    #    print 'GREEN[index,:,:,0] ',GREEN[index,:,:,0].shape
    #    print 'green_rake_00[:,2:5] ',green_rake_00[:,2:5].shape
        GREEN[index,:,:,0]=green_rake_00[:,2:5]
        GREEN[index,:,:,1]=green_rake_90[:,2:5]

###############################################################################
# CREATES OBSERVATION FILE FOR VERTICAL COMPONENT IF PROVIDED
if OBS_UP.shape[0] > 0:
###############################################################################

    print("-- Creating Green's function matrix for the vertical component")

    n_gps_up=OBS_UP.shape[0]

    GREEN_UP=np.zeros((n_dislocations, n_gps_up, 3,2))


    if TDE:
        # observation points
        SX=array_gps_up[:,0]*1.E3
        SY=array_gps_up[:,1]*1.E3
        SZ=array_gps_up[:,0]*0.0

        for index in range(SGEOMETRY.shape[0]):

            X=H_TDE[index][:,0]*1.E3
            Y=H_TDE[index][:,1]*1.E3
            Z=H_TDE[index][:,2]*1.E3
            Z=np.sqrt(Z**2)
        
            U_rake_00=tde.calc_tri_displacements(SX,SY,SZ, X, Y, Z, Poisson_Ratio, -1.0, 0.0,  0.0)
            U_rake_90=tde.calc_tri_displacements(SX,SY,SZ, X, Y, Z, Poisson_Ratio,  0.0, 0.0, -1.0)


            green_rake_00=np.zeros((array_gps_up.shape[0],5))
            green_rake_90=np.zeros((array_gps_up.shape[0],5))
            
            green_rake_00[:,0:2]=array_gps_up
            green_rake_90[:,0:2]=array_gps_up
    
    
            green_rake_00[:,2]=U_rake_00['x']
            green_rake_00[:,3]=U_rake_00['y']
            # Meade tde convention: positive Uz downward - corrected 02/05/2018 by adding -
            green_rake_00[:,4]=-U_rake_00['z']
    
            green_rake_90[:,2]=U_rake_90['x']
            green_rake_90[:,3]=U_rake_90['y']
            # Meade tde convention: positive Uz downward - corrected 02/05/2018 by adding -
            green_rake_90[:,4]=-U_rake_90['z']
    
            GREEN_UP[index,:,:,0]=green_rake_00[:,2:5]
            GREEN_UP[index,:,:,1]=green_rake_90[:,2:5]


    else:

        for index in H_fault_xy.keys():
            fault_xy=H_fault_xy[index]
            fault_lp=H_fault_lp[index]
    #        ARRAY_SOURCES[index,:]=SOURCES[index,:]
    #         # add 24/03/2016 for the weird geometry associated with the seamount from JY.Y. Collot
    #         print '!!!! JY',index,fault_lp.x,fault_lp.y,fault_lp.strike,rake
    # #         
    #         if fault_lp.strike > 90.0: 
    #             print '!!!! ',index,fault_lp.x,fault_lp.y,fault_lp.strike,rake
    #             rake=rake+180.0 
    #             print '!!!! corrected to ',index,fault_lp.x,fault_lp.y,fault_lp.strike,rake
    #         if fault_lp.strike < -90.0:
    #             print '!!!! ',index,fault_lp.x,fault_lp.y,fault_lp.strike,rake
    #             rake=rake-180.0 
    #             print '!!!! corrected to ',index,fault_lp.x,fault_lp.y,fault_lp.strike,rake
                
            green_rake_00=fault_xy.disp_slip_rake(slip,0.0,array_gps_up)
            green_rake_90=fault_xy.disp_slip_rake(slip,90.0,array_gps_up)

#            green_rake_00=fault_xy.disp_slip_rake_no_edcmp(slip,0.0,array_gps_up)
#            green_rake_90=fault_xy.disp_slip_rake_no_edcmp(slip,90.0,array_gps_up)

            
            green_en=green_rake_00[:,2:4]
            GREEN_UP[index,:,:,0]=green_rake_00[:,2:5]
            GREEN_UP[index,:,:,1]=green_rake_90[:,2:5]
    
else:
        GREEN_UP=np.zeros(0)    
    

###############################################################################
# CREATES GREEN TENSOR FOR INSAR IF PROVIDED
if OBS_INSAR.shape[0] > 0:
###############################################################################

    print("-- Creating Green's function matrix for the InSAR data ",args.insar)

    n_insar=OBS_INSAR.shape[0]

    GREEN_INSAR=np.zeros((n_dislocations, n_insar, 3,2))


    if TDE:
        # observation points
        SX=array_insar[:,0]*1.E3
        SY=array_insar[:,1]*1.E3
        SZ=array_insar[:,0]*0.0

        for index in range(SGEOMETRY.shape[0]):

            X=H_TDE[index][:,0]*1.E3
            Y=H_TDE[index][:,1]*1.E3
            Z=H_TDE[index][:,2]*1.E3
            Z=np.sqrt(Z**2)
        
            U_rake_00=tde.calc_tri_displacements(SX,SY,SZ, X, Y, Z, Poisson_Ratio, -1.0, 0.0,  0.0)
            U_rake_90=tde.calc_tri_displacements(SX,SY,SZ, X, Y, Z, Poisson_Ratio,  0.0, 0.0, -1.0)


            green_rake_00=np.zeros((array_insar.shape[0],5))
            green_rake_90=np.zeros((array_insar.shape[0],5))
            
            green_rake_00[:,0:2]=array_insar
            green_rake_90[:,0:2]=array_insar
    
    
            green_rake_00[:,2]=U_rake_00['x']
            green_rake_00[:,3]=U_rake_00['y']
            # Meade tde convention: positive Uz downward - corrected 02/05/2018 by adding -
            green_rake_00[:,4]=-U_rake_00['z']
    
            green_rake_90[:,2]=U_rake_90['x']
            green_rake_90[:,3]=U_rake_90['y']
            # Meade tde convention: positive Uz downward - corrected 02/05/2018 by adding -
            green_rake_90[:,4]=-U_rake_90['z']
    
            GREEN_INSAR[index,:,:,0]=green_rake_00[:,2:5]
            GREEN_INSAR[index,:,:,1]=green_rake_90[:,2:5]


    else:

        for index in H_fault_xy.keys():
            fault_xy=H_fault_xy[index]
            fault_lp=H_fault_lp[index]
    #        ARRAY_SOURCES[index,:]=SOURCES[index,:]
    #         # add 24/03/2016 for the weird geometry associated with the seamount from JY.Y. Collot
    #         print '!!!! JY',index,fault_lp.x,fault_lp.y,fault_lp.strike,rake
    # #         
    #         if fault_lp.strike > 90.0: 
    #             print '!!!! ',index,fault_lp.x,fault_lp.y,fault_lp.strike,rake
    #             rake=rake+180.0 
    #             print '!!!! corrected to ',index,fault_lp.x,fault_lp.y,fault_lp.strike,rake
    #         if fault_lp.strike < -90.0:
    #             print '!!!! ',index,fault_lp.x,fault_lp.y,fault_lp.strike,rake
    #             rake=rake-180.0 
    #             print '!!!! corrected to ',index,fault_lp.x,fault_lp.y,fault_lp.strike,rake
            
#            print('-- edcmp')
            t0=time()
            green_rake_00=fault_xy.disp_slip_rake(slip,0.0,array_insar)
            green_rake_90=fault_xy.disp_slip_rake(slip,90.0,array_insar)
#            print('-- ', time()-t0)
#            print('-- pyacs_okada')
#            green_rake_00=fault_xy.disp_slip_rake_no_edcmp(slip,0.0,array_insar)
#            green_rake_90=fault_xy.disp_slip_rake_no_edcmp(slip,90.0,array_insar)
#            print('-- ', time()-t0)

            #green_en=green_rake_00[:,2:4]
            GREEN_INSAR[index,:,:,0]=green_rake_00[:,2:5]
            GREEN_INSAR[index,:,:,1]=green_rake_90[:,2:5]
else:
        GREEN_INSAR=np.zeros(0)    


###############################################################################
# MAKE THE DISTANCE MATRIX Dm
# THE DISTANCE MATRIX HAS SIZE mxm (m = number of model parameters)
# Dm is useful to calculate Cm
###############################################################################

print("-- Now calculating the distance matrices")

nfaults=SGEOMETRY.shape[0]
Dm=np.zeros((nfaults,nfaults))



for i in range(nfaults):
    if verbose:print(("   - fault #%05d over %05d" %(i,nfaults)))

    (long_ref,lat_ref,depth_ref)=(SGEOMETRY.centroid_long[i],SGEOMETRY.centroid_lat[i],SGEOMETRY.centroid_depth[i])
    (x_ref,y_ref,z_ref) = pyacs.lib.coordinates.geo2xyz(np.radians(long_ref),np.radians(lat_ref),-np.fabs(depth_ref))
    
    for j in range(i):
        #print i,j
        (long_current,lat_current,depth_current)=(SGEOMETRY.centroid_long[j],SGEOMETRY.centroid_lat[j],SGEOMETRY.centroid_depth[j])
        (x_current,y_current,z_current)= pyacs.lib.coordinates.geo2xyz(np.radians(long_current),np.radians(lat_current),-np.fabs(depth_current))
    
        distance_meters=np.sqrt(np.sum( (x_ref-x_current)**2 + (y_ref-y_current)**2 + (z_ref-z_current)**2 ))
        distance_km=distance_meters / 1000.0
        Dm[i,j]=distance_km
        Dm[j,i]=distance_km

        #print i,j,distance_km

Dm_name=args.experiment+'_Dm.npy'

print('-- Creating the relative Distance matrix ',Dm_name)
np.save(Dm_name,Dm)

###############################################################################
# MAKE THE DISTANCE MATRIX D0 relative to the barycenter
# D0 is useful to calculate sigma
###############################################################################

long0=np.mean(SGEOMETRY.centroid_long)
lat0=np.mean(SGEOMETRY.centroid_lat)
depth0=np.mean(SGEOMETRY.centroid_depth)

(x0,y0,z0)= pyacs.lib.coordinates.geo2xyz(np.radians(long0),np.radians(lat0),-np.fabs(depth0))

D0=np.zeros(nfaults)
for i in range(nfaults):
    if verbose:print(("   - fault #%05d over %05d" %(i,nfaults)))

    (long_ref,lat_ref,depth_ref)=(SGEOMETRY.centroid_long[i],SGEOMETRY.centroid_lat[i],SGEOMETRY.centroid_depth[i])
    (x_ref,y_ref,z_ref)= pyacs.lib.coordinates.geo2xyz(np.radians(long_ref),np.radians(lat_ref),-np.fabs(depth_ref))

    distance_meters=np.sqrt(np.sum( (x_ref-x0)**2 + (y_ref-y0)**2 + (z_ref-z0)**2 ))
    distance_km=distance_meters / 1000.0

    D0[i]=distance_km

D0_name=args.experiment+'_D0.npy'

print('-- Creating the absolute Distance matrix ',D0_name)
np.save(D0_name,D0)


######################################################################################
# SAVE A BIG STRUCTURE ARRAY
# ORDER IS: GEOMETRY, Dm, ARRAY_SOURCES, DISLOCATIONS, GREEN, GREEN_UP,MAX_SLIP,OBS,NAME_OBS,OBS_UP,NAME_OBS_UP
######################################################################################


NPZ_OUT = ("%s_input.npz" % args.experiment)

print('-- Saving the inversion structured array %s '% NPZ_OUT)

# PYHTON 2.7
#experiment=args.experiment
# f = file(experiment+'_input.npy',"wb")
# 
# np.save(f,SGEOMETRY)
# np.save(f,Dm)
# np.save(f,GREEN)
# np.save(f,GREEN_UP)
# np.save(f,OBS)
# np.save(f,NAME_OBS)
# np.save(f,OBS_UP)
# np.save(f,NAME_OBS_UP)
# np.save(f,GREEN_INSAR)
# np.save(f,OBS_INSAR)
# f.close()
# Can be read sequentially
# np.save(f,b)
# np.save(f,c)
# f = file("tmp.bin","rb")
# aa = np.load(f)
# bb = np.load(f)
# cc = np.load(f)
# f.close()

np.save('GREEN.npy',GREEN)
np.save('GREEN_UP.npy',GREEN_UP)
np.save('GREEN_INSAR.npy',GREEN_INSAR)

# PYHTON 3.6
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

# check NPZ

npz = np.load(NPZ_OUT)

for array_name, array_value in npz.items():
    print("-- shape %s %s" % (array_name,str(array_value.shape)))

