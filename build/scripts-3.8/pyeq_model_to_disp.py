#!/usr/local/geodesy/anaconda38/bin/python

###################################################################
# SCRIPT    : pyeq_model_to_disp.py
# AUTHOR    : 
# DATE      : July 2018
# OUTPUT    : 
# NOTE      :
#         
###################################################################


###################################################################
# PYTHON MODULES IMPORT
###################################################################

import argparse, sys
import numpy as np
import math

###################################################################
# PYACS ADDITIONAL MODULES IMPORT
###################################################################

from pyacs.lib import coordinates
from pyeq.lib import eq_disloc_3d as Dislocation
import pyacs.lib.faultslip
import pyacs.lib.utils


###################################################################
# PARSE ARGUMENT LINE
###################################################################

prog_info="pyeq_model_to_disp.py : computes displacements from slip model\n"
prog_epilog="J.-M. Nocquet (Geoazur-CNRS) - May 2012"

parser = argparse.ArgumentParser(description=prog_info,epilog=prog_epilog)
parser.add_argument('-i', action='store', dest='input',required=True,help='file  with lon lat (dec.deg) where displacement will be calculated')
parser.add_argument('-model', action='store', dest='model',required=True,help='slip model')
parser.add_argument('--tde', action='count',default=0,help='triangular dislocations computed using Meade (2007)')
parser.add_argument('--verbose', '-v', action='count', default=0, help='verbose mode')
parser.add_argument('-o', action='store',type=str,required=True, dest='output',help='output file name')

args = parser.parse_args()


if (len(sys.argv)<2):parser.print_help();sys.exit()

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
# READS INPUT LOCATION FILE
###############################################################################


print("-- Reading ",args.input)

try:
    OBS = np.array(np.asmatrix(np.genfromtxt(args.input))[:,:2])
except:
    print('!!! Could not read ',args.input)

try:
    NAMES = np.genfromtxt(args.input, usecols=(7) , dtype=str )
except:
    print("[PYACS WARNING no name found in :  %s" % args.input)
    NAMES=np.arange(OBS.shape[0])
    
array_gps=np.zeros((OBS.shape[0],2))

if TDE:
    SX=np.zeros(OBS.shape[0])
    SY=np.zeros(OBS.shape[0])
    SZ=np.zeros(OBS.shape[0])

print("-- Number of points :", OBS.shape[0])

lon = OBS[:,0]
lat = OBS[:,1]
(x,y)=coordinates.geo2flat_earth(lon,lat)
array_gps = np.array( np.asmatrix(np.array( [x,y] ).T) )

if TDE:
    SX = x * 1.E3
    SY = y * 1.E3
    SZ = SY * 0.0
    
###############################################################################
# H_FAULTS INITALIZATION
###############################################################################

H_fault_lp={}
H_fault_xy={}

###############################################################################
# READS MODEL
###############################################################################

print("-- reading model ",args.model)

try:
    MODEL =  np.array(np.asmatrix(np.genfromtxt(args.model)))
    
except:
    print("!!! Could not read ",args.model)
    sys.exit()

# print some informations about the geometry
print(("  -- geometry includes %04d subfaults" % MODEL.shape[0]))
min_lon=np.min(MODEL[:,9])
max_lon=np.max(MODEL[:,9])

min_lat=np.min(MODEL[:,10])
max_lat=np.max(MODEL[:,10])

min_depth=np.min(MODEL[:,11])
max_depth=np.max(MODEL[:,11])

print(("  -- geometry bounds read from centroids %.5lf/%.5lf/%.5lf/%.5lf and depth %.5lf/%.5lf " % (min_lon,max_lon,min_lat,max_lat,min_depth,max_depth)))

RAKE_MAIN = MODEL[:,22]
RAKE_PERP = MODEL[:,24]

# removes elements with 0 slip

lindex_no_null_slip = np.where( MODEL[:,-1]==0. )
print("--" , lindex_no_null_slip[0].size, " subfaults with zero slip. Removing them.")
MM = np.delete( MODEL, lindex_no_null_slip, axis=0 )

print("-- Building dislocations from ",args.model)

###############################################################################
# TRIANGULAR DISLOCATONS
###############################################################################

if TDE:

   
    ( X_t_dis1 , Y_t_dis1 ) = coordinates.geo2flat_earth(MM[:,12],MM[:,13]); Z_t_dis1 = MM[:,14]
    ( X_t_dis2 , Y_t_dis2 ) = coordinates.geo2flat_earth(MM[:,15],MM[:,16]); Z_t_dis2 = MM[:,17]
    ( X_t_dis3 , Y_t_dis3 ) = coordinates.geo2flat_earth(MM[:,18],MM[:,19]); Z_t_dis3 = MM[:,20]
    
    DISP = np.zeros(( array_gps.shape[0],5 ))
    
    print("-- Triangular dislocation elements will be used")

    X=np.zeros(3)
    Y=np.zeros(3)
    Z=np.zeros(3)
    
    XYZ=np.zeros((3,3))

    for i in range(MM.shape[0]):
    
        if verbose:
            print('  -- ',i,' / ',MM.shape[0])
            
    
        [rdis_long,rdis_lat,rdis_depth,rdis_length,rdis_width,rdis_area,\
         ratio_rdis_tdis,\
         strike,dip,\
         centroid_long,centroid_lat,centroid_depth,\
         tdis_long1,tdis_lat1,tdis_depth1,\
         tdis_long2,tdis_lat2,tdis_depth2,\
         tdis_long3,tdis_lat3,tdis_depth3,\
         tdis_area,rake_1,slip_1,rake_2,slip2,slip] = MM[i]

        (X[0],Y[0],Z[0])=( X_t_dis1[i] , Y_t_dis1[i] , Z_t_dis1[i] )
        (X[1],Y[1],Z[1])=( X_t_dis2[i] , Y_t_dis2[i] , Z_t_dis2[i] )
        (X[2],Y[2],Z[2])=( X_t_dis3[i] , Y_t_dis3[i] , Z_t_dis3[i] )


        X = X * 1.E3
        Y = Y * 1.E3
        Z = Z * 1.E3

        Z=np.sqrt(Z**2)

        print(MM[i,:])
        print(MM[i,23], MM[i,22])
        
        ( ds_mr , ss_mr ) = pyacs.lib.faultslip.slip_rake_2_ds_ss(MM[i,23], MM[i,22])
        ( ds_cr , ss_cr ) = pyacs.lib.faultslip.slip_rake_2_ds_ss(MM[i,25], MM[i,24])
    
        ds = ds_mr + ds_cr
        ss = ss_mr + ss_cr

        print(i,ds_mr,ds)
    
        # change sign for Meade's conventions
        
        ds = -ds
        ss = -ss
        
        U = tde.calc_tri_displacements(SX,SY,SZ, X, Y, Z, Poisson_Ratio, ss, 0.0,  ds)


        DISP[:,2] = DISP[:,2] + U['x'] 
        DISP[:,3] = DISP[:,3] + U['y'] 
        DISP[:,4] = DISP[:,4] + U['z'] 
        
    DISP[:,0] = OBS[:,0]
    DISP[:,1] = OBS[:,1]

    # Meade TDE has opposite sense convention for up displacements
    DISP[:,4] = -DISP[:,4] 

    pyacs.lib.utils.save_np_array_with_string( DISP, NAMES , "%10.6lf %10.6lf %10.3lf %10.3lf %10.3lf   %s", args.output, comment=(" predictions from %s -- triangular dislocations " % args.model))
#    np.savetxt(args.output, DISP, fmt="%10.6lf %10.6lf %10.3lf %10.3lf %10.3lf ", header=(" predictions from %s -- triangular dislocations " % args.model) )


###############################################################################
# RECTANGULAR DISLOCATIONS WITH EDCMP
###############################################################################

if not TDE:

    # array fault must have the following order
    # slip,xf,yf,depth,length,width,strike,dip,rake
    # this corresponds to the following fields from pyacs/pyeq results
    #0:rdis_long     #1:rdis_lat     #2:rdis_depth     #3:rdis_length  
    #4:rdis_width    #5:rdis_area    #6:ratio_rdis_tdis 
    #7:strike        #8:dip        
    #9:centroid_long #10centroid_lat #11centroid_depth     
    #12:tdis_long1   #13:tdis_lat1   #14:tdis_depth1  
    #15:tdis_long2   #16:tdis_lat2   #17:tdis_depth2     
    #18:tdis_long3   #19:tdis_lat3   #20:tdis_depth3
    #21:tdis_area    #22:rake_1      #23:slip_1 
    #24:rake_2       #25:slip2       #26:slip

    
    array_fault_main_rake = MM[:,[23,0,1,2,3,4,7,8,22]]
    # depth must always be positive
    array_fault_main_rake[:,3] = np.sqrt( array_fault_main_rake[:,3]**2 )
    #slip,xf,yf,depth,length,width,strike,dip,rake
    (x,y)=coordinates.geo2flat_earth(MM[:,0],MM[:,1])
    
    array_fault_main_rake[:,1] = x
    array_fault_main_rake[:,2] = y

    
    array_fault_conj_rake = np.copy( array_fault_main_rake )
    array_fault_conj_rake[ :,0] = MM[ :,25]
    array_fault_conj_rake[ :,8] = MM[ :,24]

    import pyeq.lib.edcmp
    print('-- Running edcmp')
    disp_conj_rake = pyeq.lib.edcmp.disp_tilt_strain_stress_from_edcmp(array_fault_conj_rake, array_gps)
    disp_main_rake = pyeq.lib.edcmp.disp_tilt_strain_stress_from_edcmp(array_fault_main_rake, array_gps)

#    DISP = disp_main_rake + disp_conj_rake
    DISP = disp_main_rake + disp_conj_rake
    
    # fills DISP again with the geographical coordinates 
    DISP[:,0] = OBS[:,0]
    DISP[:,1] = OBS[:,1]
    
    pyacs.lib.utils.save_np_array_with_string( DISP, NAMES , "%10.6lf %10.6lf %10.3lf %10.3lf %10.3lf   %s", args.output, comment=(" predictions from %s -- rectangular dislocations " % args.model))
    
#    np.savetxt(args.output, DISP, fmt="%10.6lf %10.6lf %10.3lf %10.3lf %10.3lf ", header=(" predictions from %s -- rectangular dislocations " % args.model) )

print("-- Created %s" % args.output)
