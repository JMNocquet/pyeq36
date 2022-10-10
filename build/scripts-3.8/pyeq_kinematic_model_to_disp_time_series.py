#!/usr/local/geodesy/anaconda38/bin/python

###################################################################
# SCRIPT    : pyeq_kinematic_model_to_disp_time_series.py
# AUTHOR    : nocquet@geoazur.unice.Fr
# DATE      : June 2020
# OUTPUT    :
# NOTE      :
#
###################################################################


###################################################################
# PYTHON MODULES IMPORT
###################################################################

import argparse, sys
import numpy as np
from colors import red
import os
import pickle

###################################################################
# PYACS ADDITIONAL MODULES IMPORT
###################################################################

import pyacs.lib.coordinates
import pyacs.lib.astrotime as at
from pyacs.lib import gmtpoint as GMT_Point
from pyeq.lib import eq_disloc_3d as Dislocation
import pyeq.elastic_tensor
from pyacs.gts.Sgts import Sgts, Gts
import pyeq.lib.green_tensor

###################################################################
# PARSE ARGUMENT LINE
###################################################################

prog_info="pyeq_kinematic_model_to_disp_time_serie.py :\n"
prog_info+="                     compute the time series predicted by a kinematic model\n"

prog_epilog="J.-M. Nocquet (IRD-Geoazur) - June 2020"

parser = argparse.ArgumentParser(description=prog_info,epilog=prog_epilog)
parser.add_argument('-gps_h', action='store', dest='gmth',required=True,help='gmt psvelo files')
parser.add_argument('-pck', action='store', dest='pck',required=True,help='model pck produced by pyeq_kinematic_inversion.py')
parser.add_argument('--tde', action='count',default=0,help='triangular dislocations computed using Meade (2007)')
parser.add_argument('--verbose', '-v', action='count',default=0,help='verbose mode')
parser.add_argument('-odir', action='store',type=str,required=True, dest='odir',help='output directory')

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
    Poisson_Ratio=0.25
    TDE=True
else:
    TDE=False


###################################################################
# LOAD MODEL PCK
###################################################################

try:
    print("-- Loading %s (%.2f Gb) " % (args.pck, os.path.getsize(args.pck) / 1024 / 1024 / 1024))
    with open(args.pck, "rb") as f:
        model =  pickle.load(f)
    f.close()
except:
    print(red("[PYEQ ERROR] Could not load: %s " % (args.pck)))
    sys.exit()

###################################################################
# READ GPS FILE
###################################################################

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

###################################################################
# MAKE GREEN
###################################################################

GREEN = pyeq.elastic_tensor.make_green(model.geometry, model.sgeometry, array_gps, type='tde', verbose=verbose)
print(GREEN.shape)
print("-- Reformatting the GREEN tensor to account for the main and conjugate rake")
GREEN = pyeq.lib.green_tensor.GREEN_ENU_TO_GREEN_RAKE_MAIN_RAKE_CONJUGATE( \
    GREEN, model.geometry, model.sgeometry, model.rake_type, model.rake_value)

print(GREEN.shape)

# re-order GREEN in component, site, faults
GREEN_REORDERED = np.swapaxes( GREEN[:, :, :, 0], 0, 1).T

###################################################################
# CALCULATES PREDICTIONS
###################################################################

# CS is nfaults, np_model_date_s.shape[0]
CS = np.squeeze(model.CUMULATIVE_SLIP_PER_TIME_STEP).T

# create the interpolated cumulative slip tensor at the observation dates ICS
ICS = np.zeros((model.nfaults, model.np_obs_date_s.shape[0]))
for i in np.arange(model.nfaults):
    ICS[i, :] = np.interp(model.np_obs_date_s, model.np_model_date_s, CS[i, :])

TENSOR_MODEL_TS = np.dot(GREEN_REORDERED, ICS).T

###########################################################################
# PRINT MODEL PREDICTED TIME SERIES
###########################################################################

os.makedirs( args.odir+'/.', exist_ok=True)

Mgts = Sgts(read=False)
# save site coordinates for later use for printing displacement files

COOR = np.zeros((model.np_gps_site.shape[0], 2))

# MODEL PREDICTED DISPLACEMENT TIME SERIES
GREEN = model.green
OBS = model.obs

# TENSOR_MODEL_TS IS ORDERED: component, site_index, date
# print results
for i in np.arange(TENSOR_MODEL_TS.shape[1]):

    site = NAME_OBS[i]
    if model.verbose:
        print("-- printing modeled time series for GPS site: %s" % (site))

    TS = TENSOR_MODEL_TS[:, i, :]

    gts = Gts(code=site)
    site_number = np.argwhere(model.name_obs == site)[0, 0]
    lon = OBS[site_number, 0]
    lat = OBS[site_number, 1]
    #COOR[i, :] = np.array([lon, lat])
    he = 0.0
    X0, Y0, Z0 = pyacs.lib.coordinates.geo2xyz(lon, lat, he, unit='dec_deg')
    gts.X0 = X0
    gts.Y0 = Y0
    gts.Z0 = Z0

    gts.lon = lon
    gts.lat = lat
    gts.h = he

    # All observation dates
    gts.data = np.zeros((TS.shape[0], 10))
    gts.data[:, 0] = at.datetime2decyear(model.np_obs_datetime)
    gts.data[:, 1] = TS[:, 1]  * 1.E-3
    gts.data[:, 2] = TS[:, 0]  * 1.E-3
    gts.data[:, 3] = TS[:, 2]  * 1.E-3
    gts.data[:, 4:] = 1.E-3

    gts.write_pos(idir=args.odir + '/.', force='data', add_key='pred',
                  verbose=False)

