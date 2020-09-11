"""
Print kinematic inversion results
"""

import numpy as np

###############################################################################
def interp_slip_at_dates(SLIP, slip_dates, dates):
###############################################################################
    """
    Given a slip increment provided a a 1D vector, with dates slip_dates, calculates the slip predicted at dates using linear interpolation 
    
    :param SLIP: a 1D numpy vector of slip increments with nstep x nfaults size
    :param slip_dates: a 1D numpy vector of slip dates of size nstep + 1
    :param dates: a 1D numpy vector with interpolation dates
    
    :return: a 2D numpy vector of n_dates as row and nfault as columns
    """
    
    #print slip_dates
    #print dates
    
    # import
    
    import scipy.interpolate
    
    # nstep and nfaults
    
    nstep = slip_dates.shape[0] - 1
    nfaults = int( SLIP.shape[0] / nstep )
    # incremental slip as a 2D array
    
    ISLIP = SLIP.reshape(-1,nfaults).T
    
    # calculates the cumulative slip at slip dates
    
    CSLIP = np.zeros((nfaults,slip_dates.shape[0]))
    for  i in np.arange(slip_dates.shape[0])[1:]:
        CSLIP[:,i] = CSLIP[:,i-1] + ISLIP[:,i-1]
    
    #print 'CSLIP max ', np.max(CSLIP)
        
    # now do the interpolation
    
    f=scipy.interpolate.interp1d(slip_dates,CSLIP,axis=1)
    
    return( f(dates) )

###############################################################################
def predict_disp_from_model(SLIP, G):
###############################################################################
    """
    Given a slip model provided as a 1D vector, with Greens function G, returns the displacement prediction 
    
    np.dot(G,SLIP)
    So the result is a vector with size G.shape[0]
    if SLIP is a matrix with SLIP at different times ordered by columns
    returns np.matmul(G,SLIP)
    
    :param SLIP: a 1D or 2D numpy array of slip
    :param G: Grren's function matrix
    
    :return: a numpy array
    """



###############################################################################
def predict_time_series_from_model(SLIP, slip_dates, GREEN, NAME_OBS, gts):
###############################################################################
    """
    Given 
    - a slip increment provided as a 1D vector, with dates slip_dates, 
    - The Green's function tensor
    - the NAME_OBS vector (giving the index of a site in the Green's tensor)
    - a time series in pyacs gts format
    Returns:
    - a gts time series of model predictions at the dates of the time series
    - a gts time series of model-obs
    """
    
    # import
    
    import numpy as np
    
    # Get the model at the dates of the time series
    
    MODEL = interp_slip_at_dates(SLIP, slip_dates, gts.data[:,0])
    
    # Get the Green matrix for the gts site
    
    site_number=np.argwhere(NAME_OBS==gts.code.upper())[0,0] 
    
    # the 3 rows of GREEN at the main rake (last indice at 0) corresponding to the site 

#    print 'shape of GREEN ', GREEN.shape
    G = GREEN[site_number,:,:,0]
#    print 'shape of G ', G.shape

    # since gts is NEU and GREEN is ENU, swap 
#    print G.T[:,:10]
    NEU_G = G.T.copy()
    NEU_G[0,:] , NEU_G[1,:] = G.T[1,:] , G.T[0,:]  
#    print 'NEU_G ', NEU_G[:,:10]
    # return the product
    
#    print G.shape, MODEL.shape, np.dot(G.T,MODEL).shape
    
    return ( np.dot(NEU_G,MODEL).T )





###############################################################################
def print_kinematic_inversion_results(SLIP, input_npz, inv_CM_step, w, odir, name, np_dates, H_inversion_info,daily):
###############################################################################

    """
    Print kinematic inversion results 4
    """


    print('-- Writing solution in directory',odir)

    # IMPORT
    
    import pyacs.lib.astrotime as AT
    from datetime import datetime,timedelta
    
    from os import mkdir, path
#    import shutil
#    import pyacs.gts
    import pyacs
    from shutil import copyfile
    import copy
    import os
    import shapefile
    from pyeq.lib import eq_disloc_3d as DL
    from pyacs.lib.vel_field import Velocity_Field as VF
    from pyacs.lib.gmtpoint import GMT_Point

    # PARAMETERS

    date_tol=0.00001
    mu=3.0E10

    # DATE FORMAT
    
    s_strdate = AT.decyear2datetime(np_dates[0])
    s_doy, _  = AT.decyear2dayno(np_dates[0])

    e_strdate = AT.decyear2datetime(np_dates[-1])
    e_doy, _  = AT.decyear2dayno(np_dates[-1])
    
    time_step_diff=np.diff(np_dates) * 365.25
    
    median_time_step = np.median(time_step_diff)
    shortest_time_step = np.min(time_step_diff)
    longest_time_step = np.max(time_step_diff)


    # CREATE OUTPUT DIRECTORIES

    odir = odir + '_' + ("%d%03d" % ( int(np_dates[0]) , s_doy )) + '_' + ("%03d" % int(np.rint(shortest_time_step)) ) + 'd_' + ("%d%03d" % ( int(np_dates[-1]) , e_doy ))

    # results directory
    if not path.exists(odir):mkdir(odir)
    
    # info directory
    
    if not path.exists(odir+'/info'):mkdir(odir+'/info')
    
    # GPS time series results directory
    if not path.exists(odir+'/time_series'):mkdir(odir+'/time_series')
    if not path.exists(odir+'/time_series/obs'):mkdir(odir+'/time_series/obs')
    if not path.exists(odir+'/time_series/model'):mkdir(odir+'/time_series/model')
    if not path.exists(odir+'/time_series/res'):mkdir(odir+'/time_series/res')
    if not path.exists(odir+'/time_series/plot'):mkdir(odir+'/time_series/plot')

    # slip results directory
    if not path.exists(odir+'/slip'):mkdir(odir+'/slip')
    if not path.exists(odir+'/slip/incremental'):mkdir(odir+'/slip/incremental')
    if not path.exists(odir+'/slip/cumulative'):mkdir(odir+'/slip/cumulative')
    if not path.exists(odir+'/slip/rate'):mkdir(odir+'/slip/rate')
    
    # slip time series results directory
    if not path.exists(odir+'/slip_time_series'):mkdir(odir+'/slip_time_series')

    # stf results directory
    if not path.exists(odir+'/stf'):mkdir(odir+'/stf')

    # sum results directory
    if not path.exists(odir+'/summary'):mkdir(odir+'/summary')
    fsumn=odir+'/summary/'+name+'_sum.dat'

    # displacements directory
    if not path.exists(odir+'/displacement'):mkdir(odir+'/displacement')
    if not path.exists(odir+'/displacement/cumulative'):mkdir(odir+'/displacement/cumulative')

    # npz directory
    if not path.exists(odir+'/npz'):mkdir(odir+'/npz')




    ###################################################################
    # READING INPUT NPZ FILES
    ###################################################################
    
    import pyeq.lib.npz
    
    GEOMETRY, Dm, GREEN, GREEN_UP, OBS, NAME_OBS, OBS_UP, NAME_OBS_UP = pyeq.lib.npz.read_pyeq_input_npz( input_npz )

    # converts array to recarray
    
    names=['rdis_long','rdis_lat','rdis_depth','rdis_length','rdis_width',\
           'rdis_area','ratio_rdis_tdis','strike','dip',\
           'centroid_long','centroid_lat','centroid_depth',\
           'tdis_long1','tdis_lat1','tdis_depth1',\
           'tdis_long2','tdis_lat2','tdis_depth2',\
           'tdis_long3','tdis_lat3','tdis_depth3','tdis_area']
    
    import pyeq.lib.lib_inversion
    SGEOMETRY=pyeq.lib.lib_inversion.numpy_array_2_numpy_recarray(GEOMETRY,names)

#         try:
#             f = file(input_npy,"rb")
#             GEOMETRY=np.load(f)
#             Dm=np.load(f)
#             GREEN=np.load(f)
#             GREEN_UP=np.load(f)
#             OBS=np.load(f)
#             NAME_OBS=np.load(f)
#             OBS_UP=np.load(f)
#             NAME_OBS_UP=np.load(f)
#             f.close()
#         except:
#             print("!!! Could not read ", input_npy)
#             import sys
#             sys.exit()
#         
#         SGEOMETRY=np.rec.array(GEOMETRY)


    import pyeq.lib.lib_inversion
    GEOMETRY=pyeq.lib.lib_inversion.numpy_recarray_2_numpy_array(SGEOMETRY)
    nfaults=GEOMETRY.shape[0]
    n_fields_geometry=GEOMETRY.shape[1]

    # END INITIALIZATION
    
    ###########################################################################

    # SAVE WARNING FILE
    
    fwarning=open(odir+'/info/warning.dat','w')
    fwarning.write(H_inversion_info['warning'])
    fwarning.close()

    # SAVE SOME NPY FILES in NPY DIR
    
    np.save(odir+'/npz/sol_slip.npz',SLIP)
    np.save(odir+'/npz/np_dates.npz',np_dates)
    copyfile(input_npz, odir+'/npz/input.npz')
#    np.save(odir+'/npy/input.npy', input_npy)

    # SAVE A FILE WITH ALL USEFUL DATES IN INFO DIR

    np_dates_step=copy.deepcopy(np_dates)
    np_dates_step=np_dates_step[1:]
    np_step_duration_day = np.diff(np_dates) * 365.25

    np_dates_text      = np.zeros( (np_dates_step.shape[0] , 8 ) )
    np_dates_text[:,0] = np.arange( np_dates_step.shape[0] )
    np_dates_text[:,1] = np_dates[:-1]
    np_dates_text[:,2] = np_dates[1:]
    np_dates_text[:,3] = np_step_duration_day
    np_dates_text[:,4] = np.trunc( np_dates[:-1] )
    np_dates_text[:,5] = AT.decyear2dayno(np_dates[:-1])[0]
    np_dates_text[:,6] = np.trunc( np_dates[1:] )
    np_dates_text[:,7] = AT.decyear2dayno(np_dates[1:])[0]
    
    str_header='step st_decyear  ed_decyear delta(d) st_date_yr_doy end_date_yr_doy'
    np.savetxt(odir+'/info/dates.dat', np_dates_text, header=str_header, fmt="%05d %10.6lf %10.6lf %5.1lf      %4.0lf %3.0lf        %4.0lf %3.0lf")
    
    # SAVE GEOMETRY IN INFO DIR

    Hdtype={'names':['    rdis_long','   rdis_lat','rdis_depth',\
                    'rdis_length','rdis_width',' rdis_area','ratio_rdis_tdis',\
                    '   strike','       dip',\
                    'centroid_long','centroid_lat','centroid_depth',\
                    'tdis_long1', ' tdis_lat1','tdis_depth1', \
                    'tdis_long2',' tdis_lat2','tdis_depth2', \
                    'tdis_long3',' tdis_lat3','tdis_depth3', \
                    ' tdis_area'],\
             'formats':['f','f','f','f','f','f','f','f','f','f','f','f','f','f','f','f','f','f','f','f','f','f']}
    
    GEOMETRY_TXT=np.zeros((GEOMETRY.shape[0],GEOMETRY.shape[1]+1))
    GEOMETRY_TXT[:,0]=np.arange(GEOMETRY.shape[0])
    GEOMETRY_TXT[:,1:]=GEOMETRY
    header_cols=' '.join(Hdtype['names'])
    format_header='%04d %10.5lf %10.5lf      %6.2lf \
         %6.2lf     %6.2lf    %6.2lf          %6.2lf\
        %6.2lf %10.2lf \
       %10.5lf   %10.5lf         %6.2lf \
    %10.5lf %10.5lf      %6.2lf \
    %10.5lf %10.5lf      %6.2lf \
    %10.5lf %10.5lf      %6.2lf \
       %6.2lf '
    
    np.savetxt(odir+'/info/geometry.dat', GEOMETRY_TXT, fmt=format_header, header=header_cols)
    
    
    # DELTA_SLIP & SLIP
    
    print('-- calculating slip')

    CUMULATIVE_SLIP_TIME_STEP=np.zeros((nfaults,3))
    DELTA_SLIP_TIME_STEP=np.zeros((nfaults,3))

    DELTA_SLIP_TIME_STEP[:,0]=SGEOMETRY.centroid_long
    DELTA_SLIP_TIME_STEP[:,1]=SGEOMETRY.centroid_lat

    CUMULATIVE_SLIP_TIME_STEP[:,0]=SGEOMETRY.centroid_long
    CUMULATIVE_SLIP_TIME_STEP[:,1]=SGEOMETRY.centroid_lat

    DELTA_SLIP_PER_TIME_STEP=SLIP.reshape(-1,nfaults).T
    
    # STF & CSTF

    STF_INC=np.zeros((np_dates_step.shape[0]+1,2))
    STF=np.zeros((np_dates_step.shape[0]+1,2))
    CSTF=np.zeros((np_dates_step.shape[0]+1,2))
    
    STF[0,0]=np_dates[0]
    CSTF[0,0]=np_dates[0]
    
    # SLIP TIME SERIES

    rate_slip_max = 0.
    rate_slip_min = 1.E6
    i_rate_slip_max = 0

    inc_slip_max = 0.
    inc_slip_min = 1.E6
    i_inc_slip_max = 0


    SLIP_TIME_SERIES=np.zeros(( int(SLIP.shape[0]/(np_dates.shape[0]-1)),np_dates.shape[0]))

    date_ref=np_dates[0]

    # loop on dates    

    i=0
    for i in np.arange(np_dates_step.shape[0]):
        date=np_dates_step[i]

        (mday,month,ut)=AT.decyear2cal(date)
        (noday,ut)=AT.decyear2dayno(date)    
        date_info=("time step #%04d date %10.5lf %02d-%02d-%04d-%.1lf doy %03d | days since reference date +%4.1lf " % \
                   (i+1,date,mday,month,int(date),ut,noday,(date-date_ref)*365.25))
        
        # INCREMENTAL SLIP
        DELTA_SLIP_TIME_STEP[:,2]=DELTA_SLIP_PER_TIME_STEP[:,i]
        fname=odir+"/slip/incremental/"+("%04d_delta_slip.dat" % (i+1))
        np.savetxt(fname, DELTA_SLIP_TIME_STEP, fmt="%10.5lf %10.5lf %10.1lf", header=date_info)

        if np.max(DELTA_SLIP_PER_TIME_STEP[:,i]) > inc_slip_max:
            inc_slip_max = np.max(DELTA_SLIP_PER_TIME_STEP[:,i])
            i_inc_slip_max = i
    
        # SLIP RATE
        
        fname=odir+"/slip/rate/"+("%04d_rate_slip.dat" % (i+1))
        np.savetxt(fname, DELTA_SLIP_TIME_STEP / np_step_duration_day[i], fmt="%10.5lf %10.5lf %10.1lf", header=date_info)

        if np.max(DELTA_SLIP_PER_TIME_STEP[:,i] / np_step_duration_day[i] ) > rate_slip_max:
            rate_slip_max = np.max(DELTA_SLIP_PER_TIME_STEP[:,i]  / np_step_duration_day[i] )
            i_rate_slip_max = i
        
        # CUMULATIVE SLIP
        CUMULATIVE_SLIP_TIME_STEP[:,2]=CUMULATIVE_SLIP_TIME_STEP[:,2]+DELTA_SLIP_PER_TIME_STEP[:,i]

        fname=odir+"/slip/cumulative/"+("%04d_cumulative_slip.dat" % (i+1))
        np.savetxt(fname, CUMULATIVE_SLIP_TIME_STEP, fmt="%10.5lf %10.5lf %10.1lf", header=date_info)
        
        # SLIP TIME SERIES
        SLIP_TIME_SERIES[:,i+1]=CUMULATIVE_SLIP_TIME_STEP[:,2]
        
        # STF & CSTF
        
        STF[i+1,0]=date
        STF_INC[i+1,0]=date
        STF_INC[i+1,1]=np.sum(mu*SGEOMETRY.tdis_area*1.0E6*DELTA_SLIP_PER_TIME_STEP[:,i]*1.0E-3)
        STF[i+1,1]=np.sum(mu*SGEOMETRY.tdis_area*1.0E6*DELTA_SLIP_PER_TIME_STEP[:,i]/ np_step_duration_day[i]*1.0E-3)

        CSTF[i+1,0]=date
        CSTF[i+1,1]=np.sum(mu*SGEOMETRY.tdis_area*1.0E6*CUMULATIVE_SLIP_TIME_STEP[:,2]*1.0E-3)

    # SAVE STF & CSTF

    print('-- saving stf & cstf')

    fname=odir+"/stf/STF.dat"
    np.savetxt(fname,STF,fmt="%10.5lf %10.2E")    
        
    fname=odir+"/stf/CSTF.dat"
    np.savetxt(fname,CSTF,fmt="%10.5lf %10.2E")    

    # SAVE SLIP TIME SERIES

    print('-- saving slip time series')

    for isubfault in np.arange(SLIP_TIME_SERIES.shape[0]):
        STS=np.zeros((np_dates.shape[0],2))
        STS[:,0]=np_dates
        STS[:,1]=SLIP_TIME_SERIES[isubfault,:]
        np.savetxt(("%s/slip_time_series/%04d_ts_cumulative_slip.dat" % (odir,isubfault)), STS, fmt="%10.5lf %10.1lf")

    # REGULARIZATION TIME SERIES AND CHI2
    
    import pyeq.lib.lib_inversion
    
    chi2_regularization = 0.0
    
    ts_regularization = np.zeros( (np_dates_step.shape[0] , 12) )
    ts_regularization[:,:np_dates_text.shape[1]] = np_dates_text
    
    import pyacs.lib.glinalg
    CM = pyacs.lib.glinalg.cov_to_invcov(inv_CM_step)
    
    for i in np.arange( np_dates_step.shape[0] ):
        # calculates (m-m0).T Cm^-1*w (m-m0) for step i
        chi2_model = pyeq.lib.lib_inversion.sp(DELTA_SLIP_PER_TIME_STEP[:,i] , w[i] * inv_CM_step)
        ts_regularization[i,8] = chi2_model
        ts_regularization[i,9] = np.mean(np.sqrt(np.diag(CM)))
        ts_regularization[i,10] = np.max(DELTA_SLIP_PER_TIME_STEP[:,i])
        ts_regularization[i,11] = np.sqrt( chi2_model / float( CM.shape[0]) )
        chi2_regularization = chi2_regularization + ts_regularization[i,1]
    
    np.savetxt( odir+'/info/regularization_time_series.dat', ts_regularization , header=str_header+'      m0t Cm-1 m0  sigma  max_slip_inc  reduced_chi2', \
                fmt="%05d %10.6lf %10.6lf %5.1lf      %4.0lf %3.0lf        %4.0lf %3.0lf       %13.3lf %8.3lf     %8.3lf      %8.3lf")
    

    # MODEL AND O-C TIME SERIES
    
    import pyeq.lib.green_tensor
    
    NEW_GREEN = pyeq.lib.green_tensor.GREEN_ENU_TO_GREEN_RAKE_MAIN_RAKE_CONJUGATE(GREEN, GEOMETRY, SGEOMETRY, \
                                    H_inversion_info['rake_type'], H_inversion_info['rake_value'])

    
    
    
    from pyacs.gts.Sgts import Sgts
    from pyacs.gts.Gts import Gts
    import pyacs.lib.coordinates

    tts=Sgts(ts_dir=H_inversion_info['dir_ts'],xyz=False)

    if daily:
        ts=tts.gts('force_daily')
    else:
        ts=tts

    # populates X0, Y0, Z0 of the time series
    
    for gts in ts.lGts():
        site_number=np.argwhere(NAME_OBS==gts.code.upper())[0,0] 
        lon = OBS[site_number,0]
        lat = OBS[site_number,1]
        he  = 0.0
        X0,Y0,Z0 = pyacs.lib.coordinates.geo2xyz(lon, lat, he, unit='dec_deg')
        gts.X0 = X0 
        gts.Y0 = Y0 
        gts.Z0 = Z0 

        gts.lon = lon
        gts.lat = lat
        gts.h  = he
    
    H_code_wrms = {}

    ngps = 0
    ngps_nodata = 0
    chi2_obs = 0
    n_obs = 0
    
    bias_max =np.array([0.,0.,0.])
    bias_max_code = ''
    bias_mean =np.array([0.,0.,0.])

    wrms_max =np.array([0.,0.,0.])
    wrms_max_code = ''

    # loop on time series
    
    # np_dates in seconds
    
    lldate = []
    # get reference date in datetime format
    lref=[]
    for gts in ts.lGts():
        lref.append( gts.data[0,0] )
    
    ref_datetime = AT.decyear2datetime( np.min( np.array(lref) )  )
    
    # get an np_array of dates in seconds    
    np_dates_s = np.array( [] )
    for gts in ts.lGts():
        
        np_dates_datetime = AT.decyear2datetime(gts.data[:,0]) - ref_datetime
        np_dates_s_tmp =  np.array(list(map(int,[x.total_seconds() for x in np_dates_datetime])))
        np_dates_s = np.unique( np.hstack( (np_dates_s , np_dates_s_tmp ) ) )

        
    # list of observed displacement field
    l_obs_disp = {}
    # list of modeled displacement field
    l_model_disp = {}
    
    
    for gts in ts.lGts():
        
        print('-- printing modeled and residual time series for ',gts.code)

        new_gts = gts.extract_periods( [np_dates[0]-date_tol/365.25 , np_dates[-1] + date_tol/365.25] ,verbose = False )
        
        # time series has no data, skip
        if new_gts.data is None:
            ngps_nodata = ngps_nodata + 1
            print('! WARNING: no data found at requested time for site ',gts.code,gts.ifile)
            print('! WARNING: site ',gts.code,' will be skipped in the inversion')
            continue
        # time series has only 1 data, skip
        elif new_gts.data.shape[0] == 1:
            ngps_nodata = ngps_nodata + 1
            print('! WARNING: only one date found at requested time for site ',gts.code,gts.ifile)
            print('! WARNING: site ',gts.code,' will be skipped in the inversion')
            continue
        else:
            # set dN, dE, dU to 0.0 at the first date
            new_gts.data[:,1:4] = new_gts.data[:,1:4] - new_gts.data[0,1:4]
            ngps = ngps + 1
        
        # case up component rescaled
        
        if H_inversion_info['up'][0]:
            new_gts.data[:,6] = new_gts.data[:,6] * 0.0 + H_inversion_info['up'][1]
        
        new_gts.write_pos(idir=odir+'/time_series/obs', force='data' , add_key='obs',)
      
        MODEL_TS = predict_time_series_from_model(SLIP, np_dates, NEW_GREEN, NAME_OBS, new_gts) * 1E-3
        
        # MODEL_TS gives the prediction from the slip since the first date of the model
        # BUT, to compare to observation, we want the model since the first observation time
        
        MODEL_TS = MODEL_TS - MODEL_TS[0,:]
        
        RES_TS   = new_gts.data[:,1:4] - MODEL_TS
        RES_TS = RES_TS * 1.E3
    
        bias_tmp = np.mean(RES_TS,axis=0)
        bias_mean = bias_mean + bias_tmp

        
        if np.sum( bias_tmp**2) > np.sum( bias_max**2 ):
            bias_max = bias_tmp
            bias_max_code = gts.code
            
        
        # chi2_obs
        # if up used        
        if H_inversion_info['up'][0]:
            RES_TS[:,2] = RES_TS[:,2] / H_inversion_info['up'][1]
            chi2_obs = chi2_obs + np.sum(RES_TS**2)
            n_obs = n_obs + RES_TS.shape[0] * 3
        # else up not used
        else:
            RES_TS[:,2] = RES_TS[:,2] * 0.0
            chi2_obs = chi2_obs + np.sum(RES_TS**2)
            n_obs = n_obs + RES_TS.shape[0] * 2
            
        
        p_gts = new_gts.copy()
        p_gts.data[:,1:4] = MODEL_TS
        p_gts.write_pos(idir=odir+'/time_series/model',add_key='model',)
        
        # make png plot of obs time series with model superimposed
        png = odir+'/time_series/plot/'+gts.code+'.png'
        #new_gts.plot(superimposed=p_gts,save=png,center=False)

        # OBS AND MODELED DISPLACEMENTS PSVELO FILES
        
        print('-- making obs and modeled displacements for ' , new_gts.code)
        
        # get seconds
        
        np_dates_datetime = AT.decyear2datetime(new_gts.data[:,0]) - ref_datetime
        np_dates_seconds_gts = np.array(list(map(int,[x.total_seconds() for x in np_dates_datetime])))

        # populates Velocity Fields
        
        for i in np.arange( np_dates_seconds_gts.shape[0] ):
            
            if np_dates_seconds_gts[i] not in l_obs_disp:
                l_obs_disp[ np_dates_seconds_gts[i] ] = VF()
                 
            M = GMT_Point( code=gts.code, lon=new_gts.lon, lat=new_gts.lat, Ve=new_gts.data[i,2]*1.E3 , Vn=new_gts.data[i,1]*1.E3, Vu=new_gts.data[i,3]*1.E3 , SVe=new_gts.data[i,5]*1.E3 , SVn=new_gts.data[i,4]*1.E3, SVu=new_gts.data[i,6]*1.E3)
            l_obs_disp[ np_dates_seconds_gts[i] ].add_point( M )

            if np_dates_seconds_gts[i] not in l_model_disp:
                l_model_disp[ np_dates_seconds_gts[i] ] = VF()
                 
            M = GMT_Point( code=gts.code, lon=p_gts.lon, lat=p_gts.lat, Ve=p_gts.data[i,2]*1.E3 , Vn=p_gts.data[i,1]*1.E3 , Vu=p_gts.data[i,3]*1.E3 , SVe=p_gts.data[i,5]*1.E3 , SVn=p_gts.data[i,4]*1.E3 , SVu=p_gts.data[i,6]*1.E3)
            l_model_disp[ np_dates_seconds_gts[i] ].add_point( M )

        # residual time series                 

        p_gts.data[:,1:4] = RES_TS * 1.E-3
        p_gts.write_pos(idir=odir+'/time_series/res',add_key='res',)

        # wrms for current gts
        H_code_wrms[gts.code] = p_gts.wrms()*1000.

        # update the worst wrms site
        if np.sum( H_code_wrms[gts.code]**2) > np.sum( wrms_max**2 ):
            wrms_max = H_code_wrms[gts.code]
            wrms_max_code = gts.code
 
        # end loop on gts
    
 
    bias_mean = bias_mean / ngps
    
    # SAVE DISPLACEMENT PSVELO FILES
    i=0
    for date_seconds in sorted( l_obs_disp.keys() ):

        cal_date = datetime.isoformat( ref_datetime + timedelta(0, int( date_seconds ) ) ) 
        minutes = date_seconds / 60.
        hours = date_seconds / 3600.
        days    = date_seconds / (3600. * 24. )
        
        my_comment = ("%s | seconds=%d | minutes=%.1lf | hours=%.1f | days=%.2f" % (cal_date,date_seconds,minutes,hours,days) )
        
        print('-- writing psvelo files for ', i, ' / ', len( list(l_obs_disp.keys()) ))
        l_obs_disp[ date_seconds ].write(out_file= odir+'/displacement/cumulative/' + ("%05d_obs.gmt" % i),comment= my_comment , verbose=False )
        l_model_disp[ date_seconds ].write(out_file= odir+'/displacement/cumulative/' + ("%05d_model.gmt" % i),comment= my_comment , verbose=False  )

        l_obs_disp[ date_seconds ].write(out_file= odir+'/displacement/cumulative/' + ("%05d_obs_up.gmt" % i),comment= my_comment , up=True , verbose=False  )
        l_model_disp[ date_seconds ].write(out_file= odir+'/displacement/cumulative/' + ("%05d_model_up.gmt" % i),comment= my_comment , up=True , verbose=False )



        i = i + 1

    # SAVE PYEQ_KINEMATICS COMMAND LINE IN INFO DIR

    print('-- saving pyeq command line in ',odir+'/info/command_line.dat')
    
    fcmd=open(odir+'/info/command_line.dat','w')
    fcmd.write(H_inversion_info['cmd_line'])
    fcmd.close()


    # SAVE WRMS IF RESIDUALS IN INFO DIR

    print("-- saving wrms in ", odir+'/info/wrms.dat' )
    
    fwrms = open(odir+'/info/wrms.dat','w+')
    
    fwrms.write("code  wrms_n (mm) wrms_e (mm) wrms_u (mm)\n")
    fwrms.write("-----------------------------------------\n")
    
    wrms_n_mean = 0
    wrms_e_mean = 0
    wrms_u_mean = 0
    
    for code in sorted(H_code_wrms.keys()):
        fwrms.write("%s %10.1lf %10.1lf %10.1lf\n" % (code,H_code_wrms[code][0], H_code_wrms[code][1], H_code_wrms[code][2] ))

        wrms_n_mean += H_code_wrms[code][0]
        wrms_e_mean += H_code_wrms[code][1]
        wrms_u_mean += H_code_wrms[code][2]

    wrms_n_mean = wrms_n_mean / len(list(H_code_wrms.keys()))
    wrms_e_mean = wrms_e_mean / len(list(H_code_wrms.keys()))
    wrms_u_mean = wrms_u_mean / len(list(H_code_wrms.keys()))

    fwrms.write("-----------------------------------------\n")
    fwrms.write("mean %10.1lf %10.1lf %10.1lf\n" % (wrms_n_mean,wrms_e_mean,wrms_u_mean))

    fwrms.close()

    
    # WRITE SUMMARY

    M0=CSTF[-1,-1]
    magnitude= 2./3.*(np.log10(M0)-9.1)

    slip_max=np.max(CUMULATIVE_SLIP_TIME_STEP[:,2])
    slip_min=np.min(CUMULATIVE_SLIP_TIME_STEP[:,2])


    fsum=open(fsumn,'w')
    
    fsum.write('#--------------------------------------------------------------------------\n')
    fsum.write('# Inversion parameters: \n')
    fsum.write('#--------------------------------------------------------------------------\n')
    fsum.write("    command line                     : %s \n" % (odir+'/info/command_line.dat'))
    fsum.write("    number of parameters estimated   : %d \n" % SLIP.shape[0])
    fsum.write("    memory used in Gb Linux/Mac OS X : %.1lf / %.1lf \n" %  (H_inversion_info['memory_usage'], H_inversion_info['memory_usage']/1024) )
    fsum.write("    inversion duration (s)           : %.1lf \n" % H_inversion_info['time'])
    fsum.write("    linear system building (s)       : %.1lf \n" % H_inversion_info['time_ls'])
    fsum.write("    nnls algorithm                   : %s \n" % H_inversion_info['nnls'])
    fsum.write("    input npz file                   : %s \n" % H_inversion_info['input_npz'])
    fsum.write("    rake_type (Euler or constant)    : %s \n" % H_inversion_info['rake_type'])
    fsum.write("    rake_value (Euler pole or value) : %s \n" % H_inversion_info['rake_value'])
    fsum.write("    rake_constraint (0=fixed, float) : %s \n" % H_inversion_info['rake_constraint'])
    fsum.write("    max_slip (0=from Euler,else user): %s \n" % H_inversion_info['max_slip'])
    fsum.write("    number of subfaults              : %d \n" % SGEOMETRY.shape[0])
    fsum.write("    number of GPS sites used         : %d \n" % ngps)
    fsum.write("    number of GPS sites without data : %d \n" % ngps_nodata)
    fsum.write("    number of time steps             : %d \n" % (SLIP.shape[0] / SGEOMETRY.shape[0]))
    fsum.write("    model start date                 : %.5lf %s doy %03d\n" % (np_dates[0],s_strdate,s_doy))
    fsum.write("    model end   date                 : %.5lf %s doy %03d\n" % (np_dates[-1],e_strdate,e_doy))
    fsum.write("    median,shortest,longest time step: %.1lf %.1lf %.1lf days\n" % (median_time_step,shortest_time_step,longest_time_step))
    fsum.write("    number of GPS horizontal obs.    : %d \n" % OBS.shape[0])
    if H_inversion_info['up'][0]:
        fsum.write("    vertical obs. used / scaling     : YES / %.1lf \n" % H_inversion_info['up'][1])
    else:
        fsum.write("    vertical obs. used / scaling     : NO  /       \n")
        
    
    fsum.write('#--------------------------------------------------------------------------\n')
    fsum.write('#   Inversion results: \n')
    fsum.write('#--------------------------------------------------------------------------\n')
    fsum.write("    Cumulative Moment (N.m)          : %8.2E \n" % M0)
    fsum.write("    Equivalent Mw                    : %8.1f \n" % magnitude)
    i_mr_max = np.argmax(STF[:,1]) 
    
    fsum.write("    max moment rate (N.m / day )     : %8.2E #%05d %d_%03d-%d_%03d\n" % \
               (STF[i_mr_max,1],i_mr_max,np_dates_text[i_mr_max-1,4],np_dates_text[i_mr_max-1,5],\
                np_dates_text[i_mr_max-1,6],np_dates_text[i_mr_max-1,7]))

    magnitude_mr= 2./3.*(np.log10(STF[i_mr_max,1])-9.1)
    fsum.write("    Equivalent Mw                    : %8.1f \n" % magnitude_mr)
    
    fsum.write("    cumulative slip max  (mm)        : %8.1f\n" % slip_max)
    fsum.write("    cumulative slip min  (mm)        : %8.1f\n" % slip_min)
    fsum.write("    incremental slip max  (mm)       : %8.1f #%05d %d_%03d-%d_%03d\n" % (inc_slip_max,i_inc_slip_max,np_dates_text[i_inc_slip_max,4],np_dates_text[i_inc_slip_max,5],np_dates_text[i_inc_slip_max,6],np_dates_text[i_inc_slip_max,7]))
    fsum.write("    rate slip max  (mm/day)          : %8.1f #%05d %d_%03d-%d_%03d\n" % (rate_slip_max,i_rate_slip_max,np_dates_text[i_rate_slip_max,4],np_dates_text[i_rate_slip_max,5],np_dates_text[i_rate_slip_max,6],np_dates_text[i_rate_slip_max,7]))
    
    fsum.write('#--------------------------------------------------------------------------\n')
    fsum.write('#   Inversion statistics: \n')
    fsum.write('#--------------------------------------------------------------------------\n')
    fsum.write("    wrms mean north (mm)             : %8.1f\n" % (wrms_n_mean))
    fsum.write("    wrms mean east  (mm)             : %8.1f\n" % (wrms_e_mean))
    fsum.write("    wrms mean up    (mm)             : %8.1f\n" % wrms_u_mean)

    fsum.write("    wrms max north  (mm)             : %8.1f %s\n" % (wrms_max[0],wrms_max_code))
    fsum.write("    wrms max  east  (mm)             : %8.1f %s\n" % (wrms_max[1],wrms_max_code))
    fsum.write("    wrms max    up  (mm)             : %8.1f %s\n" % (wrms_max[2],wrms_max_code))
    
    fsum.write("    bias mean north (mm)             : %8.1f\n" % bias_mean[0])
    fsum.write("    bias mean  east (mm)             : %8.1f\n" % bias_mean[1])
    fsum.write("    bias mean    up (mm)             : %8.1f\n" % bias_mean[2])

    fsum.write("    bias max  north (mm)             : %8.1f %s\n" % (bias_max[0],bias_max_code))
    fsum.write("    bias max   east (mm)             : %8.1f %s\n" % (bias_max[1],bias_max_code))
    fsum.write("    bias max     up (mm)             : %8.1f %s\n" % (bias_max[2],bias_max_code))

    fsum.write("    chi2 obs                         : %8.1f\n" % chi2_obs)
    fsum.write("    number of observations           : %8d\n" % n_obs)
    fsum.write("    reduced chi2 obs                 : %8.2f\n" % np.sqrt(chi2_obs / float(n_obs) ))
    fsum.write("    chi2 regularization              : %8.1f\n" % chi2_regularization)
    fsum.write("    number of constraints            : %8d\n" % ( SLIP.shape[0] / SGEOMETRY.shape[0] * nfaults) )
    fsum.write("    reduced chi2 regularization      : %8.2f\n" % np.sqrt(chi2_regularization / float(SLIP.shape[0] / SGEOMETRY.shape[0] * nfaults) ))
    fsum.write("    chi2 all                         : %8.1f\n" % (chi2_obs+chi2_regularization) )
    fsum.write("    reduced chi2 all                 : %8.1f\n" % np.sqrt( (chi2_obs+chi2_regularization) / ( float(n_obs) + float(SLIP.shape[0] / SGEOMETRY.shape[0] * nfaults) ) ) )

    fsum.close()
    
    # PRINT RESULTS ON SCREEN
    
    f = open(fsumn, "r")
    text = f.read()
    print(text)
    f.close()
    
    # GENERATE SHAPEFILES OF CUMULATIVE SLIP FOR MODEL DISPLAY IN QGIS

    one_degree=111.1
    TRIANGLE=True    
    
    import glob
    
    lslip_dat=glob.glob(odir+"/slip/cumulative/*_cumulative_slip.dat")

    for slip_dat in lslip_dat:
    
        SLIP=np.genfromtxt(slip_dat,comments='#')

        lfaults=[]
        lrecord=[]
        for i in np.arange(GEOMETRY.shape[0]):
            ( rdis_long, rdis_lat, rdis_depth, rdis_length, rdis_width,  rdis_area, ratio_rdis_tdis,strike,dip,\
             centroid_long, centroid_lat, centroid_depth,\
             tdis_long1,tdis_lat1,tdis_depth1,tdis_long2,tdis_lat2,tdis_depth2,tdis_long3,tdis_lat3,tdis_depth3,tdis_area)\
             =GEOMETRY[i,:]
            depth=rdis_depth
        
    
            if TRIANGLE:
                lfaults.append([ [tdis_long1,tdis_lat1], [tdis_long2,tdis_lat2], [tdis_long3,tdis_lat3] ])
                lrecord.append([i,centroid_depth,SLIP[i,-1],0])
        
            else:
            
                # creates a dislocation object
                disloc=DL.Dislocation(i,rdis_long,rdis_lat,rdis_depth/one_degree, strike, dip, rdis_length/one_degree, rdis_width/one_degree,rdis_area, 0, 0)
                # get the corners
                (X1,X2,X3,X4)=disloc.corners()
                lfaults.append([ [X1[0],X1[1]], [X2[0],X2[1]], [X3[0],X3[1]], [X4[0],X4[1]], [X1[0],X1[1]] ])
                lrecord.append([i,depth,SLIP[i,-1],0])
    
        print("-- ",len(lfaults)," polygons read")

        ###################################################################
        # WRITES GMT PSXY FILES
        # This file can be then plotted with
        #  psxy lima_simple_sigma_010_dc_050_m0_000_sol_coupling.gmt -V -R-81/-72/-16/-7 -JM14 -L -m -W0.2/0  -Clut_coupling.cpt > test.ps
        # triangles have not been tested yet
        ###################################################################
        
        gmtfile=os.path.abspath(slip_dat).split('.')[0]+'.gmt'

        print(("-- saving gmt file %s " % gmtfile.split('/')[-1]))
        f=open(gmtfile,'w')
                
        for i in np.arange(len(lfaults)):
            [index,depth,slip,rake]=lrecord[i]
            f.write('> -Z%.3lf\n'%slip)
            fault=lfaults[i]
            for xy in fault:
                f.write("%10.3lf %10.3lf\n" %(xy[0],xy[1]))
        f.write('>\n')
        f.close()

        ###################################################################
        # WRITES SHAPEFILES
        ###################################################################

        shp_file=os.path.abspath(slip_dat).split('.')[0]
 
       ###################################################################
        # INITIALIZE SHAPEFILE
        ###################################################################
        
        # Make a polygon shapefile
        w = shapefile.Writer( shp_file ,shapeType=shapefile.POLYGON)
        w.field('ID','I','40')
        w.field('i_subfault','F','40')
        w.field('depth_top_disloc','F','40')
        w.field('slip','F','40')
        w.field('rake','F','40')

        
        ###################################################################
        # LOOP ON FAULTS
        ###################################################################
        
        for i in np.arange(len(lfaults)):
            fault=lfaults[i]
            record=lrecord[i]
            w.poly([fault])
            [index,depth,slip,rake]=lrecord[i]
            w.record(str(i),index,depth,slip,rake)
            i=i+1
        
        ###################################################################
        # SAVE SHAPEFILE
        ###################################################################
        
        print(("-- saving shapefile %s " % (shp_file.split('/')[-1]+'.shp')))
        w.close()
        
        ###################################################################
        # SAVE .PRJ FILE
        ###################################################################
        
        prj = open("%s.prj" % shp_file, "w") 
        epsg = 'GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563]],PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433]]' 
        prj.write(epsg) 
        prj.close()


###############################################################################
### END  print_kinematic_inversion_results4 ###################################
###############################################################################
