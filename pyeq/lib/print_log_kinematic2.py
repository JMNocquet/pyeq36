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
    G = GREEN[site_number,:,:,0]

    # since gts is NEU and GREEN is ENU, swap 
    NEU_G = G.T.copy()
    NEU_G[0,:] , NEU_G[1,:] = G.T[1,:] , G.T[0,:]  

    # return the product

    return ( np.dot(NEU_G,MODEL).T )





###############################################################################
def print_kinematic_inversion_results2(SLIP, 
                                       input_npz, 
                                       inv_CM_step, 
                                       w, 
                                       odir, 
                                       H_inversion_info,
                                       verbose=False):
###############################################################################

    """
    Print kinematic inversion results for the new pyeq >= 0.50.3
    """


    ###################################################################
    # IMPORT
    ###################################################################


import time

import pyacs.lib.utils
from pyacs.gts.Sgts import Sgts
    from pyacs.gts.Gts import Gts
    import pyacs.lib.coordinates
import pyeq.log.make_dir_pyeq_output
    import pyeq.plot.model2shp_gmt


    ###################################################################
    # PARAMETERS
    ###################################################################

    # shear elastic modulus for moment calculation
    mu=3.0E10


    ###################################################################
    # DATES
    # as input we have two 1D numpy array of dates
    # np_model_date_s: dates at which the model cumulative slip is estimated
    # np_obs_date_s  : dates where at least one observation is available
    # both np_date_model_s & np_date_obs_s are in integer seconds since 1980/1/1 at 00:00:00
    # for printing the results, both dates are converted:
    # delta_d: float being the time in decimal days since the first model date
    # isoformat: is a string YYYY/MM/DAY HR:MN:SC
    # slip rate, stf (moment rate), inc_stf (incremental moment every model time step) correspond to the value during one model time step
    # their date is therefore the middle of the model time step.
    # cumulative slip and and cstf are the estimated model values at the model dates
    ###################################################################
    
    # input dates in integer seconds since 1980/1/1 00:00:00

    np_model_date_s = H_inversion_info['np_model_date_s']
    np_obs_date_s = H_inversion_info['np_obs_date_s']

    # model step duration in seconds
    np_model_step_duration_s =np.diff( np_model_date_s)
    
    # middle of the model dates in seconds
    np_mid_model_date_s = ( np_model_date_s[:-1] +  np.diff( np_model_date_s)/2. ).astype( int )

    # middle of the model dates in datetime
    np_mid_model_datetime = at.seconds2datetime( np_mid_model_date_s )
  
    # middle of the model dates in isoformat
    np_mid_model_isoformat = np.array([x.isoformat(' ') for x in np_mid_model_datetime ])

    # middle of the model dates in delta_days since model first date
    np_mid_model_delta_d = ( np_mid_model_date_s - np_model_date_s[0] ) / 86400.
    
  
    # model step duration in decimal days
    np_model_step_duration_days = np_model_step_duration_s / 86400.
    
    # model dates in datetime
    np_model_datetime = at.seconds2datetime( np_model_date_s ) 

    # model dates in isoformat
    np_model_date_isoformat = np.array([x.isoformat(' ') for x in np_model_datetime ])

    # model dates in delta_days since model first date
    np_model_delta_d = ( np_model_date_s - np_model_date_s[0] ) / 86400.

    # obs dates in datetime
    np_obs_datetime = at.seconds2datetime( np_obs_date_s ) 

    # obs dates in isoformat
    np_obs_date_isoformat = np.array([x.isoformat(' ') for x in np_obs_datetime ])

    # obs dates in delta_days since model first date
    np_obs_delta_d = ( np_obs_date_s - np_model_date_s[0] ) / 86400.

    ###################################################################
    # CREATE OUTPUT DIRECTORIES
    ###################################################################

    syear, sdoy, _ut  = at.datetime2dayno( np_model_datetime[0] )
    eyear, edoy, _ut  = at.datetime2dayno( np_model_datetime[-1] )
     
    median_time_step_days   = np.median( np_model_step_duration_days )
    shortest_time_step_days = np.min( np_model_step_duration_days )
    longest_time_step_days  = np.max( np_model_step_duration_days )


    odir = odir + '_' + ("%d%03d" % ( syear , sdoy )) + '_' + ("%03d" % int(np.rint(shortest_time_step_days)) ) + 'd_' + ("%d%03d" % ( eyear , edoy ))
    print("-- Printing inversion results in: %s" % odir )
    
    
    pyeq.log.make_dir_pyeq_output.make_dir_pyeq_output(odir)

    ###################################################################
    # SAVE PYEQ_KINEMATICS COMMAND LINE IN INFO DIR
    ###################################################################

    print('-- saving pyeq command line in ',odir+'/info/command_line.dat')
    
    fcmd=open(odir+'/info/command_line.dat','w')
    fcmd.write(H_inversion_info['cmd_line'])
    fcmd.close()

    ###########################################################################
    # SAVE GEOMETRY IN INFO DIR
    ###########################################################################

    Hdtype={'names':['    rdis_long','   rdis_lat','rdis_depth',\
                    'rdis_length','rdis_width',' rdis_area','ratio_rdis_tdis',\
                    '   strike','       dip',\
                    'centroid_long','centroid_lat','centroid_depth',\
                    'tdis_long1', ' tdis_lat1','tdis_depth1', \
                    'tdis_long2',' tdis_lat2','tdis_depth2', \
                    'tdis_long3',' tdis_lat3','tdis_depth3', \
                    ' tdis_area'],\
             'formats':['f','f','f','f','f','f','f','f','f','f','f','f','f','f','f','f','f','f','f','f','f','f']}
    
    GEOMETRY = H_inversion_info['GEOMETRY']
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

    
    ###########################################################################
    # SAVE WARNING FILE
    ###########################################################################
    
    fwarning=open(odir+'/info/warning.dat','w')
    fwarning.write(H_inversion_info['warning'])
    fwarning.close()


    ###########################################################################
    # SAVE MODEL DATES IN DIR: INFO
    ###########################################################################

    
    wf = open( odir+'/info/model_dates.dat' , 'w' )
    print("-- saving model dates file in %s" % (odir+'/info/model_dates.dat') )
    
    wf.write("#step  start_date            end_date             delta_d delta_t0 start_doy   end_doy   start_decyear     end_decyear\n")
    for i in np.arange( np_model_date_s.shape[0]-1 ):
        
        sdatetime = np_model_datetime[i]
        edatetime = np_model_datetime[i+1]
        iso_sdate = np_model_date_isoformat[i]
        iso_edate = np_model_date_isoformat[i+1]
        delta_d   = ( np_model_date_s[i+1] - np_model_date_s[i] ) / 86400.
        delta_t0  = ( np_model_date_s[i+1] - np_model_date_s[0] ) / 86400.
        syear, sdoy, _ut  = at.datetime2dayno( sdatetime )
        eyear, edoy, _ut  = at.datetime2dayno( edatetime )

        sdecyear = at.datetime2decyear( sdatetime )
        edecyear = at.datetime2decyear( edatetime )

        wf.write("%04d  %s -> %s %8.3lf %8.3lf %04d %03d -> %04d %03d %15.10lf -> %15.10lf \n" % ( i, iso_sdate, iso_edate, delta_d, delta_t0, syear, sdoy, eyear, edoy, sdecyear, edecyear ))

    wf.close()

    
    ###########################################################################
    # SAVE OBSERVATION DATES IN DIR: INFO
    ###########################################################################


    wf = open( odir+'/info/obs_dates.dat' , 'w' )
    print("-- saving observation dates file in %s" % (odir+'/info/obs_dates.dat') )
    
    wf.write("#obs  date                 delta_t0  year doy      decyear\n")
    for i in np.arange( np_obs_date_s.shape[0] ):
        
        sdatetime = np_obs_datetime[i]
        iso_sdate = np_obs_date_isoformat[i]
        delta_t0  = np_obs_delta_d[i]
        syear, sdoy, _ut  = at.datetime2dayno( sdatetime )
        sdecyear = at.datetime2decyear( sdatetime )

        wf.write("%04d  %s  %8.3lf  %04d %03d %15.10lf\n" % ( i, iso_sdate, delta_t0, syear, sdoy, sdecyear ))

    wf.close()
    
    
    ###########################################################################
    # NPY TENSORS
    # RATE_SLIP, DELTA_SLIP & SLIP
    # HANDLES ORIGIN TIME CONSTANT
    # HANDLES VARIABLES RAKE CASE
    # SAVE SOME NPY FILES in NPY DIR
    ###########################################################################
    
    # SOLUTION
    
    np.save(odir+'/npy/SOL.npy',SLIP)
    
    # COPY INPUT NPZ
    copyfile(input_npz, odir+'/npy/input.npz')

    # OBS TENSOR
    ###########################################################################
    np.save(odir+'/npy/T_OBS.npy', H_inversion_info['T_OBS'])
    
    # GREEN TENSOR
    ###########################################################################
    np.save(odir+'/npy/GREEN.npy', H_inversion_info['GREEN'])

    
    print("-- saving tensors as npy in dir: %s" % (odir+'/npy') )

    # initialization
    np_rake = H_inversion_info['np_rake']
    nfaults = H_inversion_info['GREEN'].shape[1]
    n_constant = H_inversion_info['n_constant']
    
    # TENSOR RATE_SLIP_PER_TIME_STEP
    if n_constant != 0:
        RATE_SLIP_PER_TIME_STEP = SLIP[:-n_constant].reshape( -1, nfaults , np_rake.shape[0]  )
    else:
        RATE_SLIP_PER_TIME_STEP =SLIP.reshape( -1, nfaults , np_rake.shape[0] )

    # TENSOR INC_SLIP_PER_TIME_STEP
    INC_SLIP_PER_TIME_STEP = ( RATE_SLIP_PER_TIME_STEP.T * np_model_step_duration_days ).T
    
    # TENSOR CUMULATIVE_SLIP_PER_TIME_STEP
    CUMULATIVE_SLIP_PER_TIME_STEP = np.zeros( (INC_SLIP_PER_TIME_STEP.shape[0]+1, INC_SLIP_PER_TIME_STEP.shape[1], INC_SLIP_PER_TIME_STEP.shape[2] ) )
    CUMULATIVE_SLIP_PER_TIME_STEP[1:,:,:] = np.cumsum( INC_SLIP_PER_TIME_STEP , axis=0 )
    
    print("-- saving slip_rate.npy [index_date,fault,rake] tensor associated slip_rate_dates.npy (datetime format) in %s " % (odir+'/npy') )

    print("-- slip_rate.npy shape " , RATE_SLIP_PER_TIME_STEP.shape )
    np.save( odir+'/npy/slip_rate.npy' , RATE_SLIP_PER_TIME_STEP )
    np.save(  odir+'/npy/slip_rate_datetime.npy' , np_mid_model_datetime )
    np.save(  odir+'/npy/slip_rate_delta_d.npy' , np_mid_model_delta_d )

    print("-- incremental_slip.npy shape " , INC_SLIP_PER_TIME_STEP.shape )
    np.save( odir+'/npy/incremental_slip.npy' , INC_SLIP_PER_TIME_STEP )
    np.save(  odir+'/npy/incremental_slip_datetime.npy' , np_mid_model_datetime )
    np.save(  odir+'/npy/incremental_slip_delta_d.npy' , np_mid_model_delta_d )

    print("-- cumulative_slip.npy shape " , CUMULATIVE_SLIP_PER_TIME_STEP.shape )
    np.save( odir+'/npy/cumulative_slip.npy' , CUMULATIVE_SLIP_PER_TIME_STEP )
    np.save(  odir+'/npy/cumulative_slip_datetime.npy' , np_model_datetime )
    np.save(  odir+'/npy/cumulative_slip_delta_d.npy' ,  np_model_delta_d )

    # Geometry
    GEOMETRY = H_inversion_info['GEOMETRY']
    np.save(  odir+'/npy/geometry.npy' ,  GEOMETRY )


    ###########################################################################
    # SAVE SLIP TIME SERIES
    ###########################################################################

    print('-- saving rate, incremental and cumulative slip time series')

    for isubfault in np.arange(  nfaults ):
        # SLIP RATE
        DELTA_D_AND_STS_RATE = np.vstack( ( np_mid_model_delta_d , RATE_SLIP_PER_TIME_STEP[:,isubfault,:].T ) ).T
        fname_rate = ("%s/slip_time_series/rate/%04d_rate.dat" % (odir,isubfault))
        # 1 or 2 rakes
        if RATE_SLIP_PER_TIME_STEP.shape[2] == 1:
            format = "%5.2lf  %10.3lf  %s"
        else:
            format = "%5.2lf  %10.3lf  %10.3lf  %s"
            
        comment_rate = ("subfault: #%04d . units are mm/day. dates are at the middle of model time steps. Decimal days since model date start: %s" % ( isubfault, np_model_date_isoformat[0] ) )

        pyacs.lib.utils.save_np_array_with_string( DELTA_D_AND_STS_RATE , \
                                                  np_mid_model_isoformat ,\
                                                  format, \
                                                  fname_rate, \
                                                  comment_rate )


        # INCREMENTAL SLIP
        DELTA_D_AND_STS_INC = np.vstack( ( np_mid_model_delta_d , INC_SLIP_PER_TIME_STEP[:,isubfault,:].T ) ).T
        fname_inc  = ("%s/slip_time_series/incremental/%04d_incremental.dat" % (odir,isubfault))
        comment_inc  = ("subfault: #%04d . units are mm. dates are at the middle of model time steps. Decimal days since model date start: %s" % ( isubfault, np_model_date_isoformat[0] ) )
        pyacs.lib.utils.save_np_array_with_string( DELTA_D_AND_STS_INC , \
                                                  np_mid_model_isoformat ,\
                                                  format, \
                                                  fname_inc, \
                                                  comment_inc )

        # CUMULATIVE SLIP
        DELTA_D_AND_STS_CUM = np.vstack( ( np_model_delta_d , CUMULATIVE_SLIP_PER_TIME_STEP[:,isubfault,:].T ) ).T
        fname_cum  = ("%s/slip_time_series/cumulative/%04d_cumulative.dat" % (odir,isubfault))
        comment_cum  = ("subfault: #%04d . units are mm. dates are model dates. Decimal days since model date start: %s" % \
                        ( isubfault, np_model_date_isoformat[0] ) )
        pyacs.lib.utils.save_np_array_with_string( DELTA_D_AND_STS_CUM , \
                                                  np_model_date_isoformat ,\
                                                  format, \
                                                  fname_cum, \
                                                  comment_cum )


    ###########################################################################
    # STF, INC_STF, CSTF as npy
    ###########################################################################
    # area in meter-square
    AREA = H_inversion_info['SGEOMETRY'].tdis_area * 1.0E6
    
    # STF
    # norm of slip if variable rake
    NORM_RATE_SLIP_PER_TIME_STEP = np.sqrt( np.sum( RATE_SLIP_PER_TIME_STEP**2, axis=2) )
    STF = mu * 1.E-3 * np.sum( ( NORM_RATE_SLIP_PER_TIME_STEP * AREA ) , axis= 1)
    np.save( odir+'/npy/stf.npy' , STF )

    # INC_STF    
    NORM_INC_SLIP_PER_TIME_STEP = np.sqrt( np.sum( INC_SLIP_PER_TIME_STEP**2, axis=2) )
    INC_STF = mu * 1.E-3 * np.sum( ( NORM_INC_SLIP_PER_TIME_STEP * AREA ) , axis= 1)
    np.save( odir+'/npy/incremental_stf.npy' , INC_STF )

    # CSTF
    NORM_CUMULATIVE_SLIP_PER_TIME_STEP = np.sqrt( np.sum( CUMULATIVE_SLIP_PER_TIME_STEP**2, axis=2) )
    CSTF = mu * 1.E-3 * np.sum( ( NORM_CUMULATIVE_SLIP_PER_TIME_STEP * AREA ) , axis= 1)
    np.save( odir+'/npy/cstf.npy' , CSTF )

    ###########################################################################
    # STF, INC_STF, CSTF as text files
    ###########################################################################
    
    # STF
    DELTA_D_AND_STF = np.vstack( ( np_mid_model_delta_d , STF) ).T
    format = "%5.2lf  %.4E  %s"
    comment = ("dates are at the middle of model time steps. Decimal days since model date start: %s" % \
               ( np_model_date_isoformat[0] ) )

    pyacs.lib.utils.save_np_array_with_string( DELTA_D_AND_STF , \
                                               np_mid_model_isoformat ,\
                                              format, \
                                              odir+'/stf/stf.dat', \
                                              comment )
    # INC_STF
    DELTA_D_AND_INC_STF = np.vstack( ( np_mid_model_delta_d , INC_STF) ).T
    format = "%5.2lf  %.4E  %s"
    comment = ("dates are the end date of model time steps. Decimal days since model date start: %s" % ( np_model_date_isoformat[0] ) )

    pyacs.lib.utils.save_np_array_with_string( DELTA_D_AND_INC_STF , \
                                               np_mid_model_isoformat ,\
                                              format, \
                                              odir+'/stf/incremental_stf.dat', \
                                              comment )

    # CSTF
    DELTA_D_AND_CSTF = np.vstack( ( np_model_delta_d , CSTF) ).T
    format = "%5.2lf  %.4E  %s"
    comment = ("Decimal days since model date start: %s" % ( at.seconds2datetime( np_model_date_s[0] ).isoformat(' ')))

    pyacs.lib.utils.save_np_array_with_string( DELTA_D_AND_CSTF , \
                                               np_model_date_isoformat ,\
                                              format, \
                                              odir+'/stf/cstf.dat', \
                                              comment )


    
    # SLIP TIME SERIES

#    SLIP_TIME_SERIES=np.zeros(( int(SLIP.shape[0]/(np_dates.shape[0]-1)),np_dates.shape[0]))

#    rate_slip_max = 0.
#    rate_slip_min = 1.E6
#    i_rate_slip_max = 0

#    inc_slip_max = 0.
#    inc_slip_min = 1.E6
#    i_inc_slip_max = 0


#    date_ref=np_dates[0]

    ###########################################################################
    # INITIALIZE COORDINATES FROM GEOMETRY FOR PRINTING SLIP
    ###########################################################################

    RATE_SLIP_TIME_STEP = np.zeros((nfaults, 2 + np_rake.shape[0] ))
    RATE_SLIP_TIME_STEP[:,0] = H_inversion_info['SGEOMETRY'].centroid_long
    RATE_SLIP_TIME_STEP[:,1] = H_inversion_info['SGEOMETRY'].centroid_lat

    INC_SLIP_TIME_STEP = np.zeros((nfaults, 2 + np_rake.shape[0] ))
    INC_SLIP_TIME_STEP[:,0] = H_inversion_info['SGEOMETRY'].centroid_long
    INC_SLIP_TIME_STEP[:,1] = H_inversion_info['SGEOMETRY'].centroid_lat


    CUMULATIVE_SLIP_TIME_STEP=np.zeros((nfaults, 2 + np_rake.shape[0] ))
    CUMULATIVE_SLIP_TIME_STEP[:,0] = H_inversion_info['SGEOMETRY'].centroid_long
    CUMULATIVE_SLIP_TIME_STEP[:,1] = H_inversion_info['SGEOMETRY'].centroid_lat

    
    ###########################################################################
    # LOOP ON MODEL TIME STEP FOR SLIP RATE & INCREMENTAL SLIP
    ###########################################################################
    for i in np.arange( np_model_step_duration_s.shape[0] ):
        print("--- calculating incremental slip, slip rate for model time step: %d" % (i) )
        sdatetime = np_model_datetime[i]
        edatetime = np_model_datetime[i+1]
        iso_sdate = np_model_date_isoformat[i]
        iso_edate = np_model_date_isoformat[i+1]
        delta_d   = np_model_step_duration_days[i]
        delta_t0  = np_model_delta_d[i]
        syear, sdoy, _ut  = at.datetime2dayno( sdatetime )
        eyear, edoy, _ut  = at.datetime2dayno( edatetime )

        sdecyear = at.datetime2decyear( sdatetime )
        edecyear = at.datetime2decyear( edatetime )

        date_info=("step #%04d  %s -> %s %8.3lf %8.3lf %04d %03d -> %04d %03d %15.10lf -> %15.10lf \n" % \
                   ( i, iso_sdate, iso_edate, delta_d, delta_t0, syear, sdoy, eyear, edoy, sdecyear, edecyear ))
        
        # SLIP RATE
        RATE_SLIP_TIME_STEP[ : , 2: ]  =  RATE_SLIP_PER_TIME_STEP[i].reshape( nfaults, -1 )
        fname=odir+"/slip/rate/"+("%04d_slip_rate.dat" % ( i ))
        # 1 rake
        if RATE_SLIP_TIME_STEP.shape[1] == 3:
            np.savetxt(fname, RATE_SLIP_TIME_STEP, fmt="%10.5lf %10.5lf %10.3lf", header=date_info)
        # 2 rake
        else:
            np.savetxt(fname, RATE_SLIP_TIME_STEP, fmt="%10.5lf %10.5lf %10.3lf %10.3lf ", header=date_info)
            
        # INCREMENTAL SLIP
        INC_SLIP_TIME_STEP[:,2:] = INC_SLIP_PER_TIME_STEP[i].reshape( nfaults, -1 )
        fname=odir+"/slip/incremental/"+("%04d_delta_slip.dat" % ( i ))
        # 1 rake
        if INC_SLIP_TIME_STEP.shape[1] == 3:
            np.savetxt(fname, INC_SLIP_TIME_STEP, fmt="%10.5lf %10.5lf %10.3lf", header=date_info)
        # 2 rake
        else:
            np.savetxt(fname, INC_SLIP_TIME_STEP, fmt="%10.5lf %10.5lf %10.3lf %10.3lf ", header=date_info)
        
    # LOOP ON MODEL DATES FOR THE CUMULATIVE SLIP
    for i in np.arange( np_model_date_s.shape[0] ):

        isodate = np_model_date_isoformat[i]
        print("--- writing cumulative slip for model time : %04d %s" % (i ,  isodate ) )

        date_info=("# model date  %04d %s" % (i ,  isodate ) )

        # CUMULATIVE SLIP
        CUMULATIVE_SLIP_TIME_STEP[:,2:] = CUMULATIVE_SLIP_PER_TIME_STEP[i].reshape( nfaults, -1 )
        fname=odir+"/slip/cumulative/"+("%04d_cumulative_slip.dat" % (i))
        # 1 rake
        if CUMULATIVE_SLIP_TIME_STEP.shape[1] == 3:
            np.savetxt(fname, CUMULATIVE_SLIP_TIME_STEP, fmt="%10.5lf %10.5lf %10.3lf", header=date_info)
        # 2 rake
        else:
            np.savetxt(fname, CUMULATIVE_SLIP_TIME_STEP, fmt="%10.5lf %10.5lf %10.3lf %10.3lf ", header=date_info)

    ###########################################################################
    # REGULARIZATION TIME SERIES AND CHI2
    # FROM BUILD 1 - COMMENTED FOR BUILD 2
    ###########################################################################
    
#     chi2_regularization = 0.0
#     
#     ts_regularization = np.zeros( (np_dates_step.shape[0] , 12) )
#     ts_regularization[:,:np_dates_text.shape[1]] = np_dates_text
#     
#     CM = pyacs.lib.glinalg.cov_to_invcov(inv_CM_step)
#     
#     for i in np.arange( np_dates_step.shape[0] ):
#         # calculates (m-m0).T Cm^-1*w (m-m0) for step i
#         chi2_model = pyeq.lib.lib_inversion.sp(DELTA_SLIP_PER_TIME_STEP[:,i] , w[i] * inv_CM_step)
#         ts_regularization[i,8] = chi2_model
#         ts_regularization[i,9] = np.mean(np.sqrt(np.diag(CM)))
#         ts_regularization[i,10] = np.max(DELTA_SLIP_PER_TIME_STEP[:,i])
#         ts_regularization[i,11] = np.sqrt( chi2_model / float( CM.shape[0]) )
#         chi2_regularization = chi2_regularization + ts_regularization[i,1]
#     
#     np.savetxt( odir+'/info/regularization_time_series.dat', ts_regularization , header=str_header+'      m0t Cm-1 m0  sigma  max_slip_inc  reduced_chi2', \
#                 fmt="%05d %10.6lf %10.6lf %5.1lf      %4.0lf %3.0lf        %4.0lf %3.0lf       %13.3lf %8.3lf     %8.3lf      %8.3lf")
    

    ###########################################################################
    # MODEL PREDICTED TIME SERIES
    ###########################################################################
    
    Mgts = Sgts( read=False )
    # save site coordinates for later use for printing displacement files
    
    COOR = np.zeros( (  H_inversion_info['np_gps_site'].shape[0] , 2 ) )
    
    # MODEL PREDICTED DISPLACEMENT TIME SERIES
    GREEN = H_inversion_info['GREEN']
    OBS = H_inversion_info['OBS']
    
    # case main rake only
    if CUMULATIVE_SLIP_PER_TIME_STEP.shape[2] ==1:
        
        # CS is nfaults, np_model_date_s.shape[0]
        CS = np.squeeze( CUMULATIVE_SLIP_PER_TIME_STEP ).T

        # create the interpolated cumulative slip tensor at the observation dates ICS
        ICS = np.zeros( ( nfaults, np_obs_date_s.shape[0] ) )
        for i in np.arange( nfaults ):
            ICS[i,:] = np.interp( np_obs_date_s, np_model_date_s, CS[i,:] )
            
        # re-order GREEN in component, site, faults
        GREEN_REORDERED = np.swapaxes( GREEN[:,:,:,0] , 0 , 1 ).T

        # calculates prediction
        TENSOR_MODEL_TS = np.dot( GREEN_REORDERED , ICS  ).T
        
        # TENSOR_MODEL_TS IS ORDERED: component, site_index, date
        # print results
        for i in np.arange( TENSOR_MODEL_TS.shape[1] ):
            
            site = H_inversion_info['np_gps_site'][i]
            print("-- printing modeled time series for GPS site: %s" % ( site ) )
            
            TS = TENSOR_MODEL_TS[:,i,:]
            #np.savetxt(odir+'/time_series/model/'+site+'.dat' , TS, fmt="%10.5lf %10.5lf  %10.5lf ")
            gts = Gts( code = site )
            site_number=np.argwhere( H_inversion_info['NAME_OBS'] == site )[0,0]
            lon = OBS[site_number,0]
            lat = OBS[site_number,1]
            COOR[ i , : ] = np.array([ lon , lat ])
            he  = 0.0
            X0,Y0,Z0 = pyacs.lib.coordinates.geo2xyz(lon, lat, he, unit='dec_deg')
            gts.X0 = X0 
            gts.Y0 = Y0 
            gts.Z0 = Z0 
    
            gts.lon = lon
            gts.lat = lat
            gts.h  = he
            
            # All observation dates
            gts.data = np.zeros( (TS.shape[0] , 10 ) )
            gts.data[:,0] = at.datetime2decyear( np_obs_datetime ) 
            gts.data[ :,1]  = TS[:,1] *1.E-3
            gts.data[ :,2]  = TS[:,0] *1.E-3
            gts.data[ :,3]  = TS[:,2] *1.E-3
            gts.data[ :,4:] = 1.E-3
            
            gts.write_pos(idir=odir+'/time_series/model_all_dates', force='data' , add_key='model_all_dates' , verbose=False)

            # observation dates available for the current site
            TS_with_NAN = H_inversion_info['T_OBS'][:,i,:]

            # remove lines with Nan
            lindex = np.argwhere(np.isnan(TS_with_NAN[:,0]))
            TS = np.delete( TENSOR_MODEL_TS[:,i,:], lindex, axis=0)        

            # fills and write gts
            gts.data = np.zeros( (TS.shape[0] , 10 ) )
            gts.data[:,0] = np.delete( at.datetime2decyear( np_obs_datetime ) , lindex ) 
            gts.data[ :,1]  = TS[:,1] *1.E-3
            gts.data[ :,2]  = TS[:,0] *1.E-3
            gts.data[ :,3]  = TS[:,2] *1.E-3
            gts.data[ :,4:] = 1.E-3
            
            gts.write_pos(idir=odir+'/time_series/model', force='data' , add_key='model' , verbose=False)
            
            # save for later use
            Mgts.append( gts )
            
    ###########################################################################
    # OBS TIME SERIES REALLY USED IN THE INVERSION
    ###########################################################################

    Ogts = Sgts( read=False )

    
    for i in np.arange( H_inversion_info['T_OBS'].shape[1] ):
        site = H_inversion_info['np_gps_site'][i]
        print("-- writing observation time series for site: %s" % ( site ) )
        TS_with_NAN = H_inversion_info['T_OBS'][:,i,:]
        # remove lines with Nan
        lindex = np.argwhere(np.isnan(TS_with_NAN[:,0]))
        TS = np.delete(TS_with_NAN, lindex, axis=0)        
        
        # fill the gts
        gts = Gts( code = site )
        site_number=np.argwhere( H_inversion_info['NAME_OBS'] == site )[0,0]
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
        
        gts.data = np.zeros( (TS.shape[0] , 10 ) )
        gts.data[:,0] = at.datetime2decyear( np.delete(np_obs_datetime, lindex) ) 
        
        gts.data[ :,1]  = TS[:,1] *1.E-3
        gts.data[ :,2]  = TS[:,0] *1.E-3
        gts.data[ :,3]  = TS[:,2] *1.E-3
        
        gts.data[ :,4] = TS[:,4] *1.E-3
        gts.data[ :,5] = TS[:,3] *1.E-3
        gts.data[ :,6] = TS[:,5] *1.E-3
        
        gts.write_pos(idir=odir+'/time_series/obs', force='data' , add_key='obs' , verbose=False)

        Ogts.append( gts )

    ###########################################################################
    # RESIDUAL TIME SERIES
    ###########################################################################

    stats_site = np.zeros( ( H_inversion_info['np_gps_site'].shape[0] , 9 ) )

    RES = np.zeros( (0, 10 ) )

    for site in H_inversion_info['np_gps_site']:
        
        print("-- writing residual time series for site: %s "  % site )

        # copy obs to res
        res_gts = Ogts.__dict__[site].copy()

        # remove model prediction
        # caution: model prediction needs to be set to 0 for the first available date
        # this probably needs to be changed if an origin constant has been estimated
        res_gts.data[:,1:4] = res_gts.data[:,1:4] - ( Mgts.__dict__[site].data[:,1:4] - Mgts.__dict__[site].data[0,1:4] )

        # write pos file
        res_gts.write_pos(idir=odir+'/time_series/res', force='data' , add_key='res' , verbose=False)
        
        # UPDATE RES
        
        RES = np.vstack( (RES, res_gts.data * 1.E3 ) )
        
        # get some statistics
        # bias per component
        bias = np.mean( res_gts.data[:,1:4] , axis= 0)
        
        # rms per component
        rms = np.std(res_gts.data[:,1:4], axis=0 )
        
        # wrms per component
        wrms = np.sqrt( np.sum( (res_gts.data[:,1:4] / res_gts.data[:,4:7])**2 , axis=0 ) / np.sum( ( 1./ res_gts.data[:,4:7])**2 , axis=0) )
        
        # save stats for later use
        stats_site[i] = np.array([bias, rms, wrms]).flatten()

    ###########################################################################
    # MODEL/OBS/RESIDUAL DISPLACEMENT FIELDS FOR EACH OBSERVATION DATE
    ###########################################################################

    format_psvelo="%10.5lf %10.5lf %10.3lf %10.3lf   %10.3lf %10.3lf %3.1lf %s"
    
    # LOOP ON OBS DATES
    for i in np.arange( np_obs_date_s.shape[0] ):

        DISP_with_NAN = H_inversion_info['T_OBS'][i,:,:]
        lindex = np.argwhere(np.isnan(DISP_with_NAN[:,0]))
        np_gps_site_i = np.delete( H_inversion_info['np_gps_site'] , lindex )
        COOR_i = np.delete( COOR , lindex , axis=0  )

        # MODEL
        
        print("-- printing modeled displacement at obs date #%04d = %s" % ( i , np_obs_date_isoformat[i] ) )
        DISP = TENSOR_MODEL_TS[i,:,:]
        DISP = np.delete( DISP , lindex ,axis=0 )
        
        # Horizontal
        COOR_DISP_MODEL = np.zeros( ( DISP.shape[0] , 7 ) )
        COOR_DISP_MODEL[ : , :2 ] = COOR_i
        COOR_DISP_MODEL[ : , 2 ] = DISP[:, 0]
        COOR_DISP_MODEL[ : , 3 ] = DISP[:, 1]
        COOR_DISP_MODEL[ : , 4: ] = np.array([0. , 0., 0.])

        fname = ("%s/displacement/cumulative/model/%04d_model_cum_disp.dat" % (odir, i))
        comment = (" model cumulative displacement at obs date #%04d %s" % ( i, np_obs_date_isoformat[i] ))
        
        
        pyacs.lib.utils.save_np_array_with_string( COOR_DISP_MODEL , \
                                                   np_gps_site_i ,\
                                                  format_psvelo, \
                                                  fname, \
                                                  comment )

        # Up
        COOR_DISP_MODEL_UP = np.zeros( ( DISP.shape[0] , 7 ) )
        COOR_DISP_MODEL_UP[ : , :2 ] = COOR_i
        COOR_DISP_MODEL_UP[ : , 3 ] = DISP[:, 2]
        COOR_DISP_MODEL_UP[ : , 4: ] = np.array([0. , 0., 0.])

        fname = ("%s/displacement/cumulative/model/%04d_model_cum_disp_up.dat" % (odir, i))
        comment = (" model up cumulative displacement at obs date #%04d %s" % ( i, np_obs_date_isoformat[i] ))
        
        
        pyacs.lib.utils.save_np_array_with_string( COOR_DISP_MODEL_UP , \
                                                   np_gps_site_i ,\
                                                  format_psvelo, \
                                                  fname, \
                                                  comment )


        # OBSERVATION
        
        print("-- printing observation displacement at obs date #%04d = %s" % ( i , np_obs_date_isoformat[i] ) )
        

        # remove sites with Nan
        DISP = np.delete(DISP_with_NAN, lindex, axis=0)      
        
        # Horizontal
        COOR_DISP_OBS = np.zeros( ( DISP.shape[0] , 7 ) )
        COOR_DISP_OBS[ : , :2 ] = COOR_i
        COOR_DISP_OBS[ : , 2 ] = DISP[:, 0]
        COOR_DISP_OBS[ : , 3 ] = DISP[:, 1]
        COOR_DISP_OBS[ : , 4 ] = DISP[:, 3]
        COOR_DISP_OBS[ : , 5 ] = DISP[:, 4]

        fname = ("%s/displacement/cumulative/obs/%04d_obs_cum_disp.dat" % (odir, i))
        comment = (" observed cumulative displacement at obs date #%04d %s" % ( i, np_obs_date_isoformat[i] ))
        
        
        pyacs.lib.utils.save_np_array_with_string( COOR_DISP_OBS , \
                                                   np_gps_site_i ,\
                                                  format_psvelo, \
                                                  fname, \
                                                  comment )

        # Up
        COOR_DISP_OBS_UP = np.zeros( ( DISP.shape[0] , 7 ) )
        COOR_DISP_OBS_UP[ : , :2 ] = COOR_i
        COOR_DISP_OBS_UP[ : , 3 ] = DISP[:, 2]
        COOR_DISP_OBS_UP[ : , 5 ] = DISP[:, 5]

        fname = ("%s/displacement/cumulative/obs/%04d_obs_cum_disp_up.dat" % (odir, i))
        comment = (" observed up cumulative displacement at obs date #%04d %s" % ( i, np_obs_date_isoformat[i] ))
        
        
        pyacs.lib.utils.save_np_array_with_string( COOR_DISP_OBS_UP , \
                                                   np_gps_site_i ,\
                                                  format_psvelo, \
                                                  fname, \
                                                  comment )

        # RESIDUALS

        print("-- printing residual displacement at obs date #%04d = %s" % ( i , np_obs_date_isoformat[i] ) )

        # Horizontal        
        COOR_DISP_RES = COOR_DISP_OBS - COOR_DISP_MODEL
        COOR_DISP_RES[:,:2] = COOR_DISP_OBS[:,:2]
    
        fname = ("%s/displacement/cumulative/res/%04d_res_cum_disp.dat" % (odir, i))
        comment = (" observed residual cumulative displacement at obs date #%04d %s" % ( i, np_obs_date_isoformat[i] ))
        
        pyacs.lib.utils.save_np_array_with_string( COOR_DISP_RES , \
                                                   np_gps_site_i ,\
                                                  format_psvelo, \
                                                  fname, \
                                                  comment )

        # Up        
        COOR_DISP_RES_UP = COOR_DISP_OBS_UP - COOR_DISP_MODEL_UP
        COOR_DISP_RES_UP[:,:2] = COOR_DISP_OBS_UP[:,:2]
    
        fname = ("%s/displacement/cumulative/res/%04d_res_cum_disp_up.dat" % (odir, i))
        comment = (" observed up residual cumulative displacement at obs date #%04d %s" % ( i, np_obs_date_isoformat[i] ))
        
        pyacs.lib.utils.save_np_array_with_string( COOR_DISP_RES_UP , \
                                                   np_gps_site_i ,\
                                                  format_psvelo, \
                                                  fname, \
                                                  comment )



#     # SAVE WRMS IF RESIDUALS IN INFO DIR
# 
#     print("-- saving wrms in ", odir+'/info/wrms.dat' )
#     
#     fwrms = open(odir+'/info/wrms.dat','w+')
#     
#     fwrms.write("code  wrms_n (mm) wrms_e (mm) wrms_u (mm)\n")
#     fwrms.write("-----------------------------------------\n")
#     
#     wrms_n_mean = 0
#     wrms_e_mean = 0
#     wrms_u_mean = 0
#     
#     for code in sorted(H_code_wrms.keys()):
#         fwrms.write("%s %10.1lf %10.1lf %10.1lf\n" % (code,H_code_wrms[code][0], H_code_wrms[code][1], H_code_wrms[code][2] ))
# 
#         wrms_n_mean += H_code_wrms[code][0]
#         wrms_e_mean += H_code_wrms[code][1]
#         wrms_u_mean += H_code_wrms[code][2]
# 
#     wrms_n_mean = wrms_n_mean / len(list(H_code_wrms.keys()))
#     wrms_e_mean = wrms_e_mean / len(list(H_code_wrms.keys()))
#     wrms_u_mean = wrms_u_mean / len(list(H_code_wrms.keys()))
# 
#     fwrms.write("-----------------------------------------\n")
#     fwrms.write("mean %10.1lf %10.1lf %10.1lf\n" % (wrms_n_mean,wrms_e_mean,wrms_u_mean))
# 
#     fwrms.close()

    
    ###########################################################################
    # WRITE SUMMARY
    ###########################################################################


    M0=CSTF[-1]
    magnitude= 2./3.*(np.log10(M0)-9.1)

#    slip_max=np.max(CUMULATIVE_SLIP_TIME_STEP[:,2])
#    slip_min=np.min(CUMULATIVE_SLIP_TIME_STEP[:,2])

    fsumn=odir+'/summary/sum.dat'


    fsum=open(fsumn,'w')
    
    fsum.write('#--------------------------------------------------------------------------\n')
    fsum.write('# Inversion speed and memory usage: \n')
    fsum.write('#--------------------------------------------------------------------------\n')
    fsum.write("    number of parameters estimated   : %d \n" % SLIP.shape[0])
    fsum.write("    memory used in Gb Linux/Mac OS X : %.1lf / %.1lf \n" %  (H_inversion_info['memory_usage'], H_inversion_info['memory_usage']/1024) )
    fsum.write("    inversion start date & time      : %s \n" % time.strftime("%Y/%m/%d  %H:%M:%S", time.localtime(H_inversion_info['start_time']) ) )
    fsum.write("    inversion end   date & time      : %s \n" % time.strftime("%Y/%m/%d  %H:%M:%S", time.localtime(time.time()) ) )
    fsum.write("    time to build linear system (s)  : %.1lf \n" % H_inversion_info['time_ls'])
    fsum.write("    time nnls (s)                    : %.1lf \n" % H_inversion_info['time'])
    fsum.write("    nnls algorithm                   : %s \n" % H_inversion_info['nnls'])

    
    fsum.write('#--------------------------------------------------------------------------\n')
    fsum.write('# Inversion argument: \n')
    fsum.write('#--------------------------------------------------------------------------\n')
    fsum.write("    command line                     : %s \n" % (odir+'/info/command_line.dat'))
    fsum.write("    input npz file                   : %s \n" % H_inversion_info['input_npz'])

    fsum.write('#--------------------------------------------------------------------------\n')
    fsum.write('# Geometry information: \n')
    fsum.write('#--------------------------------------------------------------------------\n')
    fsum.write("    number of subfaults              : %d \n" % H_inversion_info['SGEOMETRY'].shape[0])

    fsum.write('#--------------------------------------------------------------------------\n')
    fsum.write('# Rake information: \n')
    fsum.write('#--------------------------------------------------------------------------\n')
    fsum.write("    rake_type (Euler or constant)    : %s \n" % H_inversion_info['rake_type'])
    fsum.write("    rake_value (Euler pole or value) : %s \n" % H_inversion_info['rake_value'])
    fsum.write("    rake_constraint (0=fixed, float) : %s \n" % H_inversion_info['rake_constraint'])
    fsum.write("    max_slip (0=from Euler,else user): %s \n" % H_inversion_info['max_slip'])

    fsum.write('#--------------------------------------------------------------------------\n')
    fsum.write('# Observation information: \n')
    fsum.write('#--------------------------------------------------------------------------\n')
    fsum.write("    time series are from             : %s \n" % H_inversion_info['input_obs'])
    fsum.write("    number of input GPS time series  : %d \n" %  H_inversion_info['n_site_input_ts'] )
    fsum.write("    number of input time series epoch: %d \n" %  H_inversion_info['n_date_input_ts'] )
    fsum.write("    number of GPS sites without data : %d \n" %  H_inversion_info['n_site_no_data'] )
    fsum.write("    number of input Green functions  : %d \n" %  H_inversion_info['n_site_in_green'] )
    fsum.write("    number of GPS sites used         : %d \n" %  H_inversion_info['T_OBS'].shape[1] )
    fsum.write("    number of epochs used            : %d \n" %  H_inversion_info['T_OBS'].shape[0] )
    fsum.write("    uncertainties E & N              : %s \n" %  H_inversion_info['h_uncertainty'] )
    fsum.write("    up component                     : %s \n" %  H_inversion_info['up'] )
    fsum.write("    uncertainty up rescaled by       : %s \n" %  H_inversion_info['s_up'] )

    fsum.write('#--------------------------------------------------------------------------\n')
    fsum.write('# Green tensor information: \n')
    fsum.write('#--------------------------------------------------------------------------\n')
    fsum.write("    Green tensor used for inversion  : %s \n" %  str( H_inversion_info['GREEN'].shape ) )


    fsum.write('#--------------------------------------------------------------------------\n')
    fsum.write('# Date information: \n')
    fsum.write('#--------------------------------------------------------------------------\n')
    fsum.write("    date rounding option             : %s \n" % ( H_inversion_info['rounding'] ))
    fsum.write("    number of model dates            : %d \n" % ( np_model_date_s.shape[0] ))
    fsum.write("    number of model time steps       : %d \n" % ( np_model_step_duration_s.shape[0] ))
    fsum.write("    model start date                 : %s doy %03d %13.8lf\n" % \
               ( np_model_date_isoformat[0] , at.datetime2dayno( np_model_datetime[0] )[1] , at.datetime2decyear(np_model_datetime[0])))
    fsum.write("    model end  date                  : %s doy %03d %13.8lf\n" % \
               ( np_model_date_isoformat[-1] , at.datetime2dayno( np_model_datetime[-1] )[1] , at.datetime2decyear(np_model_datetime[-1])))
    fsum.write("    median,shortest,longest time step: %.2lf %.2lf %.2lf days\n" % \
               (median_time_step_days,shortest_time_step_days,longest_time_step_days))
        
    
    fsum.write('#--------------------------------------------------------------------------\n')
    fsum.write('#   Inversion results: \n')
    fsum.write('#--------------------------------------------------------------------------\n')
    fsum.write("    Cumulative Moment (N.m)          : %8.2E \n" % M0)
    fsum.write("    Equivalent Mw                    : %8.1f \n" % magnitude)

    i_mr_max = np.argmax( STF ) 
    
    fsum.write("    max moment rate (N.m / day )     : %8.2E #%05d %s-%s\n" % \
               ( STF[i_mr_max], i_mr_max, np_model_date_isoformat[i_mr_max], np_model_date_isoformat[i_mr_max+1]))

    magnitude_mr= 2./3.*(np.log10(STF[ i_mr_max ])-9.1)
    fsum.write("    Equivalent Mw / day              : %8.1f \n" % magnitude_mr)

    i_cs_max = lindex = np.unravel_index(np.argmax(NORM_CUMULATIVE_SLIP_PER_TIME_STEP), NORM_CUMULATIVE_SLIP_PER_TIME_STEP.shape)
    slip_max = NORM_CUMULATIVE_SLIP_PER_TIME_STEP[ i_cs_max ]
    
    fsum.write("    cumulative slip max  (mm)        : %8.1f fault #%d %s\n" % ( slip_max , i_cs_max[1] , np_model_date_isoformat[ i_cs_max[0] ] ) )


    i_is_max = lindex = np.unravel_index(np.argmax(NORM_INC_SLIP_PER_TIME_STEP), NORM_INC_SLIP_PER_TIME_STEP.shape)
    inc_slip_max = NORM_INC_SLIP_PER_TIME_STEP[ i_is_max ]
    
    
    fsum.write("    incremental slip max  (mm)       : %8.1f fault #%04d time_step %04d %s -> %s\n" % \
               (inc_slip_max,i_is_max[1],i_is_max[0],np_model_date_isoformat[i_is_max[0]],np_model_date_isoformat[i_is_max[0]+1]))

    i_rs_max = lindex = np.unravel_index(np.argmax(NORM_RATE_SLIP_PER_TIME_STEP), NORM_RATE_SLIP_PER_TIME_STEP.shape)
    rate_slip_max = NORM_INC_SLIP_PER_TIME_STEP[ i_rs_max ]
    
    
    fsum.write("    slip rate max  (mm/day)          : %8.1f fault #%04d time_step %04d %s -> %s\n" % \
               (rate_slip_max,i_rs_max[1],i_rs_max[0],np_model_date_isoformat[i_rs_max[0]],np_model_date_isoformat[i_rs_max[0]+1]))
    
    
    # STATISTICS
    [rms_e, rms_n, rms_u] = np.std( RES[:,1:4] , axis=0 )
    [bias_e, bias_n, bias_u] = np.mean( RES[:,1:4] , axis=0 )
    
    num = np.sum( (RES[:,1:4]/RES[:,4:7])**2 , axis=0 )
    denom = np.sum( 1./RES[:,4:7]**2, axis=0  )
    [ wrms_e, wrms_n, wrms_u ] = np.sqrt( num / denom ) 

    [ mbias_e, mbias_n, mbias_u, mrms_e, mrms_n, mrms_u, mwrms_e, mwrms_n, mwrms_u ] = \
    np.mean( stats_site , axis=0 )
    
    i_max_wrms = np.argmax( stats_site[:,6]**2+stats_site[:,7] )
    
    chi2_obs = np.sum( num )
    
    
    fsum.write('#--------------------------------------------------------------------------\n')
    fsum.write('#   Inversion statistics: \n')
    fsum.write('#--------------------------------------------------------------------------\n')

    fsum.write("    rms  overall east  (mm)          : %8.1f\n" % ( rms_e ))
    fsum.write("    rms  overall north (mm)          : %8.1f\n" % ( rms_n ))
    fsum.write("    rms  overall up    (mm)          : %8.1f\n" % ( rms_u ))

    fsum.write("    wrms overall east  (mm)          : %8.1f\n" % ( wrms_e ))
    fsum.write("    wrms overall north (mm)          : %8.1f\n" % ( wrms_n ))
    fsum.write("    wrms overall up    (mm)          : %8.1f\n" % ( wrms_u ))

    fsum.write("    bias overall  east (mm)          : %8.1f\n" % bias_e)
    fsum.write("    bias overall north (mm)          : %8.1f\n" % bias_n)
    fsum.write("    bias overall    up (mm)          : %8.1f\n" % bias_u)

    fsum.write("    wrms mean east  (mm)             : %8.1f\n" % (mwrms_e))
    fsum.write("    wrms mean north (mm)             : %8.1f\n" % (mwrms_n))
    fsum.write("    wrms mean up    (mm)             : %8.1f\n" % mwrms_u)

    fsum.write("    bias mean  east (mm)             : %8.1f\n" % mbias_e)
    fsum.write("    bias mean north (mm)             : %8.1f\n" % mbias_n)
    fsum.write("    bias mean    up (mm)             : %8.1f\n" % mbias_u)


    fsum.write("    wrms max  east  (mm)             : %8.1f %s\n" % (stats_site[i_max_wrms,6] , H_inversion_info['np_gps_site'][i_max_wrms]))
    fsum.write("    wrms max north  (mm)             : %8.1f %s\n" % (stats_site[i_max_wrms,7] , H_inversion_info['np_gps_site'][i_max_wrms]))
    fsum.write("    wrms max    up  (mm)             : %8.1f %s\n" % (stats_site[i_max_wrms,8] , H_inversion_info['np_gps_site'][i_max_wrms]))

    fsum.write("    chi2 obs                         : %8.1f\n" %  chi2_obs )
    
    n_component = 2
    if H_inversion_info['up'] != 'not used':
        n_component = 3 
    
    fsum.write("    number of observations           : %8d\n" % RES.shape[0] * n_component )

    fsum.write("    reduced chi2 obs                 : %8.2f\n" % np.sqrt(chi2_obs / ( RES.shape[0] * n_component )  ))
#    fsum.write("    chi2 regularization              : %8.1f\n" % chi2_regularization)
#    fsum.write("    number of constraints            : %8d\n" % ( SLIP.shape[0] / SGEOMETRY.shape[0] * nfaults) )
#    fsum.write("    reduced chi2 regularization      : %8.2f\n" % np.sqrt(chi2_regularization / float(SLIP.shape[0] / SGEOMETRY.shape[0] * nfaults) ))
#    fsum.write("    chi2 all                         : %8.1f\n" % (chi2_obs+chi2_regularization) )
#    fsum.write("    reduced chi2 all                 : %8.1f\n" % np.sqrt( (chi2_obs+chi2_regularization) / ( float(n_obs) + float(SLIP.shape[0] / SGEOMETRY.shape[0] * nfaults) ) ) )

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
    pyeq.log.model2shp_gmt.model2shp_gmt(GEOMETRY, 'tde', lslip_dat, out_dir_shp=odir + '/shapefile/slip_cumulative', out_dir_gmt=odir + '/gmt/slip_cumulative', verbose=verbose)

    lslip_dat=glob.glob(odir+"/slip/incremental/*_delta_slip.dat")
    pyeq.log.model2shp_gmt.model2shp_gmt(GEOMETRY, 'tde', lslip_dat, out_dir_shp=odir + '/shapefile/slip_incremental', out_dir_gmt=odir + '/gmt/slip_incremental', verbose=verbose)


#     sys.exit()
# 
#     for slip_dat in lslip_dat:
#     
#         SLIP=np.genfromtxt(slip_dat,comments='#')
# 
#         lfaults=[]
#         lrecord=[]
#         for i in np.arange(GEOMETRY.shape[0]):
#             ( rdis_long, rdis_lat, rdis_depth, rdis_length, rdis_width,  rdis_area, ratio_rdis_tdis,strike,dip,\
#              centroid_long, centroid_lat, centroid_depth,\
#              tdis_long1,tdis_lat1,tdis_depth1,tdis_long2,tdis_lat2,tdis_depth2,tdis_long3,tdis_lat3,tdis_depth3,tdis_area)\
#              =GEOMETRY[i,:]
#             depth=rdis_depth
#         
#     
#             if TRIANGLE:
#                 lfaults.append([ [tdis_long1,tdis_lat1], [tdis_long2,tdis_lat2], [tdis_long3,tdis_lat3] ])
#                 lrecord.append([i,centroid_depth,SLIP[i,-1],0])
#         
#             else:
#             
#                 # creates a dislocation object
#                 disloc=DL.Dislocation(i,rdis_long,rdis_lat,rdis_depth/one_degree, strike, dip, rdis_length/one_degree, rdis_width/one_degree,rdis_area, 0, 0)
#                 # get the corners
#                 (X1,X2,X3,X4)=disloc.corners()
#                 lfaults.append([ [X1[0],X1[1]], [X2[0],X2[1]], [X3[0],X3[1]], [X4[0],X4[1]], [X1[0],X1[1]] ])
#                 lrecord.append([i,depth,SLIP[i,-1],0])
#     
#         print("-- ",len(lfaults)," polygons read")
# 
#         ###################################################################
#         # WRITES GMT PSXY FILES
#         # This file can be then plotted with
#         #  psxy lima_simple_sigma_010_dc_050_m0_000_sol_coupling.gmt -V -R-81/-72/-16/-7 -JM14 -L -m -W0.2/0  -Clut_coupling.cpt > test.ps
#         # triangles have not been tested yet
#         ###################################################################
#         
#         gmtfile=os.path.abspath(slip_dat).split('.')[0]+'.gmt'
# 
#         print(("-- saving gmt file %s " % gmtfile.split('/')[-1]))
#         f=open(gmtfile,'w')
#                 
#         for i in np.arange(len(lfaults)):
#             [index,depth,slip,rake]=lrecord[i]
#             f.write('> -Z%.3lf\n'%slip)
#             fault=lfaults[i]
#             for xy in fault:
#                 f.write("%10.3lf %10.3lf\n" %(xy[0],xy[1]))
#         f.write('>\n')
#         f.close()
# 
#         ###################################################################
#         # WRITES SHAPEFILES
#         ###################################################################
# 
#         shp_file=os.path.abspath(slip_dat).split('.')[0]
#  
#        ###################################################################
#         # INITIALIZE SHAPEFILE
#         ###################################################################
#         
#         # Make a polygon shapefile
#         w = shapefile.Writer( shp_file ,shapeType=shapefile.POLYGON)
#         w.field('ID','I','40')
#         w.field('i_subfault','F','40')
#         w.field('depth_top_disloc','F','40')
#         w.field('slip','F','40')
#         w.field('rake','F','40')
# 
#         
#         ###################################################################
#         # LOOP ON FAULTS
#         ###################################################################
#         
#         for i in np.arange(len(lfaults)):
#             fault=lfaults[i]
#             record=lrecord[i]
#             w.poly([fault])
#             [index,depth,slip,rake]=lrecord[i]
#             w.record(str(i),index,depth,slip,rake)
#             i=i+1
#         
#         ###################################################################
#         # SAVE SHAPEFILE
#         ###################################################################
#         
#         print(("-- saving shapefile %s " % (shp_file.split('/')[-1]+'.shp')))
#         w.close()
#         
#         ###################################################################
#         # SAVE .PRJ FILE
#         ###################################################################
#         
#         prj = open("%s.prj" % shp_file, "w") 
#         epsg = 'GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563]],PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433]]' 
#         prj.write(epsg) 
#         prj.close()


###############################################################################
### END  print_kinematic_inversion_results4 ###################################
###############################################################################
