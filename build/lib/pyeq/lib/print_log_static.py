###############################################################################
def print_static_inversion_results(M, GG, d, INPUT_NPZ , RAKE, MAX_SLIP, out_dir , SQRT_Wm, SQRT_Wd, m0, H_inversion_info):
###############################################################################

    """
    Print static inversion results
    """


    ###################################################################
    # IMPORT
    ###################################################################

    import numpy as np
    import os
    from pyeq.lib import lib_inversion
    import pyacs.lib.faultslip


    ###################################################################
    # INITIALIZATION
    ###################################################################

    name = out_dir
    
    big = H_inversion_info['big'] 

    if H_inversion_info['insar']:
        # The last 3 parameters are the coefficient of the plane correction for InSAR shift by an amount of 'big'
        
        [a,b,c] = M[-3:] - big
        M[-3:] = np.array([a,b,c])
        
        SLIP = M[:-3] 
        G = GG[:,0:-3]
        
        m0 = m0[:-3]
        m0_a_b_c = m0[-3:] - big
        
        # Get SQRT_Wm
        
        Wmm = np.dot(SQRT_Wm.T,SQRT_Wm) # we get CM^{-1}, which should be a block diagonal matrix 
        
        sigma_a_b_c = 1. / np.sqrt(np.diag( Wmm[-3:,-3:] )) 
    
    else:
        SLIP = M
        G = GG
        [a,b,c] = [0.,0.,0.]

    mu=3.0E10 # shear modulus in Pa

    ###################################################################
    # PYTHON 3.6 READING INPUT NPZ FILE
    ###################################################################
    
    INPUT = np.load(INPUT_NPZ)
    
    GEOMETRY    = INPUT["GEOMETRY"]
    Dm          = INPUT["Dm"]
    GREEN       = INPUT["GREEN"]
    GREEN_UP    = INPUT["GREEN_UP"]
    OBS         = INPUT["OBS"]
    NAME_OBS    = INPUT["NAME_OBS"]
    OBS_UP      = INPUT["OBS_UP"]
    NAME_OBS_UP = INPUT["NAME_OBS_UP"]
    GREEN_INSAR = INPUT["GREEN_INSAR"]
    OBS_INSAR   = INPUT["OBS_INSAR"]

    if np.count_nonzero( OBS[:,4:6] ) == 0:
        print("-- warning: data have no uncertainty. Setting to 1.")
        OBS[:, 4:6] = 1.


    # converts array to recarray
    
    names=['rdis_long','rdis_lat','rdis_depth','rdis_length','rdis_width',\
           'rdis_area','ratio_rdis_tdis','strike','dip',\
           'centroid_long','centroid_lat','centroid_depth',\
           'tdis_long1','tdis_lat1','tdis_depth1',\
           'tdis_long2','tdis_lat2','tdis_depth2',\
           'tdis_long3','tdis_lat3','tdis_depth3','tdis_area']
    
    SGEOMETRY=lib_inversion.numpy_array_2_numpy_recarray(GEOMETRY,names)

    
    
    ###################################################################
    # END PYTHON 3.6 READING INPUT NPZ FILE
    ###################################################################
    

    nfaults = GEOMETRY.shape[0]
    n_fields_geometry = len(GEOMETRY[0])
    
    
    ngps=OBS.shape[0]
    ngps_up=OBS_UP.shape[0]
    ninsar=OBS_INSAR.shape[0]

    # nn=1 principal rake only, nn=2 rake_principal and rake_conjugate
    
    if SLIP.shape[0] == 2 * nfaults:
        nn = 2
    else:
        nn = 1  

    ###################################################################
    def save_np_array_with_string_2_psvelo_format(A,S,fmt,name):
    ###################################################################
        """
        saves a numpy array with the last column being a 1D numpy array of dtype string
        """
        # Create a structured numpy array
        nrow=A.shape[0]
        ncol=A.shape[1]+1
        
        Snp = np.zeros(nrow, dtype=[('lon', 'float64'), ('lat', 'float64'), ('ve', 'float64'), ('vn', 'float64'), ('sve', 'float64'), ('svn', 'float64'), ('sven', 'float64'),('code', 'U255'), ])
        lfield=list(Snp.dtype.names)
        lfield.pop(-1)
        
        i=0
        for field in lfield:
            Snp[field]=A[:,i] 
            i=i+1
           
        Snp['code'] = S

        np.savetxt(name, Snp, fmt=fmt)

    ###################################################################

    print('-- Writing solutions for inversion in directory',out_dir)

    ###################################################################
    # FILE NAMES FOR OUTPUT RESULTS
    ###################################################################


    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    else:
        import shutil
        shutil.rmtree(out_dir)        
        os.mkdir(out_dir)

    fslip     = out_dir+'/'+name+'_sol_slip.dat'
    fcoupling = out_dir+'/'+name+'_sol_coupling.dat'
    fpred     = out_dir+'/'+name+'_model.dat'
    fpred_up  = out_dir+'/'+name+'_model_up.dat'
    fpred_insar  = out_dir+'/'+name+'_model_insar.dat'
    fres      = out_dir+'/'+name+'_residuals.dat'
    fres_up   = out_dir+'/'+name+'_residuals_up.dat'
    fres_insar= out_dir+'/'+name+'_residuals_insar.dat'
    fsumn     = out_dir+'/'+name+'_sum.dat'
    fslip_dir = out_dir+'/'+name+'_slip_dir.dat'
    f_geometry= out_dir+'/'+name+'_geometry.dat'
    f_obs     = out_dir+'/'+name+'_obs.dat'
    f_obsup   = out_dir+'/'+name+'_obs_up.dat'


    ###########################################################################    
    # MODEL
    ###########################################################################    

    import pyeq.lib.lib_inversion
    GEOMETRY = pyeq.lib.lib_inversion.numpy_recarray_2_numpy_array(SGEOMETRY)
    
    MODEL=np.zeros((nfaults,GEOMETRY.shape[1]+5)) # geometry_fields,rake1,slip1,rake2,slip2,slip_total
    MODEL[:,:n_fields_geometry]=GEOMETRY
    MODEL[:,n_fields_geometry]=RAKE # rake principal
    MODEL[:,n_fields_geometry+1]=SLIP[:nfaults] # rake principal
    MODEL[:,n_fields_geometry+2]=RAKE+90.0 # rake conjugate
    
    if SLIP.shape[0] == 2 * nfaults:
        # variable rake case
        MODEL[:,n_fields_geometry+3]=SLIP[nfaults:] # rake conjugate case
    else:
        MODEL[:,n_fields_geometry+3]=SLIP*0.0 # rake principal only case
        
    MODEL[:,n_fields_geometry+4]=np.sqrt(MODEL[:,n_fields_geometry+1]**2+MODEL[:,n_fields_geometry+3]**2)
    
    header_model="rdis_long   rdis_lat rdis_depth     rdis_length rdis_width  rdis_area ratio_rdis_tdis        strike       dip\
      centroid_long centroid_lat centroid_depth     tdis_long1\
  tdis_lat1 tdis_depth1     tdis_long2  tdis_lat2 tdis_depth2     tdis_long3  tdis_lat3 tdis_depth3      tdis_area\
 rake_1 slip_1 rake_2  slip2   slip"

    header_model="rdis_long   rdis_lat rdis_depth     rdis_length rdis_width  rdis_area ratio_rdis_tdis        strike       dip\
      centroid_long centroid_lat centroid_depth     tdis_long1\
  tdis_lat1 tdis_depth1     tdis_long2  tdis_lat2 tdis_depth2     tdis_long3  tdis_lat3 tdis_depth3      tdis_area\
 rake_1 slip_1 rake_2  slip2  ISC(%)"

    format_model=\
    '%10.5lf %10.5lf      %6.2lf \
         %6.2lf     %6.2lf    %6.2lf          %6.2lf\
        %6.2lf %10.2lf \
       %10.5lf   %10.5lf         %6.2lf \
    %10.5lf %10.5lf      %6.2lf \
    %10.5lf %10.5lf      %6.2lf \
    %10.5lf %10.5lf      %6.2lf \
       %6.2lf   %6.2lf %6.2lf %6.2lf %6.2lf %6.2lf'

    np.savetxt(fslip,MODEL,fmt=format_model,header=header_model)
    
    ###########################################################################
    # MOMENT
    ###########################################################################

    SLIP_POTENCY=MODEL[:,-1]*1.0E-3 * MODEL[:,5]*1.0E6

    M0=mu*np.sum(SLIP_POTENCY)
    
    magnitude= 2./3.*(np.log10(M0)-9.1)
    
    ###########################################################################
    # COUPLING
    ###########################################################################
    
    COUPLING=np.copy(MODEL)
    COUPLING[:,n_fields_geometry+4]=MODEL[:,n_fields_geometry+4]/MAX_SLIP*100.0
    np.savetxt(fcoupling,COUPLING,fmt=format_model,header=header_model)
    
    ###########################################################################
    # MODEL PREDICTIONS
    ###########################################################################

    PREDICTION_ALL=np.dot(GG,M)

    # horizontal GPS
    
    PREDICTION_EAST=PREDICTION_ALL[0:ngps]
    PREDICTION_NORTH=PREDICTION_ALL[ngps:2*ngps]

    # fills the horizontal model prediction GPS file for writing in gmt psvelo format
    PREDICTION_GPS=np.copy(OBS)
    PREDICTION_GPS[:,2]=PREDICTION_EAST
    PREDICTION_GPS[:,3]=PREDICTION_NORTH

    PREDICTION_GPS[:,4]=0.0
    PREDICTION_GPS[:,5]=0.0
    PREDICTION_GPS[:,6]=0.0
    
    save_np_array_with_string_2_psvelo_format(PREDICTION_GPS,NAME_OBS,'%10.5lf  %10.5lf %10.5lf  %10.5lf %3.1lf %3.1lf %3.1lf %s',fpred)

    # vertical GPS
    if ngps_up>0:

        PREDICTION_UP=PREDICTION_ALL[2*ngps:2*ngps+ngps_up]
        PREDICTION_GPS_UP=np.copy(OBS_UP)
        PREDICTION_GPS_UP[:,2]=PREDICTION_UP
        PREDICTION_GPS_UP[:,3]=0.0
        
        np.savetxt(fpred_up, PREDICTION_GPS_UP, fmt='%10.5lf  %10.5lf %10.5lf  %10.5lf', delimiter=' ')
    
    # InSAR

    if H_inversion_info['insar']:

        MODEL_PREDICTION_INSAR=PREDICTION_ALL[2*ngps+ngps_up:2*ngps+ngps_up+ninsar]
        PREDICTION_INSAR=np.copy(OBS_INSAR[:,0:4])

        PREDICTION_INSAR[:,2] = MODEL_PREDICTION_INSAR
        PREDICTION_INSAR[:,3]=0.0
        
        np.savetxt(fpred_insar, PREDICTION_INSAR, fmt='%10.5lf  %10.5lf %10.5lf  %10.5lf', delimiter=' ')
    
    ###########################################################################
    # MODEL RESIDUALS
    ###########################################################################
    
    # horizontal GPS
    
    RESIDUALS_GPS=np.copy(OBS)
    RESIDUALS_GPS[:,2]=OBS[:,2]-PREDICTION_GPS[:,2]
    RESIDUALS_GPS[:,3]=OBS[:,3]-PREDICTION_GPS[:,3]

    save_np_array_with_string_2_psvelo_format(RESIDUALS_GPS,NAME_OBS,'%10.5lf  %10.5lf %10.5lf  %10.5lf %3.1lf %3.1lf %3.1lf %s',fres)

    # vertical GPS
    if ngps_up>0:
        RESIDUALS_GPS_UP=np.copy(OBS_UP)
        RESIDUALS_GPS_UP[:,2]=OBS_UP[:,2]-PREDICTION_GPS_UP[:,2]
        np.savetxt(fres_up, RESIDUALS_GPS_UP, fmt='%10.5lf  %10.5lf %10.5lf  %10.5lf', delimiter=' ')

    # insar
    
    if H_inversion_info['insar']:
        RESIDUALS_INSAR = np.copy(OBS_INSAR[:,0:4])
        RESIDUALS_INSAR[:,2]=OBS_INSAR[:,2] - PREDICTION_INSAR[:,2]
        np.savetxt(fres_insar, RESIDUALS_INSAR, fmt='%10.5lf  %10.5lf %10.5lf  %10.5lf', delimiter=' ')
        
    
    ###########################################################################
    # STATISTICS SLIP MIN & MAX
    ###########################################################################
    
    slip_max=np.max(MODEL[:,n_fields_geometry+4])
    slip_max_rake_principal=np.max(MODEL[:,n_fields_geometry+1])
    slip_max_rake_conjugate=np.sqrt(np.max(MODEL[:,n_fields_geometry+3])**2)

    slip_min=np.min(MODEL[:,n_fields_geometry+4])
    slip_min_rake_principal=np.min(MODEL[:,n_fields_geometry+1])
    slip_min_rake_conjugate=np.sqrt(np.min(MODEL[:,n_fields_geometry+3])**2)


    ###########################################################################
    # STATISTICS RMS (NO WEIGHT ACCOUNTED FOR)
    ###########################################################################
    
    # rms horizontal GPS
    
    rms_obs  = np.sqrt( np.sum( RESIDUALS_GPS[:,2]**2 + RESIDUALS_GPS[:,3]**2 ) / (2*ngps) ) 
    
    # rms vertical GPS
    if ngps_up>0:
        rms_obs_up= np.sqrt( np.sum( (RESIDUALS_GPS_UP[:,2]**2) ) / ngps_up )
    else:
        rms_obs_up = 0.

    # rms insar
    if H_inversion_info['insar']:
        rms_insar = np.sqrt( np.sum( RESIDUALS_INSAR[:,2]**2) / ninsar )
    else:
        rms_insar = 0.

    # overall rms
    
    rms_obs_all = np.sqrt( ( rms_obs**2 * (2*ngps) + rms_obs_up**2 * ngps_up + rms_insar**2 * ninsar ) / ( 2*ngps + ngps_up + ninsar ) )

    
    ###########################################################################
    # STATISTICS CHI2
    ###########################################################################
    
    # chi2 horizontal GPS
    
    chi2_obs  = np.sum((RESIDUALS_GPS[:,2]/RESIDUALS_GPS[:,4])**2 + (RESIDUALS_GPS[:,3]/RESIDUALS_GPS[:,5])**2)
    
    
    # chi2 vertical GPS
    if ngps_up>0:
        chi2_obs_up= np.sum((RESIDUALS_GPS_UP[:,2]/RESIDUALS_GPS_UP[:,3])**2)
    else:
        chi2_obs_up = 0.

    # chi2 insar
    if H_inversion_info['insar']:
        chi2_insar= np.sum((RESIDUALS_INSAR[:,2]/RESIDUALS_INSAR[:,3])**2)
    else:
        chi2_insar = 0.

    # overall chi2
    
    chi2_obs_all = chi2_obs
    
    if ngps_up>0:
        chi2_obs_all = chi2_obs_all + chi2_obs_up
        
    if H_inversion_info['insar']:
        chi2_obs_all = chi2_obs_all + chi2_insar

    ###########################################################################
    # STATISTICS REDUCED CHI2
    ###########################################################################

    red_chi2_obs = np.sqrt( chi2_obs / ( 2*ngps ) )
    
    if ngps_up > 0 :
        red_chi2_up = np.sqrt( chi2_obs_up / ngps_up )
    else:
        red_chi2_up = 0.

    if H_inversion_info['insar']:
        red_chi2_insar = np.sqrt( chi2_insar / ninsar )
    else:
        red_chi2_insar = 0.

    red_chi2_all = np.sqrt( chi2_obs_all / (2*ngps+ngps_up+ninsar) )

    ###########################################################################
    # STATISTICS WRMS
    ###########################################################################

    # wrms horizontal GPS        
    denom_obs = np.sum(1./RESIDUALS_GPS[:,4]**2) + np.sum(1./RESIDUALS_GPS[:,5]**2)
    wrms_obs=np.sqrt( chi2_obs /  denom_obs )
    
    # wrms vertical GPS
    if ngps_up>0:
        denom_obs_up = np.sum(1./RESIDUALS_GPS_UP[:,3]**2)
        wrms_obs_up  = np.sqrt( chi2_obs_up / denom_obs_up )
    else:
        wrms_obs_up  = 0.0
        denom_obs_up = 0.0
    
    # wrms insar
    if H_inversion_info['insar']:
        denom_insar = np.sum(1./RESIDUALS_INSAR[:,3]**2)
        wrms_insar  = np.sqrt( chi2_insar / denom_insar )
    else:
        wrms_insar  = 0.0
        denom_insar = 0.0

    # wrms all obs
    
    wrms_all = np.sqrt( chi2_obs_all / (denom_obs + denom_obs_up + denom_insar) )

    ###########################################################################
    # STATISTICS VARIANCE REDUCTION
    ###########################################################################
    
    # variance reduction horizontal GPS
    
    chi2_obs_prefit  = np.sum((OBS[:,2]/OBS[:,4])**2 + (OBS[:,3]/OBS[:,5])**2)
    redvar_obs = ( chi2_obs_prefit - chi2_obs ) / chi2_obs_prefit * 100.0

    # variance reduction vertical GPS
    if ngps_up>0:
        chi2_obs_up_prefit  = np.sum((OBS_UP[:,2]/OBS_UP[:,3])**2)
        redvar_obs_up = ( chi2_obs_up_prefit - chi2_obs_up ) / chi2_obs_up_prefit * 100.0
    else:
        redvar_obs_up = 0.0
        chi2_obs_up_prefit  = 0.0
             
    # variance reduction insar
    if H_inversion_info['insar']:
        chi2_insar_prefit  = np.sum((OBS_INSAR[:,2]/OBS_INSAR[:,3])**2)
        redvar_insar = ( chi2_insar_prefit - chi2_insar ) / chi2_insar_prefit * 100.0
    else:
        redvar_insar = 0.0
        chi2_insar_prefit  = 0.0
    
    # variance reduction all
    
    chi2_prefit  = chi2_obs_prefit + chi2_obs_up_prefit +chi2_insar_prefit
    chi2_postfit = chi2_obs 
    if ngps_up > 0:
        chi2_postfit = chi2_postfit + chi2_obs_up
    if H_inversion_info['insar']:
        chi2_postfit = chi2_postfit + chi2_insar
         
    redvar_all = ( chi2_prefit - chi2_postfit ) / chi2_prefit * 100.0


    ###########################################################################
    # STATISTICS FULL CHI2
    # This is now removed - I'll see in future whether it is useful or not
    ###########################################################################

    # complete chi2 S(m)=(Gm-d).T Cd-1 (Gm-d) + (m-m0).T Cm-1 (m-m0)
    # S(m)=(Gm-d).T SQRT_Wd.T SQRT_Wd (Gm-d) + (m-m0).T SQRT_Wm.T SQRT_Wm (m-m0)
    
#     SQRT_Wd_Gm_d=np.dot(SQRT_Wd,np.dot(G,SLIP)-d)
# 
#     if H_inversion_info['insar']:
#         SQRT_Wm_m_m0=np.dot(SQRT_Wm[0:-3,0:-3],SLIP-m0)
#     else:
#         SQRT_Wm_m_m0=np.dot(SQRT_Wm,SLIP-m0)
#     
#     chi2_1=np.dot(SQRT_Wd_Gm_d.T,SQRT_Wd_Gm_d)
#     chi2_2=np.dot(SQRT_Wm_m_m0.T,SQRT_Wm_m_m0)
#     
#     chi2_all=0.5*(chi2_1+chi2_2)

    ###########################################################################
    # STATISTICS MEAN BIAS
    ###########################################################################
    
    bias_east=np.mean(RESIDUALS_GPS[:,2])
    bias_north=np.mean(RESIDUALS_GPS[:,3])
    if ngps_up>0:
        bias_up=np.mean(RESIDUALS_GPS_UP[:,2])
    else:
        bias_up=0.0

    if H_inversion_info['insar']:
        bias_insar=np.mean(RESIDUALS_INSAR[:,2])
    else:
        bias_insar=0.0

    ###########################################################################
    # WORST RESIDUALS
    ###########################################################################
    
    # hor. GPS
    id_worst_gps = np.argmax( RESIDUALS_GPS[:,2]**2 + RESIDUALS_GPS[:,3]**2 )
    worst_gps, wgps_east, wgps_north = NAME_OBS[id_worst_gps],RESIDUALS_GPS[id_worst_gps,2],RESIDUALS_GPS[id_worst_gps,3]
    # up GPS
    if ngps_up > 0:
        id_worst_gps_up = np.argmax( RESIDUALS_GPS_UP[:,2]**2)
        worst_gps_up, wgps_up = NAME_OBS_UP[id_worst_gps_up],RESIDUALS_GPS_UP[id_worst_gps_up,2]
    else:
        worst_gps_up, wgps_up = 'NA',0.
        
    # insar
    if H_inversion_info['insar']:
        id_worst_insar = np.argmax( RESIDUALS_INSAR[:,2]**2)
        [worst_res_insar_lon,worst_res_insar_lat,worst_res_insar] = RESIDUALS_INSAR[id_worst_insar,:3]
    else:
        [worst_res_insar_lon,worst_res_insar_lat,worst_res_insar] = [0.,0.,0.]
        
    ###########################################################################
    # LARGEST SLIP
    ###########################################################################
    
    id_largest_slip = np.argmax( MODEL[:,-1] )
    [lon_ls,lat_ls] = MODEL[id_largest_slip,:2]
    
    ###########################################################################
    # SUMARY FILE
    ###########################################################################

    
    fsum=open(fsumn,'w')
    
    fsum.write('#--------------------------------------------------------------------------\n')
    fsum.write('# Input file: \n')
    fsum.write('#--------------------------------------------------------------------------\n')
    fsum.write("    command line                     : %s \n" % H_inversion_info['cmd_line'])
    fsum.write("    input npy file                   : %s \n" % H_inversion_info['input_npz'])
    fsum.write('#--------------------------------------------------------------------------\n')
    fsum.write('# Inversion type: \n')
    fsum.write('#--------------------------------------------------------------------------\n')
    fsum.write("    rake_type (Euler or constant)    : %s \n" % H_inversion_info['rake_type'])
    fsum.write("    rake_value (Euler pole or value) : %s \n" % H_inversion_info['rake_value'])
    fsum.write("    rake_constraint (0=fixed, float) : %s \n" % H_inversion_info['rake_constraint'])
    fsum.write("    max_slip (0=from Euler,else user): %s \n" % H_inversion_info['max_slip'])
    fsum.write('#--------------------------------------------------------------------------\n')
    fsum.write('# Model and observation parameters: \n')
    fsum.write('#--------------------------------------------------------------------------\n')
    fsum.write("    number of subfaults              : %d \n" % SGEOMETRY.shape[0])
    fsum.write("    number of parameters estimated   : %d \n" % GG.shape[1])
    fsum.write("    number of total observations     : %d \n" % GG.shape[0])
    fsum.write("    number of GPS horizontal obs.    : %d \n" % OBS.shape[0])
    fsum.write("    number of GPS vertical   obs.    : %d \n" % OBS_UP.shape[0])
    fsum.write("    number of InSAR          obs.    : %d \n" % ninsar )
    fsum.write("    scaling factor for h.GPS unc.    : %4.1lf \n" % H_inversion_info['wh'] )
    fsum.write("    scaling factor for u.GPS unc.    : %4.1lf \n" % H_inversion_info['wu'] )
    fsum.write("    scaling factor for InSAR         : %4.1lf \n" % H_inversion_info['wi'] )
    fsum.write("    constraints on InSAR corrections : %s \n" % H_inversion_info['ci'] )
    fsum.write('#--------------------------------------------------------------------------\n')
    fsum.write("    inversion duration (s)           : %.1lf \n" % H_inversion_info['time'])
    fsum.write('#--------------------------------------------------------------------------\n')
    fsum.write('#--------------------------------------------------------------------------\n')
    fsum.write('#   Inversion results: \n')
    fsum.write('#--------------------------------------------------------------------------\n')
    fsum.write("    Moment (N.m)                     : %8.2E \n" % M0)
    fsum.write("    Equivalent Mw                    : %8.2f \n" % magnitude)
    fsum.write("    slip max (unit of input gps file): %8.1f (%4.2f,%4.2f)\n" % (np.max(MODEL[:,n_fields_geometry+4]),lon_ls,lat_ls))
    fsum.write("    slip min (unit of input gps file): %8.1f\n" % np.min(MODEL[:,n_fields_geometry+4]))
    fsum.write("    InSAR correction a*x+b*y+c       : %5.1f %5.1f %5.1f\n" %  (a,b,c) )
    fsum.write('#--------------------------------------------------------------------------\n')
    fsum.write('#   Inversion statistics: \n')
    fsum.write('#--------------------------------------------------------------------------\n')
    fsum.write("    rms h_gps (unit input gps file)  : %8.2f \n" % rms_obs)
    fsum.write("    rms up (unit of input gps file)  : %8.2f \n" % rms_obs_up)
    fsum.write("    rms InSAR                        : %8.2f \n" % rms_insar)
    fsum.write("    rms all                          : %8.2f \n" % rms_obs_all)
    fsum.write('#--------------------------------------------------------------------------\n')
    fsum.write("    wrms h_gps (unit input gps file) : %8.2f \n" % wrms_obs)
    fsum.write("    wrms up (unit of input gps file) : %8.2f \n" % wrms_obs_up)
    fsum.write("    wrms InSAR                       : %8.2f \n" % wrms_insar)
    fsum.write("    wrms all                         : %8.2f \n" % wrms_all)
    fsum.write('#--------------------------------------------------------------------------\n')
    fsum.write("    chi2 prefit gps horizontal       : %8.1f \n" % chi2_obs_prefit)
    fsum.write("    chi2 prefit gps vertical         : %8.1f \n" % chi2_obs_up_prefit)
    fsum.write("    chi2 prefit InSAR                : %8.1f \n" % chi2_insar_prefit)
    fsum.write("    chi2 prefit all  (only diagonal) : %8.1f \n" % chi2_prefit)
    fsum.write('#--------------------------------------------------------------------------\n')
    fsum.write("    chi2 postfit gps horizontal      : %8.1f \n" % chi2_obs)
    fsum.write("    chi2 postfit gps vertical        : %8.1f \n" % chi2_obs_up)
    fsum.write("    chi2 postfit InSAR               : %8.1f \n" % chi2_insar)
    fsum.write("    chi2 all obs (only diagonal)     : %8.1f \n" % chi2_obs_all)
    fsum.write('#--------------------------------------------------------------------------\n')
    fsum.write("    reduced chi2 hor. GPS            : %8.1f \n" % red_chi2_obs )
    fsum.write("    reduced chi2 up   GPS            : %8.1f \n" % red_chi2_up )
    fsum.write("    reduced chi2 InSAR               : %8.1f \n" % red_chi2_insar )
    fsum.write("    reduced chi2 all obs             : %8.1f \n" % red_chi2_all )
    fsum.write('#--------------------------------------------------------------------------\n')
    fsum.write("    variance reduction hor. GPS      : %8.1f %%\n" % redvar_obs)
    fsum.write("    variance reduction up   GPS      : %8.1f %%\n" % redvar_obs_up)
    fsum.write("    variance reduction InSAR         : %8.1f %%\n" % redvar_insar)
    fsum.write("    variance reduction all obs       : %8.1f %%\n" % redvar_all)
    fsum.write('#--------------------------------------------------------------------------\n')
    fsum.write("    bias east                        : %8.1f \n" % bias_east)
    fsum.write("    bias north                       : %8.1f \n" % bias_north)
    fsum.write("    bias up                          : %8.1f \n" % bias_up)
    fsum.write("    bias InSAR                       : %8.1f \n" % bias_insar)
    fsum.write('#--------------------------------------------------------------------------\n')
    fsum.write("    worst hor. GPS                   : %s %5.1f %5.1f\n" %  (worst_gps , wgps_east, wgps_north) )
    fsum.write("    worst up. GPS                    : %s %5.1f \n" %  (worst_gps_up , wgps_up) )
    fsum.write("    worst InSAR                      : %5.1f at (%4.2f,%4.2f)\n" %  (worst_res_insar,worst_res_insar_lon,worst_res_insar_lat) )
    fsum.write('#--------------------------------------------------------------------------\n')
    fsum.close()
    
        
    # print results on screen
    
    f = open(fsumn, "r")
    text = f.read()
    print(text)
    f.close()

    # slip direction file
    
    import pyacs.lib.coordinates as C
    import pyacs.lib

    MODEL[:,:n_fields_geometry]=GEOMETRY
    MODEL[:,n_fields_geometry]=RAKE # rake principal
    MODEL[:,n_fields_geometry+1]=SLIP[:nfaults] # rake principal
    MODEL[:,n_fields_geometry+2]=RAKE+90.0 # rake conjugate

    SLIP_DIR=np.zeros((MODEL.shape[0],8))
    
    for i in range(MODEL.shape[0]):
        # centroid
        SLIP_DIR[i,0]=MODEL[i,9]
        SLIP_DIR[i,1]=MODEL[i,10]
        # strike & dip
        strike=MODEL[i,7]
        dip=MODEL[i,8]
        # rake_1 & rake_2
        rake_1=MODEL[i,22]
        rake_2=MODEL[i,24]
        # slip_1 & slip_2
        slip_1=MODEL[i,23]
        slip_2=MODEL[i,25]
        # get direction
        direction_1=pyacs.lib.faultslip.strike_dip_rake_to_dir(strike, dip, rake_1)
        (e1,n1)=C.azimuth_to_en(direction_1)
        direction_2=pyacs.lib.faultslip.strike_dip_rake_to_dir(strike, dip, rake_2)
        (e2,n2)=C.azimuth_to_en(direction_2)
        
        e=e1*slip_1+e2*slip_2
        n=n1*slip_1+n2*slip_2
    
        SLIP_DIR[i,2]=e
        SLIP_DIR[i,3]=n
        
    SLIP_DIR[:,7]=list(range(MODEL.shape[0]))
        
    np.savetxt(fslip_dir, SLIP_DIR, fmt="%10.5lf %10.5lf %10.3lf %10.3lf %10.3lf %10.3lf %10.3lf %0d", header='# Slip direction projected onto the surface')


    # save geometry

    print("-- Writing geometry file ",f_geometry)
    
    names=['rdis_long','rdis_lat','rdis_depth','rdis_length','rdis_width',\
           'rdis_area','ratio_rdis_tdis','strike','dip',\
           'centroid_long','centroid_lat','centroid_depth',\
           'tdis_long1','tdis_lat1','tdis_depth1',\
           'tdis_long2','tdis_lat2','tdis_depth2',\
           'tdis_long3','tdis_lat3','tdis_depth3','tdis_area']

    Hdtype={'names':['    rdis_long','   rdis_lat','rdis_depth',\
                    'rdis_length','rdis_width',' rdis_area','ratio_rdis_tdis',\
                    '   strike','       dip',\
                    'centroid_long','centroid_lat','centroid_depth',\
                    'tdis_long1', ' tdis_lat1','tdis_depth1', \
                    'tdis_long2',' tdis_lat2','tdis_depth2', \
                    'tdis_long3',' tdis_lat3','tdis_depth3', \
                    ' tdis_area'],\
             'formats':['f','f','f','f','f','f','f','f','f','f','f','f','f','f','f','f','f','f','f','f','f','f']}
   
    
    SGEOMETRY=lib_inversion.numpy_array_2_numpy_recarray(GEOMETRY,names)
    
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
    
    np.savetxt(f_geometry, GEOMETRY_TXT, fmt=format_header, header=header_cols)

    # save obs

    save_np_array_with_string_2_psvelo_format(OBS, NAME_OBS, '%10.5lf  %10.5lf %10.5lf  %10.5lf %3.1lf %3.1lf %3.1lf %s', f_obs)

    # save obs up
    if ngps_up>0:
        
        MOBSUP = np.zeros( (OBS_UP.shape[0],7) )
        MOBSUP[:,0] = OBS_UP[:,0]
        MOBSUP[:,1] = OBS_UP[:,1]
        MOBSUP[:,3] = OBS_UP[:,2]
        MOBSUP[:,5] = OBS_UP[:,3]
        
        save_np_array_with_string_2_psvelo_format(MOBSUP, NAME_OBS_UP, '%10.5lf  %10.5lf %10.5lf  %10.5lf %3.1lf %3.1lf %3.1lf %s', f_obsup)
