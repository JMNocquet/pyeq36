"""
Routines to build the linear system
"""

def build1(lGts,np_dates,NAME_OBS,GREEN, UP, s_up, date_tol=0.1, verbose=True,debug=False):
    """
    
    :param lGts: list of Gts
    :param np_dates: numpy 1D array of dates
    :param NAME_OBS: numpy 1D array of site name
    :param GREEN: Green tensor as a 4D numpy array
    :param UP: boolean, True means up component will be included
    :param date_tol: tolerance (unit in fractional day) on dates to search for dates in time series corresponding to the modeling period. Use 0.1 days for daily time series.
    :param verbose: verbose mode
    
    :return ATPA, ATPB: LHS matrix and RHS vector ready for inversion algorithm
    """

    print('-- Building the linear system using routine build_linear_system_kinematic_no_uncertainty')
    
    # debug - verbose
    if debug:
        verbose=True
    
    # import
    
    import numpy as np
    import pyacs.lib.astrotime as at
    import resource
    
#    from numba import autojit    
    
    # initialize the list of components to be included in the linear system
    if UP:
        lcomponents=[0,1,2]
    else:
        lcomponents=[0,1]
    
    # shapes of the linear system to be built
    nfaults=GREEN.shape[1]
    nstep = np_dates.shape[0]-1
    nparameters=nfaults*nstep
    
    print('-- nfaults = ',nfaults)
    print('-- number of time steps ', nstep)
    print('-- nparameters ',nparameters)
    print('-- np_date ', np_dates)
    
    # START DATES SECTION
    
    # converts modeling dates to integer seconds
    np_dates_datetime = at.decyear2datetime(np_dates) - at.decyear2datetime(np_dates[0])
    np_dates_s =  np.array(list(map(int,[x.total_seconds() for x in np_dates_datetime])))

    np_step_duration_s = np.diff(np_dates_s)
    np_step_start_s = np_dates_s[0:-1] 
    np_step_end_s   = np_dates_s[1:] 

    nstep = np_step_duration_s.shape[0]

    if verbose:    
        print('-- converted dates of modeling time converted from decimal years to integer seconds')
        print('-- original decimal years modeling dates')
        print(np_dates)
        print('-- modeling dates in seconds wrt reference modeling date')
        print(np_dates_s)
        print('-- modeling start dates in seconds wrt reference modeling date')
        print(np_step_start_s)
        print('-- modeling end   dates in seconds wrt reference modeling date')
        print(np_step_end_s)
        print('-- modeling duration in seconds for each modeling time step')
        print(np_step_duration_s)
        print('-- number of modeling time steps : ', nstep)

    # shrink Gts to the period encompassed by np_dates
    
    wlgts = []
    l_end_date_ts = []
    l_start_date_ts = []
    for gts in lGts:
        
        if verbose:
            print('-- Extracting dates for site ',gts.code)
        
        new_gts = gts.extract_periods( [np_dates[0]-date_tol/365.25 , np_dates[-1] + date_tol/365.25] ,verbose = False )


        
        # time series has no data, skip
        if new_gts.data is None:
            print('! WARNING: no data found at requested time for site ',gts.code,gts.ifile)
            print('! WARNING: site ',gts.code,' will be skipped in the inversion')
            continue
        # time series has only 1 data, skip
        elif new_gts.data.shape[0] == 1:
            print('! WARNING: only one date found at requested time for site ',gts.code,gts.ifile)
            print('! WARNING: site ',gts.code,' will be skipped in the inversion')
            continue
        else:
            # set dN, dE, dU to 0.0 at the first date
            if debug:
                print('-- set dN, dE, dU to 0.0 at the first date ', new_gts.data[0,0], " = ", at.decyear2cal(new_gts.data[0,0]), " = ", at.decyear2dayno(new_gts.data[0,0]))
            new_gts.data[:,1:4] = new_gts.data[:,1:4] - new_gts.data[0,1:4]
            wlgts.append(new_gts)
            # record first and last date of time series
            l_end_date_ts.append(new_gts.data[-1,0])
            l_start_date_ts.append(new_gts.data[0,0])

    # set last date for modeling as min( np_dates[-1], max(l_end_date_ts) ) 

    print('-- first date in np_dates ' ,    np_dates[0],at.decyear2cal(np_dates[0]))
    print('-- first date in time series ' , np.max(l_start_date_ts), at.decyear2cal(np.min(l_start_date_ts)))
    print('-- last  date in np_dates ' ,    np_dates[-1],at.decyear2cal(np_dates[-1]))
    print('-- last  date in time series ' , np.max(l_end_date_ts), at.decyear2cal(np.max(l_end_date_ts)))
    
    np_dates[-1] = np.max(l_end_date_ts)
    np_dates[0]  = np.min(l_start_date_ts)
    
    print('-- start date for modeling is now set to ', np_dates[0], " = ", at.decyear2cal(np_dates[0]), " = ", at.decyear2dayno(np_dates[0]))
    print('-- last  date for modeling is now set to ', np_dates[-1], " = ", at.decyear2cal(np_dates[-1]), " = ", at.decyear2dayno(np_dates[-1]))

    # END DATES SECTION

    # INITIALIZE
    memory_usage = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024 / 1024
    print("-- memory usage before ATPB allocation: %.1lf " % (memory_usage) )
   
    ATB = np.zeros((nfaults*nstep))

    memory_usage = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024 / 1024
    print("-- memory usage after ATPB allocation: %.1lf " % (memory_usage) )


    L_X = {}
#    L_XT = {}
    L_ALPHA_T_ALPHA = {}


    
    # START LOOP ON TIME SERIES
    
    for gts in wlgts:

        IATB = np.zeros((nfaults*nstep))

        print('-- building observation equations for site ',gts.code)

        # index of site for Green's function reading        
        site_number=np.argwhere(NAME_OBS==gts.code.upper())[0,0] 

        # initiate G, should be a 1D numpy array of length nstep
        # GREEN(i,j,k,l) is the prediction for dislocation j at site i component k for rake l
        # k=0,1,2 = east, north, up
        # l=0,1 : rake_principal & rake_conjugate
        # for the moment, we do not handle rake perpendicular component

        G   = {}
        GTG = {}
        
        # BUILD THE ELEMENT ARRAY 
        X = np.zeros((nfaults,nfaults))
        XT = np.zeros((nfaults,1))
        for i_component in lcomponents:
        
            G[i_component] = GREEN[site_number,:,i_component,0].reshape(1,-1)

            # the square matrix having the Green's functions for the given site/component/all faults
            GTG[i_component] = np.dot(G[i_component].T,G[i_component])
            X  += GTG[i_component]
            XT += G[i_component].T
            
        G[0],G[1] = G[1],G[0] 
        L_X[gts.code] = X
        
        # handle dates: all time series dates are converted to integer seconds since the modeling reference date
    
        np_t_datetime = at.decyear2datetime( gts.data[:,0] ) - at.decyear2datetime( np_dates[0] )
        np_t_s =  np.array(list(map(int,[x.total_seconds() for x in np_t_datetime])))

        if verbose:    
            if debug:
                print('-- Time series dates in decimal years converted to integer seconds elapsed since modeling reference date')
                print(gts.data[:,0])
                print(np_dates_s)
            print('-- start time for model ',np_step_start_s[0])
            print('-- start time for observation ',np_t_s[0])
            print('-- end   time for observation ',np_t_s[-1])

    
        # alpha is the time step vector
        
        alpha_tref = np.zeros(nstep)

        # alpha_tref at start observation time
        # this is required because some time series may start later than the reference modeling date
        
        date_s = np_t_s[0]
        lindex = np.where( ( (np_step_end_s >= date_s ) & (np_step_start_s < date_s ) ) )

        if lindex[0].size > 0:
            index_ref = lindex[0][0]
        else:
            index_ref = 0 
    
        alpha_tref[index_ref] = ( float(date_s) - np_step_start_s[index_ref]) / np_step_duration_s[index_ref]
        if index_ref >0:
            alpha_tref[:index_ref] = 1.
        
        # LOOP ON DATES OF THE TIME SERIES
        
        ALPHA_T_ALPHA = np.zeros((nstep,nstep))
        
        for i in np.arange(np_t_s.shape[0])[1:]:

            ALPHA_T = np.zeros((nstep,1))
            
            date_s = np_t_s[i]
#            print '-- ',gts.code,' observation at ',gts.data[i,0],' = ',date_s,' s since the reference modeling date ',np_dates[0]
            if debug:
                print(("-- %s  observation at %15.10lf  = %010d s since the reference modeling date %15.10lf" % (gts.code,gts.data[i,0],date_s,np_dates[0])))
    
            # alpha for the current observation date
            alpha      = np.zeros(nstep)

            index = np.where( ( (np_step_end_s >= date_s ) & (np_step_start_s < date_s ) ) )[0][0]
            alpha[index] = ( float(date_s) - np_step_start_s[index]) / np_step_duration_s[index]
            if index >0:
                alpha[:index] = 1.

            # correct alpha to account for possible late start of observation
            alpha = alpha - alpha_tref

            ALPHA_T_ALPHA += np.dot(alpha.reshape(1,-1).T,alpha.reshape(1,-1))
            ALPHA_T       += alpha.reshape(1,-1).T
            
            # deal with ATPB
            
            for i_component in lcomponents:
                if i_component == 2:
                    IATB += (ALPHA_T.reshape(-1,1) * G[i_component]).flatten() * gts.data[i,i_component+1]*1000.0 * 1./s_up**2
                else:
                    IATB += (ALPHA_T.reshape(-1,1) * G[i_component]).flatten() * gts.data[i,i_component+1]*1000.0
            
                if np.isnan(IATB).any():
                    import sys
                    print('-- IATB has Nan. Exiting.')
                    sys.exit()
            
            # end loop dates of the time series
        
        L_ALPHA_T_ALPHA[gts.code] = ALPHA_T_ALPHA
        
        # fills ATPB
        
        print('-- Filling ATPB RHS vector for site ',gts.code)

        ATB += IATB
        
        # end processing time series

    # end loop dates in time series

    

        # COMBINE X and ALPHA_T_ALPHA

#        @autojit        
    def fill_from_X_1(L_X,L_ALPHA_T_ALPHA, nstep,nfaults):

        # init
        memory_usage = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024 / 1024
        print("-- memory usage before ATPA allocation: %.1lf " % (memory_usage) )
        IATA = np.zeros((nstep*nfaults,nstep*nfaults))
        memory_usage = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024 / 1024
        print("-- memory usage after ATPA allocation: %.1lf " % (memory_usage) )
        
        # loop on time steps
        for l in np.arange(nstep):

            print(("-- filling time step # %04d  over %04d time steps" %(l+1, nstep) ))

        # diagonal blocks
            row_sindex = l * nfaults
            row_eindex = row_sindex + nfaults
            col_sindex = l * nfaults
            col_eindex = col_sindex + nfaults
            
            XX = np.zeros((nfaults,nfaults))
            for code, X in L_X.items():
                XX = XX + X * L_ALPHA_T_ALPHA[code][l,l]
            IATA[row_sindex:row_eindex, col_sindex:col_eindex] = XX

        # off diagonal blocks
            for m in np.arange(l+1,nstep):
                XX = np.zeros((nfaults,nfaults))
                for code, X in L_X.items():
                    XX = XX+ X * L_ALPHA_T_ALPHA[code][l,m] 
                
                row_sindex = l * nfaults
                row_eindex = row_sindex + nfaults
                col_sindex = m * nfaults
                col_eindex = col_sindex + nfaults
                
                IATA[row_sindex:row_eindex, col_sindex:col_eindex] = XX
                IATA[col_sindex:col_eindex, row_sindex:row_eindex] = XX
                
        return(IATA)

    
    print('-- Filling ATA with nstep =',nstep, ' x nfaults =', nfaults,' == ', nstep*nfaults, 'parameters')
    ATA = fill_from_X_1(L_X,L_ALPHA_T_ALPHA, nstep,nfaults)

    # end loop time series

    return(ATA,ATB)

#######################################################################################
def G_d(NEW_GREEN,NAME_OBS,UP,ts,np_dates, s_up, tol_date=5, window=0.1, verbose=False):
#######################################################################################
    """
    Builds the model matrix G and observation vector d
    """
    
    import numpy as np
    import copy
    from pyacs.lib import astrotime as AT 
    
    # init
    
    nfaults=NEW_GREEN.shape[1]
    nparameters=nfaults*(np_dates.shape[0]-1)

    G=np.empty((0,nparameters),dtype=np.float32)
    d=np.empty((0),dtype=np.float32)
    diag_Cd=np.empty((0),dtype=np.float32)

    # dictionary that records the index 
    H_site_index_in_G={}
    H_site_index_ldate={}

    j = 0

    warning = ''

    # loop on time series    
    for gts in ts.lGts():
        print("-- Creating G and d for site ",gts.code)
        
        H_site_index_ldate[gts.code]=[]
    
        # find its number in the Green tensor
        if UP:
            H_site_index_in_G[gts.code]=np.empty((0,3),dtype=int)
        else:
            H_site_index_in_G[gts.code]=np.empty((0,2),dtype=int)
            
        site_number=np.argwhere(NAME_OBS==gts.code.upper())[0,0] 
        
        # find the reference date for this time series
        
        date_ref_site=None
        np_dates_site=copy.deepcopy(np_dates)
        
        for i in np.arange(np_dates.shape[0]):
            # first test the date is possible
            test_date=np_dates[i]
            (mday,month,ut)=AT.decyear2cal(test_date)
            (noday,ut)=AT.decyear2dayno(test_date) 
            if verbose:   
                print(("    -- testing date for reference   : #date in ts %04d %10.5lf +%4.1lf days %02d-%02d-%04d-%.1lf %03d" % (i,test_date,(test_date-np_dates[0])*365.25,mday,month,int(test_date),ut,noday)))
            
            if gts.extract_periods([[test_date-window/365.25,test_date+window/365.25]]).data is None:
                print("!!! No data available at initial date ",test_date," window=",window," days, for site ",gts.code)
                np_dates_site=np_dates_site[1:]
            else:
                date_ref_site=test_date
                H_site_index_ldate[gts.code].append(date_ref_site)
    
               
            if date_ref_site != None:break
        
        if date_ref_site == None:
            print("! WARNING: No usable dates found for time series ",gts.code)
            warning+=("! WARINNG No usable dates found for time series: %s \n" % gts.code)
        
        else:
            (mday,month,ut)=AT.decyear2cal(date_ref_site)
            (noday,ut)=AT.decyear2dayno(date_ref_site)
            if verbose:    
                print(("    -- Reference date for site %s : #date in ts %04d date %10.5lf +%4.1lf days %02d-%02d-%04d-%.1lf %03d" % (gts.code,i,test_date,(test_date-np_dates[0])*365.25,mday,month,int(test_date),ut,noday)))
            if i==0:
                warning+=("--  Reference date for site %s : step# %04d date %10.5lf %02d-%02d-%04d-%.1lf %03d +%4.1lf days since reference date \n" % (gts.code,i,test_date,mday,month,int(test_date),ut,noday,(test_date-np_dates[0])*365.25))
            else:
                warning+=("! WARNING: Reference date for site %s : step# %04d date %10.5lf %02d-%02d-%04d-%.1lf %03d +%4.1lf days since reference date \n" % (gts.code,i,test_date,mday,month,int(test_date),ut,noday,(test_date-np_dates[0])*365.25))
            start_delta_m=i
            np_dates_site=np_dates_site[1:]
            # the 3 rows corresponding to the observation 
            if UP:
                g=np.zeros((3,nparameters))
            else:
                g=np.zeros((2,nparameters))
    
            # BUILDS THE LINEAR SYSTEM
            for i in np.arange(np_dates_site.shape[0]):
                date=np_dates_site[i]
                # step count
                # format dates
                (mday,month,ut)=AT.decyear2cal(date)
                (noday,ut)=AT.decyear2dayno(date)
                if verbose:    
                    print(("    -- testing date for site %s : inversion step %04d date %10.5lf +%4.1lf days %02d-%02d-%04d-%.1lf %03d" % (gts.code,start_delta_m,date,(date-np_dates[0])*365.25,mday,month,int(date),ut,noday)))
                # since observation is N,E,U, g should also also N,E,U g[0,] corresponds to north
                g[0,nfaults*start_delta_m:nfaults*(start_delta_m+1)]=NEW_GREEN[site_number,:,1,0]
                g[1,nfaults*start_delta_m:nfaults*(start_delta_m+1)]=NEW_GREEN[site_number,:,0,0]
                if UP:g[2,nfaults*start_delta_m:nfaults*(start_delta_m+1)]=NEW_GREEN[site_number,:,2,0]
                
                displacement=gts.displacement(sdate=date_ref_site, edate=date, window=window, method='median',verbose=False)
    
                if isinstance(displacement,np.ndarray):
                        if verbose:
                            print("    -- Adding observation for site ",gts.code," at date ",date)
                        G = np.append(G, g, axis=0)
                        if UP:
                            d=np.append(d,displacement[0:3].reshape(3,1)*1.0E3)
                            # rescale up uncertainty
                            displacement[-1]=displacement[-1] * s_up
                            diag_Cd=np.append(diag_Cd,displacement[3:].reshape(3,1)*1.0E3)
                        else:
                            d=np.append(d,displacement[0:2].reshape(2,1)*1.0E3)
                            diag_Cd=np.append(diag_Cd,displacement[3:5].reshape(2,1)*1.0E3)
    
    
                        # we need to record the lines corresponding to the site
                        if UP:
                            H_site_index_in_G[gts.code]=np.append(H_site_index_in_G[gts.code],np.array([[j,j+1,j+2]]),axis=0)
                        else:
                            H_site_index_in_G[gts.code]=np.append(H_site_index_in_G[gts.code],np.array([[j,j+1]]),axis=0)
                            
                        if UP:
                            j=j+3
                        else:
                            j=j+2
                        start_delta_m=start_delta_m+1
                        
                        H_site_index_ldate[gts.code].append(date)
    
                else:
                    print("! WARNING: No observation for site ",gts.code," at inversion step #",start_delta_m," date ",date)
                    
                    (mday,month,ut)=AT.decyear2cal(date)
                    (noday,ut)=AT.decyear2dayno(date)    
                    warning+=("! WARNING: No observation for site %s at inversion step #%d date %10.5lf +%4.1lf days %02d-%02d-%04d-%.1lf %03d\n"  % (gts.code,start_delta_m,date,(date-np_dates[0])*365.25,mday,month,int(date),ut,noday))
                    start_delta_m=start_delta_m+1
               
        print("-- Number of time steps for site ",gts.code,":",len(H_site_index_ldate[gts.code]))
    
    for code in sorted(H_site_index_ldate.keys()):
        print("-- Number of time steps for site ",code,":",len(H_site_index_ldate[code]))
        warning+=("-- Number of time steps for site %s : %d over %d requested by user\n" % (code,len(H_site_index_ldate[code]),np_dates.shape[0]-1))
    
    
    print("-- Number of observations ",d.shape[0])
    if np.isnan(diag_Cd).any():
        print("!!! ERROR building Cd from time series (Nan)")
        print(diag_Cd)
        import sys
        sys.exit()
    
    if np.isnan(d).any():
        print("!!! ERROR building d from time series (Nan)")
        print(d)
        sys.exit()
    
    return(G,d,diag_Cd)

# #######################################################################################
# def G_d_static(NEW_GREEN,NAME_OBS,UP,OBS,OBS_UP, np_dates , s_up, verbose=False):
# #######################################################################################
#     """
#     Builds the model matrix G and observation vector d for the static case
#     """
#     
#     import numpy as np
#     import copy
#     from pyacs.lib import astrotime as AT 
#     
#     # init
#     
#     nfaults=NEW_GREEN.shape[1]
#     nparameters=nfaults*(np_dates.shape[0]-1)
# 
#     G=np.empty((0,nparameters),dtype=np.float32)
#     d=np.empty((0),dtype=np.float32)
#     diag_Cd=np.empty((0),dtype=np.float32)
# 
#     # dictionary that records the index 
#     H_site_index_in_G={}
#     H_site_index_ldate={}
# 
#     j = 0
# 
#     warning = ''
# 
#     # loop on sites
#     for i in np.arange(NAME_OBS.shape[0]):
#         
#         [lon,lat,ve,vn,sve,svn,sven] = OBS[i]
#         
#         if UP:
#             [lon,lat,vu,svu] = OBS_UP[i]
#             
#         
#         displacement = np.array(vn,)
#         
#         code = NAME_OBS[i]
#         print "-- Creating G and d for site " , code
#         
#         H_site_index_ldate[code]=[]
#     
#         # find its number in the Green tensor
#         if UP:
#             H_site_index_in_G[code]=np.empty((0,3),dtype=int)
#         else:
#             H_site_index_in_G[code]=np.empty((0,2),dtype=int)
#             
#         site_number=np.argwhere(NAME_OBS==code.upper())[0,0] 
#         
#         
#         # the 3 rows corresponding to the observation 
#         if UP:
#             g=np.zeros((3,nparameters))
#         else:
#             g=np.zeros((2,nparameters))
# 
#         # BUILDS THE LINEAR SYSTEM
# 
#         # since observation is N,E,U, g should also also N,E,U g[0,] corresponds to north
#         g[0,0:nfaults]=NEW_GREEN[site_number,:,1,0]
#         g[1,0:nfaults]=NEW_GREEN[site_number,:,0,0]
# 
#         if UP:
#             g[2,0:nfaults]=NEW_GREEN[site_number,:,2,0]
# 
#         G = np.append(G, g, axis=0)
# 
#         if UP:
#             d=np.append(d,displacement[0:3].reshape(3,1)*1.0E3)
#             # rescale up uncertainty
#             displacement[-1]=displacement[-1] * s_up
#             diag_Cd=np.append(diag_Cd,displacement[3:].reshape(3,1)*1.0E3)
#         else:
#             d=np.append(d,displacement[0:2].reshape(2,1)*1.0E3)
#             diag_Cd=np.append(diag_Cd,displacement[3:5].reshape(2,1)*1.0E3)
# 
# 
#             # we need to record the lines corresponding to the site
#             if UP:
#                 H_site_index_in_G[gts.code]=np.append(H_site_index_in_G[gts.code],np.array([[j,j+1,j+2]]),axis=0)
#             else:
#                 H_site_index_in_G[gts.code]=np.append(H_site_index_in_G[gts.code],np.array([[j,j+1]]),axis=0)
#                 
#             if UP:
#                 j=j+3
#             else:
#                 j=j+2
#             start_delta_m=start_delta_m+1
#             
#             H_site_index_ldate[gts.code].append(date)
# 
#         else:
#             print "! WARNING: No observation for site ",gts.code," at inversion step #",start_delta_m," date ",date
#             
#             (mday,month,ut)=AT.decyear2cal(date)
#             (noday,ut)=AT.decyear2dayno(date)    
#             warning+=("! WARNING: No observation for site %s at inversion step #%d date %10.5lf +%4.1lf days %02d-%02d-%04d-%.1lf %03d\n"  % (gts.code,start_delta_m,date,(date-np_dates[0])*365.25,mday,month,int(date),ut,noday))
#             start_delta_m=start_delta_m+1
#            
#     print "-- Number of time steps for site ",gts.code,":",len(H_site_index_ldate[gts.code])
#     
#     for code in sorted(H_site_index_ldate.keys()):
#         print "-- Number of time steps for site ",code,":",len(H_site_index_ldate[code])
#         warning+=("-- Number of time steps for site %s : %d over %d requested by user\n" % (code,len(H_site_index_ldate[code]),np_dates.shape[0]-1))
#     
#     
#     print "-- Number of observations ",d.shape[0]
#     if np.isnan(diag_Cd).any():
#         print "!!! ERROR building Cd from time series (Nan)"
#         print diag_Cd
#         import sys
#         sys.exit()
#     
#     if np.isnan(d).any():
#         print "!!! ERROR building d from time series (Nan)"
#         print d
#         import sys
#         sys.exit()
#     
#     return(G,d,diag_Cd)

