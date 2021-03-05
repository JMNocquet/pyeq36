"""
Handle the date argument from pyeq_kinematic_inversion.py
"""

def get_np_dates_from_arg_old(arg, ts, verbose):
    
    ###################################################################
    # DEALING WITH DATES
    ###################################################################
    
    """
    returns array of dates from arg
    
    :param arg: string controlling date generation
    :param ts: an Sgts instance of gps time series (gts instance)
    :param verbose: verbose mode
    
    :return np_dates: 1D numpy array of dates for subsequent modeling
    """
    
    import os
    import numpy as np
    import pyacs
    from pyacs.lib import astrotime as at
    
    if verbose:
        print("-- Dealing with the dates")
    
    # tolerance for dates
    tolerance=0.1/365.25
    
    if os.path.isfile( arg ):
        # dates are specified through a file
        
        np_dates=np.genfromtxt( arg ,dtype=str)
        l_available_date=list(map(at.guess_date,np_dates))
        #np_step_date=map(int,np_dates[:,1])
    
    else:
        # dates are specified through an argument string
    
            # Case all dates used
            if  arg =='all':
                print("-- date option:all. Using all dates available")
                l_available_date=[]
                for gts in ts.lGts():
                    l_available_date=l_available_date+gts.data[:,0].tolist()
                l_available_date=np.array(sorted(list(set(sorted(l_available_date)))))
            
            # Case all dates between two dates
            if 'all' in  arg  and '[' in  arg  and ']' in  arg :
                print("-- date option:all observations within a period.")
                tmp_str= arg .replace('[','')
                tmp_str=tmp_str.replace(']','')
                ldate=tmp_str.split(',')
                sdate = at.guess_date(ldate[0])
                edate = at.guess_date(ldate[-1])
                
                l_available_date=[]
                for gts in ts.lGts():
                    l_available_date=l_available_date+gts.data[:,0].tolist()
                l_available_date=np.array(sorted(list(set(sorted(l_available_date)))))
                lindex=np.where((l_available_date>sdate-tolerance) & (l_available_date<edate+tolerance))
                l_available_date=np.array( sorted(set(sorted(  l_available_date[lindex].tolist() ))))
            
            
            if 'all' not in  arg :
                tmp_str= arg .replace('[','')
                tmp_str=tmp_str.replace(']','')
                ldate=tmp_str.split(',')
                sdate = at.guess_date(ldate[0])
                edate = at.guess_date(ldate[-1])
            
            # Case time step provided in days (d option)
                
                tol_day = 0.1 / 365.25
                
                if 'd' in ldate[1]:
                    print("    -- time step provided in days ",ldate[1])
                    mjd_step= float(ldate[1].replace('d',''))
                    
                    l_available_date=at.mjd2decyear( np.arange( \
                        at.decyear2mjd(at.guess_date(sdate)) , at.decyear2mjd(at.guess_date(edate)) , mjd_step) )
    
                    # np.arange does include the last date, so append forces the last date to be included
                    if l_available_date[-1]< edate and (edate - l_available_date[-1]) > 0.1:
                        l_available_date = np.append(l_available_date, edate + tol_day)
                    
                    l_available_date=l_available_date.tolist()
    
                    print("    -- re-arranging dates")
                    l_date_mid_day=[]
                    for date in sorted(l_available_date):
                        date_mid_day_decyear=at.mjd2decyear(int(at.decyear2mjd(date))+0.5)
                        l_date_mid_day.append(date_mid_day_decyear)
                    l_available_date=sorted(l_date_mid_day)
    
                else:
            
                    if len(ldate)==3 :
                # Case time step within a period 
                        if 'y' in ldate[1]:
                            print("    -- time step provided in years ",ldate[1])
                            decyr_step= float(ldate[1].replace('y',''))
                            l_available_date=np.arange(at.guess_date(sdate),at.guess_date(edate),decyr_step)
                            
                            # np.arange does include the last date, so append forces the last date to be included
                            if l_available_date[-1]<edate:
                                l_available_date = np.append(l_available_date, edate)
                            l_available_date=l_available_date.tolist()
            
                            print("    -- re-arranging dates")
                            l_date_mid_day=[]
                            for date in sorted(l_available_date):
                                date_mid_day_decyear=at.mjd2decyear(int(at.decyear2mjd(date))+0.5)
                                l_date_mid_day.append(date_mid_day_decyear)
                            l_available_date=sorted(l_date_mid_day)
                            
                            
                        if 'y' not in ldate[1] and 'd' not in ldate[1]:
                            print("    -- time step provided in number of time steps ",ldate[1])
                            nstep=int(ldate[1])+1
                            l_available_date=np.linspace(sdate, edate, num=nstep).tolist()
            
                            print("    -- re-arranging dates")
                            l_date_mid_day=[]
                            for date in sorted(l_available_date):
                                date_mid_day_decyear=at.mjd2decyear( at.decyear2mjd(date) )
                                #date_mid_day_decyear=at.mjd2decyear(int(at.decyear2mjd(date))+0.5)
                                l_date_mid_day.append(date_mid_day_decyear)
                            l_available_date=sorted(l_date_mid_day)
            
                        
                    else:
                # case list of dates provided
                        l_available_date=list(map( at.guess_date , ldate ))
    
    print(" -- dates used for inversion:")
    
    # NP_DATES ARRAY
    np_dates=np.array(sorted(l_available_date))
    date_ref=np_dates[0]
    
    for i in np.arange(np_dates.shape[0]):
        date=np_dates[i]
        (mday,month,ut)=at.decyear2cal(date)
        (noday,ut)=at.decyear2dayno(date)    
        print(("%04d %10.5lf +%4.1lf days %02d-%02d-%04d-%.1lf %03d" % (i,date,(date-date_ref)*365.25,mday,month,int(date),ut,noday)))
    
    return(np_dates)
