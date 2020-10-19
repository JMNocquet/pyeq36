"""
Handle the date argument from pyeq_kinematic_inversion.py
"""


def get_np_dates_from_arg(arg, np_obs_date_s, rounding='day', verbose=False):

    """
    returns array of dates in integer seconds since 1980/1/1 00:00:00
    
    :param arg: string controlling date generation
    :param np_obs_date_s : 1D numpy array of all available bservation dates in integer seconds since 1980/1/1 00:00:00 
    :param rounding : controls rounding. 'day' means at 12:00:00, 
    'hour' at 00:00, 'minute' at the current minute with 00 seconds, 'second' at the integer of the current second.
    :param verbose: verbose mode
    
    :return: 1D integer (np.int64) numpy array 
    """

    ###########################################################################    
    def round_np_datetime(np_datetime, rounding ):
    ###########################################################################    

        from datetime import timedelta

        for i in np.arange( np_datetime.shape[0] ):
            if rounding == 'day':
                np_datetime[i] = np_datetime[i].replace(hour= 12, minute = 0, second = 0, microsecond = 0)
            if rounding == 'hour':
                if np_datetime[i].minute >= 30:
                    np_datetime[i] = np_datetime[i] + timedelta( minutes=30 ) 
                np_datetime[i] = np_datetime[i].replace( minute = 0, second = 0, microsecond = 0)
            if rounding == 'minute':
                if np_datetime[i].second >= 30:
                    np_datetime[i] = np_datetime[i] + timedelta( seconds=30 ) 
                np_datetime[i] = np_datetime[i].replace( second = 0, microsecond = 0)
            if rounding == 'second':
                np_datetime[i] = np_datetime[i].replace( microsecond = 0)
        
        return np_datetime
    ###########################################################################    

    
    import os
    import numpy as np
    from pyacs.lib import astrotime as at
    from datetime import datetime
    import pandas

    # reference date is 1980.0
    ref_date_time = datetime(1980, 1, 1, 0, 0, 0)
    
    if verbose:
        print("-- Dealing with the dates using pyeq.lib.dateget_np_dates_from_arg using rounding=%s" , rounding )

    ###################################################################
    # DATES PROVIDED AS PANDAS.DATE_RANGE COMMAND
    # THE COMMAND IS USE TO GENERATE A FILE EXECUTED BY PYTHON
    # THAT WILL CREATE A DATE FILE
    ###################################################################
    
    if 'pandas.date_range' in arg:
        # creates python code file
        with open('user_requested_model_dates.py', 'w') as file:
            file.write("#!/usr/bin/env python \n")
            file.write("import pandas\n")
            file.write("import numpy as np\n")
            file.write("np_datetime = %s.to_pydatetime()\n" % arg )
            file.write("np.save('np_datetime.npy' , np_datetime)")
        
        # execute code
        from pyacs.lib import syscmd
        syscmd.getstatusoutput('python user_requested_model_dates.py')
        
        # load
        np_datetime = np.load( 'np_datetime.npy' , allow_pickle=True )
        # convert to seconds
        np_date_s =  np.array(list(map(int,[x.total_seconds() for x in ( np_datetime - ref_date_time) ])), dtype=np.int64)
        
        # clean files
        os.remove( 'np_datetime.npy' )
        os.remove( 'user_requested_model_dates.py' )
        
        # return
        return np_date_s
    
    ###################################################################
    # DATES PROVIDED AS A NPY FILE
    ###################################################################
    
    if os.path.isfile( arg ) and ( arg[-4:] == '.npy' ) :

        np_datetime = np.load( arg , allow_pickle=True )

        # handle rounding
        np_datetime = round_np_datetime(np_datetime, rounding)
        
        # convert to timedelta and then to seconds
        np_date_s =  np.array(list(map(int,[x.total_seconds() for x in ( np_datetime - ref_date_time) ])), dtype=np.int64)
        
        print("-- dates read from npy file: %s" % arg )
        return np_date_s
    
    ###################################################################
    # DATES PROVIDED AS A TEXT FILE
    # FORMAT IS EXPECTED TO BE
    #2016-04-17 12:00:00
    #2016-04-18 12:00:00
    #2016-04-19 12:00:00
    # ....
    ###################################################################

    if os.path.isfile( arg ) and ( arg[-4:] != '.npy' ) :
    
        # load file
        d = np.genfromtxt( arg ,delimiter=',',dtype=str) 
        # parser
        parse = lambda x: pandas.datetime.strptime(x, '%Y-%m-%d %H:%M:%S') 
        # get np_datetime
        np_datetime = np.array(list(map(parse, d)))
        # convert to seconds
        np_date_s =  np.array(list(map(int,[x.total_seconds() for x in ( np_datetime - ref_date_time) ])), dtype=np.int64)
        
        print("-- dates read from text file: %s" % arg )
        return np_date_s

    ###################################################################
    # DATES PROVIDED AS A ARGUMENT STRING
    ###################################################################

    # dates are specified through an argument string

    # Case all dates used
    
    if  arg =='all':
        print("--- date option:all. Using all available observation dates. Rounding ignored.")
        return np.copy( np_obs_date_s )
    
    # Case all dates between two dates
    if 'all' in  arg  and '[' in  arg  and ']' in  arg :
        print("-- date option: all observations within a period. Rounding ignored.")
        tmp_str= arg .replace('[','')
        tmp_str=tmp_str.replace(']','')
        ldate=tmp_str.split(',')
        sdate = at.guess_date(ldate[0])
        edate = at.guess_date(ldate[-1])
        
        sdate_s = at.decyear2seconds( sdate )
        edate_s = at.decyear2seconds( edate )

        lindex = np.where( (np_obs_date_s >= sdate_s ) &  (np_obs_date_s <= edate_s ) )
        np_date_s = np_obs_date_s[ lindex ]

        return np_date_s
    
    # Case without keyword 'all'         

    if 'all' not in  arg and '[' in  arg  and ']' in  arg :
        tmp_str= arg .replace('[','')
        tmp_str=tmp_str.replace(']','')
        ldate=tmp_str.split(',')
        sdate = at.guess_date(ldate[0])
        edate = at.guess_date(ldate[-1])
    
        sdate_s = at.decyear2seconds( sdate , rounding=rounding )
        edate_s = at.decyear2seconds( edate , rounding=rounding )
    
    # Case time step provided in days (d option)
        
        if 'd' in ldate[1]:
            print("    -- time step provided in days ",ldate[1])
            # define the associated delta in seconds
            delta_s = int( float(ldate[1].replace('d','')) * 86400. )
            np_date_s = np.arange(sdate_s,edate_s+delta_s,step=delta_s) 
            # round again
            np_datetime = at.seconds2datetime( np_date_s )
            np_datetime = round_np_datetime(np_datetime, rounding)

            # convert to timedelta and then to seconds
            np_date_s =  np.array(list(map(int,[x.total_seconds() for x in ( np_datetime - ref_date_time) ])), dtype=np.int64)
            
            return np_date_s

    # Case time step provided in years (y option)

        if 'y' in ldate[1]:
            print("    -- time step provided in years ",ldate[1])
            # define the associated delta in seconds
            delta_s = int( float(ldate[1].replace('y','')) * 86400. * 365.25 )
            np_date_s = np.arange(sdate_s,edate_s+delta_s,step=delta_s) 
            # round again
            np_datetime = at.seconds2datetime( np_date_s )
            np_datetime = round_np_datetime(np_datetime, rounding)

            # convert to timedelta and then to seconds
            np_date_s =  np.array(list(map(int,[x.total_seconds() for x in ( np_datetime - ref_date_time) ])), dtype=np.int64)
            
            return np_date_s

    # Case number of time step provided between two dates

        if ldate[1].isdigit():
            
            np_date_s = np.linspace( sdate_s, edate_s , int( ldate[1] ) ).astype( int )
            # round again
            np_datetime = at.seconds2datetime( np_date_s )
            np_datetime = round_np_datetime(np_datetime, rounding)

            # convert to timedelta and then to seconds
            np_date_s =  np.array(list(map(int,[x.total_seconds() for x in ( np_datetime - ref_date_time) ])), dtype=np.int64)
            
            return np_date_s

    # ERROR
    print("!!!ERROR in pyeq.lib.date_args2.args2")
    print("!!!Could not decipher argument |%s| as a date argument")
    import sys
    sys.exit()
           
            