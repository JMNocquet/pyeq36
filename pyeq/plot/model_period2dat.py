def model_period2dat( model , sdate , edate, dat_file, rate=False ):
  """
  Creates a date file corresponding to the cumulated slip model between two dates
  :param model: a pyeq model instance
  :param sdate, edate: start & end date as a datetime.datetime instance
  :dat_file: output shapefile name
  """
  # import
  import numpy as np
  
  # 
  
  # get the index for sdate and edate from model.np_model_datetime
  i = 0
  while model.np_model_datetime[i] <= sdate:
    i = i+1
    pass
  idx_sdate = i-1

  i = 0
  while model.np_model_datetime[i] <= edate:
    i = i+1
    pass
  idx_edate = i-1
  
  # making the cumulated slip vector
  CUMULATIVE_SLIP_TIME_STEP=np.zeros(( model.nfaults, 2 + model.np_rake.shape[0] ))
  CUMULATIVE_SLIP_TIME_STEP[:,0] = model.sgeometry.centroid_long
  CUMULATIVE_SLIP_TIME_STEP[:,1] = model.sgeometry.centroid_lat
  CUMULATIVE_SLIP_TIME_STEP[:,2:] = model.CUMULATIVE_SLIP_PER_TIME_STEP[idx_edate].reshape( model.nfaults, -1 ) -model.CUMULATIVE_SLIP_PER_TIME_STEP[idx_sdate].reshape( model.nfaults, -1 )

  
  
  # save dat file
  fname=("%s.dat" % ( dat_file ))
  
  # average rate
  if rate:
    timedelta = model.np_model_datetime[idx_edate] - model.np_model_datetime[idx_sdate]
    delta_day = timedelta.days + timedelta.seconds / (24*60*60)
    date_info = ("average slip rate for period from %s to %s " % ( model.np_model_date_isoformat[idx_sdate] , model.np_model_date_isoformat[idx_edate] ) )
    CUMULATIVE_SLIP_TIME_STEP[:,2] = CUMULATIVE_SLIP_TIME_STEP[:,2] / delta_day
  # cumulated slip
  else:
    date_info = ("cumulated slip for period from %s to %s " % ( model.np_model_date_isoformat[idx_sdate] , model.np_model_date_isoformat[idx_edate] ) )
  
  np.savetxt(fname, CUMULATIVE_SLIP_TIME_STEP, fmt="%10.5lf %10.5lf %10.3lf", header=date_info)
  
  return
