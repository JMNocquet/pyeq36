def print_modeled_time_series( model, save=True ):
    """
    Computes model predicted time series
    modeled time series are stored as pyacs Sgts object and pyacs time series tensors
    optionally, model time series are printing as pos files if write_pos is True
    """

    ###########################################################################
    # IMPORT
    ###########################################################################

    import numpy as np

    from pyacs.lib.gmtpoint import GMT_Point
    import pyacs.lib.utils
    import pyacs.lib.glinalg
    from pyacs.gts.Sgts import Sgts
    from pyacs.gts.Gts import Gts
    import pyacs.lib.coordinates
    import pyacs.lib.astrotime as at

    import pyeq.message.message as MESSAGE
    import pyeq.message.verbose_message as VERBOSE
    import pyeq.message.error as ERROR
    import pyeq.message.warning as WARNING
    import pyeq.message.debug_message as DEBUG

    ###########################################################################
    # PRINT MODEL PREDICTED TIME SERIES
    ###########################################################################
    
    Mgts = Sgts( read=False )
    # save site coordinates for later use for printing displacement files
    
    COOR = np.zeros( (  model.np_gps_site.shape[0] , 2 ) )
    
    # MODEL PREDICTED DISPLACEMENT TIME SERIES
    GREEN = model.green
    OBS = model.obs
    
    # case main rake only
    if model.CUMULATIVE_SLIP_PER_TIME_STEP.shape[2] ==1:
        
        # CS is nfaults, np_model_date_s.shape[0]
        CS = np.squeeze( model.CUMULATIVE_SLIP_PER_TIME_STEP ).T

        # create the interpolated cumulative slip tensor at the observation dates ICS
        ICS = np.zeros( ( model.nfaults, model.np_obs_date_s.shape[0] ) )
        for i in np.arange( model.nfaults ):
            ICS[i,:] = np.interp( model.np_obs_date_s, model.np_model_date_s, CS[i,:] )
            
        # re-order GREEN in component, site, faults
        GREEN_REORDERED = np.swapaxes( model.green[:,:,:,0] , 0 , 1 ).T

        # calculates prediction
        TENSOR_MODEL_TS = np.dot( GREEN_REORDERED , ICS  ).T
        
        # TENSOR_MODEL_TS IS ORDERED: component, site_index, date
        # print results
        for i in np.arange( TENSOR_MODEL_TS.shape[1] ):
            
            site = model.np_gps_site[i]
            DEBUG("printing modeled time series for GPS site: %s" % ( site ) )
            
            TS = TENSOR_MODEL_TS[:,i,:]

            gts = Gts( code = site )
            site_number=np.argwhere( model.name_obs == site )[0,0]
            lon = model.obs[site_number,0]
            lat = model.obs[site_number,1]
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
            gts.data[:,0] = at.datetime2decyear( model.np_obs_datetime ) 
            gts.data[ :,1]  = ( TS[:,1] + model.estimated_offset_per_site[ site ][1] ) *1.E-3
            gts.data[ :,2]  = ( TS[:,0] + model.estimated_offset_per_site[ site ][0] ) *1.E-3
            gts.data[ :,3]  = ( TS[:,2] + model.estimated_offset_per_site[ site ][2] ) *1.E-3
            gts.data[ :,4:] = 1.E-3

            if save:
                gts.write_pos(idir=model.odir+'/time_series/model_all_dates', force='data' , add_key='model_all_dates' , verbose=False)

            # observation dates available for the current site
            TS_with_NAN = model.t_obs[:,i,:]

            # remove lines with Nan
            lindex = np.argwhere(np.isnan(TS_with_NAN[:,0]))
            TS = np.delete( TENSOR_MODEL_TS[:,i,:], lindex, axis=0)        

            # fills and write gts
            gts.data = np.zeros( (TS.shape[0] , 10 ) )
            gts.data[:,0] = np.delete( at.datetime2decyear( model.np_obs_datetime ) , lindex ) 
            gts.data[ :,1]  = TS[:,1] *1.E-3
            gts.data[ :,2]  = TS[:,0] *1.E-3
            gts.data[ :,3]  = TS[:,2] *1.E-3
            gts.data[ :,4:] = 1.E-3

            # set the time series with 0 as the first obs and adds the origin time offsets

            gts.data[ :,1] = gts.data[ :,1 ] - gts.data[ 0,1 ] + model.estimated_offset_per_site[ site ][1] *1.E-3
            gts.data[ :,2] = gts.data[ :,2 ] - gts.data[ 0,2 ] + model.estimated_offset_per_site[ site ][0] *1.E-3
            gts.data[ :,3] = gts.data[ :,3 ] - gts.data[ 0,3 ] + model.estimated_offset_per_site[ site ][2] *1.E-3

            if save:
                gts.write_pos(idir=model.odir+'/time_series/model', force='data' , add_key='model' , verbose=False)
            
            # save for later use
            Mgts.append( gts )
            
    model.Mgts = Mgts
    model.coor = COOR
    model.t_mod = TENSOR_MODEL_TS

    return model
