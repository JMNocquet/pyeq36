def print_observed_time_series( model ):
    """
    Print the time series used in the inversion
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

    ###########################################################################
    # OBS TIME SERIES REALLY USED IN THE INVERSION
    ###########################################################################

    Ogts = Sgts( read=False )

    
    for i in np.arange( model.t_obs.shape[1] ):
        site = model.np_gps_site[i]
        if model.verbose:
            print("-- writing observation time series for site: %s" % ( site ) )
        TS_with_NAN = model.t_obs[:,i,:]
        # remove lines with Nan
        lindex = np.argwhere(np.isnan(TS_with_NAN[:,0]))
        TS = np.delete(TS_with_NAN, lindex, axis=0)        
        
        # fill the gts
        gts = Gts( code = site )
        site_number=np.argwhere( model.name_obs == site )[0,0]
        lon = model.obs[site_number,0]
        lat = model.obs[site_number,1]
        he  = 0.0
        X0,Y0,Z0 = pyacs.lib.coordinates.geo2xyz(lon, lat, he, unit='dec_deg')
        gts.X0 = X0 
        gts.Y0 = Y0 
        gts.Z0 = Z0 

        gts.lon = lon
        gts.lat = lat
        gts.h  = he
        
        gts.data = np.zeros( (TS.shape[0] , 10 ) )
        gts.data[:,0] = at.datetime2decyear( np.delete( model.np_obs_datetime, lindex) ) 
        
        gts.data[ :,1]  = TS[:,1] *1.E-3
        gts.data[ :,2]  = TS[:,0] *1.E-3
        gts.data[ :,3]  = TS[:,2] *1.E-3
        
        gts.data[ :,4] = TS[:,4] *1.E-3
        gts.data[ :,5] = TS[:,3] *1.E-3
        gts.data[ :,6] = TS[:,5] *1.E-3
        
        gts.write_pos(idir=model.odir+'/time_series/obs', force='data' , add_key='obs' , verbose=False)

        Ogts.append( gts )

    model.Ogts = Ogts
    
    return model


    