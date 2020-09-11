def print_displacement_fields( model ):
    """
    Print displacement fields
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
    # MODEL/OBS/RESIDUAL DISPLACEMENT FIELDS FOR EACH OBSERVATION DATE
    ###########################################################################

    format_psvelo="%10.5lf %10.5lf %10.3lf %10.3lf   %10.3lf %10.3lf %3.1lf %s"

    print("-- printing modeled / observed / residuals displacements" )
    
    # LOOP ON OBS DATES
    for i in np.arange( model.np_obs_date_s.shape[0] ):

        DISP_with_NAN = np.copy( model.t_obs[i,:,:] )
        lindex = np.argwhere(np.isnan(DISP_with_NAN[:,0]))
        np_gps_site_i = np.delete( model.np_gps_site , lindex )
        COOR_i = np.delete( model.coor , lindex , axis=0  )

        # MODEL
        if model.verbose:
            print("-- printing modeled displacement at obs date #%04d = %s" % ( i , model.np_obs_date_isoformat[i] ) )
        DISP = np.copy( model.t_mod[i,:,:] )
        DISP = np.delete( DISP , lindex ,axis=0 )
        
        # Horizontal
        COOR_DISP_MODEL = np.zeros( ( DISP.shape[0] , 7 ) )
        COOR_DISP_MODEL[ : , :2 ] = COOR_i
        COOR_DISP_MODEL[ : , 2 ] = DISP[:, 0]
        COOR_DISP_MODEL[ : , 3 ] = DISP[:, 1]
        COOR_DISP_MODEL[ : , 4: ] = np.array([0. , 0., 0.])

        fname = ("%s/displacement/cumulative/model/%04d_model_cum_disp.dat" % (model.odir, i))
        comment = (" model cumulative displacement at obs date #%04d %s" % ( i, model.np_obs_date_isoformat[i] ))
        
        
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

        fname = ("%s/displacement/cumulative/model/%04d_model_cum_disp_up.dat" % (model.odir, i))
        comment = (" model up cumulative displacement at obs date #%04d %s" % ( i, model.np_obs_date_isoformat[i] ))
        
        
        pyacs.lib.utils.save_np_array_with_string( COOR_DISP_MODEL_UP , \
                                                   np_gps_site_i ,\
                                                  format_psvelo, \
                                                  fname, \
                                                  comment )


        # OBSERVATION
        if model.verbose:
            print("-- printing observation displacement at obs date #%04d = %s" % ( i , model.np_obs_date_isoformat[i] ) )
        

        # remove sites with Nan
        DISP = np.delete(DISP_with_NAN, lindex, axis=0)      
        
        # Horizontal
        COOR_DISP_OBS = np.zeros( ( DISP.shape[0] , 7 ) )
        COOR_DISP_OBS[ : , :2 ] = COOR_i
        COOR_DISP_OBS[ : , 2 ] = DISP[:, 0]
        COOR_DISP_OBS[ : , 3 ] = DISP[:, 1]
        COOR_DISP_OBS[ : , 4 ] = DISP[:, 3]
        COOR_DISP_OBS[ : , 5 ] = DISP[:, 4]

        fname = ("%s/displacement/cumulative/obs/%04d_obs_cum_disp.dat" % (model.odir, i))
        comment = (" observed cumulative displacement at obs date #%04d %s" % ( i, model.np_obs_date_isoformat[i] ))
        
        
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

        fname = ("%s/displacement/cumulative/obs/%04d_obs_cum_disp_up.dat" % (model.odir, i))
        comment = (" observed up cumulative displacement at obs date #%04d %s" % ( i, model.np_obs_date_isoformat[i] ))
        
        
        pyacs.lib.utils.save_np_array_with_string( COOR_DISP_OBS_UP , \
                                                   np_gps_site_i ,\
                                                  format_psvelo, \
                                                  fname, \
                                                  comment )

        # RESIDUALS
        if model.verbose:
            print("-- printing residual displacement at obs date #%04d = %s" % ( i , model.np_obs_date_isoformat[i] ) )

        # Horizontal        
        COOR_DISP_RES = COOR_DISP_OBS - COOR_DISP_MODEL
        COOR_DISP_RES[:,:2] = COOR_DISP_OBS[:,:2]
    
        fname = ("%s/displacement/cumulative/res/%04d_res_cum_disp.dat" % (model.odir, i))
        comment = (" observed residual cumulative displacement at obs date #%04d %s" % ( i, model.np_obs_date_isoformat[i] ))
        
        pyacs.lib.utils.save_np_array_with_string( COOR_DISP_RES , \
                                                   np_gps_site_i ,\
                                                  format_psvelo, \
                                                  fname, \
                                                  comment )

        # Up        
        COOR_DISP_RES_UP = COOR_DISP_OBS_UP - COOR_DISP_MODEL_UP
        COOR_DISP_RES_UP[:,:2] = COOR_DISP_OBS_UP[:,:2]
    
        fname = ("%s/displacement/cumulative/res/%04d_res_cum_disp_up.dat" % (model.odir, i))
        comment = (" observed up residual cumulative displacement at obs date #%04d %s" % ( i, model.np_obs_date_isoformat[i] ))
        
        pyacs.lib.utils.save_np_array_with_string( COOR_DISP_RES_UP , \
                                                   np_gps_site_i ,\
                                                  format_psvelo, \
                                                  fname, \
                                                  comment )
