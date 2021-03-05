def print_residual_time_series( model ):
    """
    Print residual time series
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
    # RESIDUAL TIME SERIES
    ###########################################################################

    Rgts = Sgts( read=False )


    stats_site = np.zeros( ( model.np_gps_site.shape[0] , 9 ) )

    RES = np.zeros( (0, 10 ) )

    for site in model.np_gps_site:
        DEBUG("writing residual time series for site: %s "  % site )

        # copy obs to res
        res_gts = model.Ogts.__dict__[site].copy()

        # remove model prediction
        # caution: model prediction needs to be set to 0 for the first available date
        # this probably needs to be changed if an origin constant has been estimated
#        res_gts.data[:,1:4] = res_gts.data[:,1:4] - ( model.Mgts.__dict__[site].data[:,1:4] - model.Mgts.__dict__[site].data[0,1:4] )
        res_gts.data[:,1:4] = res_gts.data[:,1:4] - model.Mgts.__dict__[site].data[:,1:4]

        # write pos file
        res_gts.write_pos(idir=model.odir+'/time_series/res', force='data' , add_key='res' , verbose=False)
        Rgts.append( res_gts )
        
    model.Rgts = Rgts
    
    return model
