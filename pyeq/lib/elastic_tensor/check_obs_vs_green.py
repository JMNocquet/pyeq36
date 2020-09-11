def check_obs_vs_green( model ):
    
    ###########################################################################
    # MAKING OBSERVATION AND GREEN TENSOR CONSISTENT
    ###########################################################################
    
    # import
    import numpy as np
    import textwrap
    
    print("-- Reordering and extracting Observation and Green tensors")

    # find common gps sites
    np_gps_site , np_common_index_obs_tensor, np_common_index_elastic_tensor = \
    np.intersect1d( model.np_names_t_obs, model.name_obs, assume_unique=True, return_indices=True)

    if model.verbose:
        print('-- List of GPS site both in the Green and Observation Tensor')
        print("%s" % ( textwrap.fill(" ".join( np_gps_site.tolist() ) , width=80 ) ))

    # remove excluded sites
    np_gps_exclude, np_index_exclude, _ = np.intersect1d( np_gps_site, np.array( model.lexclude_gps ), assume_unique=True, return_indices=True )

    if model.verbose:
        print('-- User defined GPS site to be excluded from inversion')
        print("%s" % ( textwrap.fill(" ".join( np_gps_exclude.tolist() ) , width=80 ) ))
    
    np_gps_site = np.delete( np_gps_site, np_index_exclude )
    np_common_index_obs_tensor = np.delete( np_common_index_obs_tensor, np_index_exclude )
    np_common_index_elastic_tensor = np.delete( np_common_index_elastic_tensor, np_index_exclude )

    # making extraction    
    model.t_obs_raw = model.t_obs_raw[ : , np_common_index_obs_tensor, :]
    model.green     = model.green[np_common_index_elastic_tensor]
   
    model.np_gps_site = np_gps_site

    return model 