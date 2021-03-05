def check_obs_vs_green( model ):
    
    ###########################################################################
    # MAKING OBSERVATION AND GREEN TENSOR CONSISTENT
    ###########################################################################
    
    ###########################################################################
    # import
    ###########################################################################
    import numpy as np
    import textwrap
    import pyeq.message.message as MESSAGE
    import pyeq.message.verbose_message as VERBOSE
    import pyeq.message.error as ERROR
    import pyeq.message.warning as WARNING
    import pyeq.message.debug_message as DEBUG

    DEBUG("Reordering and extracting Observation and Green tensors")

    # find common gps sites
    np_gps_site , np_common_index_obs_tensor, np_common_index_elastic_tensor = \
    np.intersect1d( model.np_names_t_obs, model.name_obs, assume_unique=True, return_indices=True)

    DEBUG("List of GPS site both in the Green and Observation Tensor")
    DEBUG("%s" % ( textwrap.fill(" ".join( np_gps_site.tolist() ) , width=80 ) ))

    # remove excluded sites
    np_gps_exclude, np_index_exclude, _ = np.intersect1d( np_gps_site, np.array( model.lexclude_gps ), assume_unique=True, return_indices=True )

    DEBUG("User defined GPS site to be excluded from inversion")
    DEBUG("%s" % ( textwrap.fill(" ".join( np_gps_exclude.tolist() ) , width=80 ) ))

    model.warning +="# GPS sites excluded from the inversion in elastic.tensor.chck_obs_vs_green\n"
    model.warning +=("%s" % ( textwrap.fill(" ".join( np_gps_exclude.tolist() ) , width=80 ) ))
    model.warning += "\n"

    np_gps_site = np.delete( np_gps_site, np_index_exclude )
    np_common_index_obs_tensor = np.delete( np_common_index_obs_tensor, np_index_exclude )
    np_common_index_elastic_tensor = np.delete( np_common_index_elastic_tensor, np_index_exclude )

    # making extraction    
    model.t_obs_raw = model.t_obs_raw[ : , np_common_index_obs_tensor, :]
    model.green     = model.green[np_common_index_elastic_tensor]
   
    model.np_gps_site = np_gps_site

    return model 