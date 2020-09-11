def print_geometry( model ):
    """
    print geometry
    """

    # import
    import numpy as np
    
    ###########################################################################
    # SAVE GEOMETRY IN INFO DIR
    ###########################################################################

    Hdtype={'names':['    rdis_long','   rdis_lat','rdis_depth',\
                    'rdis_length','rdis_width',' rdis_area','ratio_rdis_tdis',\
                    '   strike','       dip',\
                    'centroid_long','centroid_lat','centroid_depth',\
                    'tdis_long1', ' tdis_lat1','tdis_depth1', \
                    'tdis_long2',' tdis_lat2','tdis_depth2', \
                    'tdis_long3',' tdis_lat3','tdis_depth3', \
                    ' tdis_area'],\
             'formats':['f','f','f','f','f','f','f','f','f','f','f','f','f','f','f','f','f','f','f','f','f','f']}
    
    GEOMETRY = model.geometry
    GEOMETRY_TXT=np.zeros((GEOMETRY.shape[0],GEOMETRY.shape[1]+1))
    GEOMETRY_TXT[:,0]=np.arange(GEOMETRY.shape[0])
    GEOMETRY_TXT[:,1:]=GEOMETRY
    header_cols=' '.join(Hdtype['names'])
    format_header='%04d %10.5lf %10.5lf      %6.2lf \
         %6.2lf     %6.2lf    %6.2lf          %6.2lf\
        %6.2lf %10.2lf \
       %10.5lf   %10.5lf         %6.2lf \
    %10.5lf %10.5lf      %6.2lf \
    %10.5lf %10.5lf      %6.2lf \
    %10.5lf %10.5lf      %6.2lf \
       %6.2lf '
    
    np.savetxt( model.odir+'/info/geometry.dat', GEOMETRY_TXT, fmt=format_header, header=header_cols)
    
    
    
    return model