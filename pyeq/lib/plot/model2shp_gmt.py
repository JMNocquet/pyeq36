"""
Converts slip models from pyeq to shapfile and gmt
"""

def model2shp_gmt(geometry, type_dis, lslip, out_dir_shp=None, out_dir_gmt=None , verbose=True ):
    """
    :param geometry: geometry as either .dat, .npy, or 2D numpy array
    :param type_dis: dislocation type either 'tde' or 'rde'
    :param lslip: list of slip files
    :param out_dir_shp: output directory for shapefiles
    :param out_dir_gmt: output directory for gmt files
    :param verbose: verbose mode
    """
    
    ###########################################################################
    # import
    ###########################################################################

    import shapefile
    from pyeq.lib import eq_disloc_3d as DL
    import numpy as np 
    
    ###########################################################################
    # type_dis
    ###########################################################################
    
    TRIANGLE = False
    if type_dis == 'tde':
        TRIANGLE = True
    
    ###########################################################################
    # geometry
    ###########################################################################

    GEOMETRY = None
    
    if isinstance( geometry , str ):
        if geometry[-3:] == 'dat':
            # reads the dat text file
            import sys
            print('!!!TODO: geometry as dat file not implemented yet.')
            sys.exit()
        if geometry[-3:] == 'npy':
            # reads the npy
            import pyeq.lib.geometry.to_np_array
            GEOMETRY = pyeq.lib.geometry.to_np_array.npy_geometry_to_array_and_recarray( geometry , verbose=verbose )
        
    if isinstance( geometry, np.ndarray ):
        GEOMETRY = geometry
    
    if GEOMETRY is None:
        print('!!!ERROR: Could not understand argument: geometry')
        import sys
        sys.exit()
    
    
    ###########################################################################
    # loop on slip files
    ###########################################################################
    for slip_dat in lslip:
    
        SLIP=np.genfromtxt(slip_dat,comments='#')

        lfaults=[]
        lrecord=[]
        for i in np.arange(GEOMETRY.shape[0]):
            ( rdis_long, rdis_lat, rdis_depth, rdis_length, rdis_width,  rdis_area, ratio_rdis_tdis,strike,dip,\
             centroid_long, centroid_lat, centroid_depth,\
             tdis_long1,tdis_lat1,tdis_depth1,tdis_long2,tdis_lat2,tdis_depth2,tdis_long3,tdis_lat3,tdis_depth3,tdis_area)\
             =GEOMETRY[i,:]
            depth=rdis_depth
        
    
            if TRIANGLE:
                lfaults.append([ [tdis_long1,tdis_lat1], [tdis_long2,tdis_lat2], [tdis_long3,tdis_lat3] ])
                lrecord.append([i,centroid_depth,SLIP[i,-1],0])
        
            else:
            
                # creates a dislocation object
                disloc=DL.Dislocation(i,rdis_long,rdis_lat,rdis_depth/one_degree, strike, dip, rdis_length/one_degree, rdis_width/one_degree,rdis_area, 0, 0)
                # get the corners
                (X1,X2,X3,X4)=disloc.corners()
                lfaults.append([ [X1[0],X1[1]], [X2[0],X2[1]], [X3[0],X3[1]], [X4[0],X4[1]], [X1[0],X1[1]] ])
                lrecord.append([i,depth,SLIP[i,-1],0])
    
        slip_file_basename = slip_dat.split('/')[-1].split('.')[-2]

        ###################################################################
        # WRITES GMT PSXY FILES
        # This file can be then plotted with
        #  psxy lima_simple_sigma_010_dc_050_m0_000_sol_coupling.gmt -V -R-81/-72/-16/-7 -JM14 -L -m -W0.2/0  -Clut_coupling.cpt > test.ps
        # triangles have not been tested yet
        ###################################################################
        
        if out_dir_gmt is not None:
        
        #gmtfile=os.path.abspath(slip_dat).split('.')[0]+'.gmt'
            gmtfile = ("%s/%s.gmt" % ( out_dir_gmt , slip_file_basename ) )
            
            if verbose:
                print(("-- saving gmt file %s " % gmtfile ))

            f=open(gmtfile,'w')
            for i in np.arange(len(lfaults)):
                [index,depth,slip,rake]=lrecord[i]
                f.write('> -Z%.3lf\n'%slip)
                fault=lfaults[i]
                for xy in fault:
                    f.write("%10.3lf %10.3lf\n" %(xy[0],xy[1]))
            f.write('>\n')
            f.close()

        ###################################################################
        # WRITES SHAPEFILES
        ###################################################################

        if out_dir_shp is not None:
    
            shp_file = ("%s/%s" % ( out_dir_shp , slip_file_basename ) )
#            shp_file=os.path.abspath(slip_dat).split('.')[0]
     
           ###################################################################
            # INITIALIZE SHAPEFILE
            ###################################################################
            
            # Make a polygon shapefile
            w = shapefile.Writer( shp_file ,shapeType=shapefile.POLYGON)
            w.field('ID','I','40')
            w.field('i_subfault','F','40')
            w.field('depth_top_disloc','F','40')
            w.field('slip','F', decimal=5 )
            w.field('rake','F',decimal=1 )
    
            
            ###################################################################
            # LOOP ON FAULTS
            ###################################################################
            
            for i in np.arange(len(lfaults)):
                fault=lfaults[i]
                record=lrecord[i]
                w.poly([fault])
                [index,depth,slip,rake]=lrecord[i]
                w.record(str(i),index,depth,slip,rake)
                i=i+1
            
            ###################################################################
            # SAVE SHAPEFILE
            ###################################################################
            if verbose:
                print(("-- saving shapefile %s " % (shp_file.split('/')[-1]+'.shp')))
            w.close()
            
            ###################################################################
            # SAVE .PRJ FILE
            ###################################################################
            
            prj = open("%s.prj" % shp_file, "w") 
            epsg = 'GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563]],PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433]]' 
            prj.write(epsg) 
            prj.close()
        