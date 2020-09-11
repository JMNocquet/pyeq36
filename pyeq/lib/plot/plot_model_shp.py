###################################################################
# PLOT MODEL SHP
###################################################################
def plot_model_shp( shp , \
                    cmap='jet' ,\
                    ncolor=None ,\
                    contour=None,\
                    crs=None,\
                    log=False ,\
                    title='' ,\
                    bounds=None ,\
                    outdir='.',\
                    outfile=None ,\
                    shp_poly=None,\
                    shp_line=None ,\
                    shp_point=None ):
  """
  Load a model polygon shapefile with slip attributes and produces a pdf map
  :param shp: shapefile to be plotted
  :param cmap: matplotlib color map
  :param ncolor: number of discrete colors, None for continuous colors
  :param contour: if None, then colored polygons will be used, number will produce a contour map, refined according to number (usually 2 is fine).
  :param crs: crs for map projection. If None, no projection is used. use crs="EPSG:3395" for a Mercator projection
  :param log: logarithmic map.
  :param title: title string
  :param bounds: bound values for the colormap
  :param outdir: output directory
  :param outfile: outmap name
  :param shp_poly: additional polygon shapefile to be plotted
  :param shp_line: additional line shapefile to be plotted. Can be a list.
  :param shp_point: additional point shapefile to be plotted. Can be a list.
  """
  
  # import 
  import numpy as np
  from os import mkdir, path
  import geopandas
  from geopandas.plotting import plot_polygon_collection, plot_linestring_collection, plot_point_collection    
  import matplotlib as mpl
  import matplotlib.pylab as plt
  import matplotlib.colors as colors
  from matplotlib.tri import Triangulation, TriAnalyzer, UniformTriRefiner, LinearTriInterpolator
  from matplotlib import ticker
  from matplotlib.ticker import LogFormatter        

  # cmap
  cmap=plt.cm.__dict__[cmap]
  
  # read shapefile
  slip = geopandas.read_file( shp )
  
  # converts to projection
  # Mercator
  #slip = slip.to_crs("EPSG:3395")
  if crs is not None:
      slip = slip.to_crs( crs )
  
  # converts slip attribute values to numpy array
  listarray=[]
  for s in slip.slip: 
    listarray.append([s]) 
  z = np.array(listarray).flatten()
  if log:
    # modifies z
    z[np.where(z==0)] = 0.01
  
  # bounds
  if bounds is None:
    bounds = [np.min(z),np.max(z)]
  
  # start plot
  fig, ax = plt.subplots()
  ax.set_aspect('equal')
  
  
  # discrete/continuous color map
  if ncolor is None:
    ncolor=1000

  # bound_colors
  bound_colors = np.linspace(bounds[0], bounds[-1] , ncolor)
  norm = mpl.colors.BoundaryNorm(bound_colors, cmap.N)
  # extract all colors from cmap
  cmaplist = [cmap(i) for i in range(cmap.N)]
  # force the first color entry to be grey
  cmaplist[0] = (1., 1., 1., 1. )
  # create the new map
  cmap = mpl.colors.LinearSegmentedColormap.from_list('my_cmap', cmaplist, cmap.N)

  # plot the grid (polygon and contour map cases)
  plot_polygon_collection(ax, slip['geometry'], edgecolor='k', facecolor='w', linewidth=.1)

  # model polygon map case
  if contour is None:
#    plot_polygon_collection(ax, slip['geometry'], edgecolor='k',facecolor='w',linewidth=.1)
    if log:
        cs = plot_polygon_collection(ax, slip['geometry'][slip.slip>0], slip['slip'][slip.slip>0],cmap=cmap, norm=colors.LogNorm()  )
        cbar = plt.colorbar( cs , ticks=[0.01,0.1,1,10,100,1000] )
        cbar.set_label('slip  mm', rotation=270)
    else:
        cs = plot_polygon_collection(ax, slip['geometry'][slip.slip>0], slip['slip'][slip.slip>0],cmap=cmap, norm=norm )
        cbar = plt.colorbar( cs  )
        cbar.set_label('slip  mm', rotation=270)
  
  else:
  # contour map case
    #plot_polygon_collection(ax, slip['geometry'], edgecolor='k', linewidth=.1)
    # get centroid
    listarray=[]
    for pp in slip.centroid : 
        listarray.append([pp.x, pp.y]) 
    coor = np.array(listarray) 
    # meshing with Delaunay triangulation
    tri = Triangulation( coor[:,0], coor[:,1])
    ntri = tri.triangles.shape[0]
    # masking badly shaped triangles at the border of the triangular mesh.
    min_circle_ratio = .1  # Minimum circle ratio - border triangles with circle
                            # ratio below this will be masked if they touch a
                            # border. Suggested value 0.01; use -1 to keep
                            # all triangles.

    mask = TriAnalyzer(tri).get_flat_tri_mask(min_circle_ratio)
    tri.set_mask(mask)
    # refining the data
    refiner = UniformTriRefiner(tri)
    tri_refi, z_test_refi = refiner.refine_field( z , triinterpolator=LinearTriInterpolator(tri,z), subdiv=int(contour))
    if log:
        z_test_refi[np.where(z_test_refi<0.011)] = 0.01
    else:
        z_test_refi[np.where(z_test_refi<0)] = 0.0

    # make actual plot
    if log:
        # Alternatively, you can manually set the levels
        # and the norm:
#        lev_exp = np.arange(np.floor(np.log10(z.min())), np.ceil(np.log10(bounds[-1])+1) , 0.1 )
        lev_exp = np.arange(np.floor(np.log10(0.01)), np.ceil(np.log10(bounds[-1])+1) , 0.1 )
        levs = np.power(10, lev_exp)
        cs = plt.tricontourf(tri_refi, z_test_refi, locator=ticker.LogLocator() , levels=levs, norm=colors.LogNorm() , cmap=cmap )
        cbar = plt.colorbar( cs , ticks=[0.01,0.1,1,10,100,1000] )
        #cbar.set_label('slip  mm', rotation=270)
    else:
        levels = np.linspace(bounds[0], bounds[-1], ncolor )        
        plt.tricontourf(tri_refi, z_test_refi, levels=levels, cmap=cmap , norm=norm )
        # colorbar
        sm = plt.cm.ScalarMappable(cmap=cmap , norm=norm )
        plt.colorbar(sm)
  
  
  # fix map bounds
  xlim = plt.gca().get_xlim()
  ylim = plt.gca().get_ylim()
  

  # additional map element to plot
  if shp_poly is not None:
    if not isinstance(shp_poly,list):
      shp_poly = [shp_poly]
    for shp_name in shp_poly:
        pshp = geopandas.read_file( shp_name )
        if crs is not None:
            pshp = slip.to_crs( pshp )
        pshp.boundary.plot( edgecolor='k', linewidth=.5, ax=ax)

        #plot_polygon_collection(ax, pshp['geometry'], facecolor=None, edgecolor='k',linewidth=.5)

  if shp_line is not None:
    if not isinstance(shp_line,list):
      shp_line = [shp_line]
    for shp_name in shp_line:  
      print("-- reading " , shp_name )
      shpl = geopandas.read_file( shp_name )
      if crs is not None:
          shpl = slip.to_crs( shpl )
      print('-- now plotting ', shp_line)
      plot_linestring_collection(ax, shpl['geometry'], edgecolor='k',linewidth=.5)

      print('-- OK ', shp_line)



  if shp_point is not None:
    if not isinstance(shp_point,list):
      shp_point = [shp_point]
    for shp_name in shp_point:  
      shpl = geopandas.read_file( shp_name )
      if crs is not None:
          shpl = slip.to_crs( shpl )
      plot_point_collection(ax, shpl['geometry'], marker='o' , markersize=2,color='r')
   
  # force map limits to fit model
  plt.gca().set_xlim( xlim )
  plt.gca().set_ylim( ylim )
    
  # title
  title = ("%s - min/max : %.1f/%.1f" % (title , z.min(),z.max() ))
  plt.title( title , fontsize=8 )
  
  # save
  plot_name = ("%s/%s.pdf" % ( outdir , shp.split('/')[-1][:-4] ) )
  plt.savefig( plot_name )
  
  # close
  plt.close( 'all' )
###################################################################


