###################################################################
# PLOT MODEL SHP
###################################################################
def plot_model_gmt(model,
                   cmap='hot',
                   ncolor=None,
                   contour=None,
                   interpolate='cubic',
                   proj="M15c",
                   log=False,
                   title='',
                   bounds=None,
                   outdir='.',
                   outfile=None,
                   line=None,
                   point=None,
                   trench=None,
                   disp=None,
                   disp_scale=0.1):
    """
    Load a pyaks model produces a pdf map
    :param model: shapefile to be plotted
    :param cmap: gmt color map
    :param ncolor: number of discrete colors, None for continuous colors
    :param contour: if None, then colored polygons will be used, number will produce a contour map, refined according to number (usually 2 is fine).
    :param proj: GMT map projection. If None, no projection is used. use crs="EPSG:3395" for a Mercator projection
    :param log: logarithmic map.
    :param title: title string
    :param bounds: bound values for the colormap
    :param outdir: output directory
    :param outfile: outmap name
    :param line: additional line shapefile to be plotted. Can be a list.
    :param point: additional point shapefile to be plotted. Can be a list.
    :param trench: plot trench
    :param disp: displacement/velocity field
    :param disp_scale: scale for disp
    """

    #############################################################################
    # import
    #############################################################################

    import numpy as np
    import pygmt
    from os import mkdir, path

    import pyeq.message.message as MESSAGE
    import pyeq.message.verbose_message as VERBOSE
    import pyeq.message.error as ERROR
    import pyeq.message.warning as WARNING
    import pyeq.message.debug_message as DEBUG

    #############################################################################
    # cmap
    #############################################################################
    #cmap = plt.cm.__dict__[cmap]

    #############################################################################
    # read model
    #############################################################################
    model = np.genfromtxt(model_name)

    #############################################################################
    # converts slip attribute values to numpy array
    #############################################################################
    z = model[:,2]
    if log:
        # modifies z
        z[np.where(z == 0)] = 0.01

    #############################################################################
    # color map bounds
    #############################################################################
    if bounds is None:
        bounds = [np.min(z), np.max(z)]

    #############################################################################
    # start plot
    #############################################################################

    fig = pygmt.Figure()
    pygmt.config(MAP_FRAME_TYPE="plain")
    pygmt.config(FORMAT_GEO_MAP="ddd.xx")

    #############################################################################
    # discrete/continuous color map
    #############################################################################

    pygmt.makecpt(cmap=cmap, series=bounds, reverse=True)

    #############################################################################
    # model polygon map case
    #############################################################################
    if contour is None:

        if log:
            cs = plot_polygon_collection(ax, slip['geometry'][slip.slip > 0], slip['slip'][slip.slip > 0], cmap=cmap,
                                         norm=colors.LogNorm())
            cbar = plt.colorbar(cs, ticks=[0.01, 0.1, 1, 10, 100, 1000])
            cbar.set_label('slip  mm', rotation=270)
        else:
            cs = plot_polygon_collection(ax, slip['geometry'][slip.slip > 0], slip['slip'][slip.slip > 0], cmap=cmap,
                                         norm=norm)
            cbar = plt.colorbar(cs)
            cbar.set_label('slip  mm', rotation=270)

    else:
        #############################################################################
        # contour map case
        #############################################################################
        # get centroid
        listarray = []
        #    for pp in slip.centroid :
        #        listarray.append([pp.x, pp.y])
        records = slip_as_shp.records()
        for record in records:
            listarray.append([record.clon, record.clat])
        coor = np.array(listarray)
        # meshing with Delaunay triangulation
        tri = Triangulation(coor[:, 0], coor[:, 1])
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
        if interpolate not in ['cubic', 'linear']:
            ERROR("interpolate must be linear or cubic. Using linear.")
            interpolate = 'linear'
        if interpolate == 'linear':
            tri_refi, z_test_refi = refiner.refine_field(z, triinterpolator=LinearTriInterpolator(tri, z),
                                                         subdiv=int(contour))
        if interpolate == 'cubic':
            tri_refi, z_test_refi = refiner.refine_field(z, triinterpolator=CubicTriInterpolator(tri, z),
                                                         subdiv=int(contour))
        if log:
            z_test_refi[np.where(z_test_refi < 0.011)] = 0.01
        else:
            z_test_refi[np.where(z_test_refi < 0)] = 0.0

        #############################################################################
        # make actual plot
        #############################################################################
        if log:
            # log plot
            # Alternatively, you can manually set the levels
            # and the norm:
            #        lev_exp = np.arange(np.floor(np.log10(z.min())), np.ceil(np.log10(bounds[-1])+1) , 0.1 )

            # lev_exp = np.arange(np.floor(np.log10(0.01)), np.ceil(np.log10(bounds[-1])+1) , 0.1 )
            # lev_exp = np.arange(np.floor(np.log10(0.01)), np.ceil( bounds[-1] ) , 0.1 )
            lev_exp = np.arange(bounds[0], bounds[1] * 1.1, (bounds[1] * 1.1 - bounds[0]) / ncolor)

            levs = np.power(10, lev_exp)
            cs = plt.tricontourf(tri_refi, z_test_refi, locator=ticker.LogLocator(), levels=levs, norm=colors.LogNorm(),
                                 cmap=cmap)
            cbar = plt.colorbar(cs, ticks=[0.01, 0.05, 0.1, 0.5, 1, 5, 10, 50, 100, 1000])
            # cbar.set_label('slip  mm', rotation=270)
        else:
            # linear plot
            levels = np.linspace(bounds[0], bounds[-1], ncolor)
            plt.tricontourf(tri_refi, z_test_refi, levels=levels, cmap=cmap, norm=norm)
            # colorbar
            sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
            plt.colorbar(sm)

    # fix map bounds
    xlim = plt.gca().get_xlim()
    ylim = plt.gca().get_ylim()

    #############################################################################
    # additional map element to plot
    #############################################################################
    # polygons
    if shp_poly is not None:
        if not isinstance(shp_poly, list):
            shp_poly = [shp_poly]
        for shp_name in shp_poly:
            pshp = geopandas.read_file(shp_name)
            if crs is not None:
                pshp = slip.to_crs(pshp)
            pshp.boundary.plot(edgecolor='k', linewidth=.5, ax=ax)

    # lines
    if shp_line is not None:
        if not isinstance(shp_line, list):
            shp_line = [shp_line]
        for shp_name in shp_line:
            shpl = geopandas.read_file(shp_name)
            if crs is not None:
                shpl = slip.to_crs(shpl)
            plot_linestring_collection(ax, shpl['geometry'], edgecolor='k', linewidth=.5)

    # points
    if shp_point is not None:
        if not isinstance(shp_point, list):
            shp_point = [shp_point]
        for shp_name in shp_point:
            try:
                shpl = geopandas.read_file(shp_name)
                if crs is not None:
                    shpl = slip.to_crs(shpl)
                plot_point_collection(ax, shpl['geometry'], marker='o', markersize=2, color='r')
            except:
                pass

    # displacements as arrows
    if disp is not None:
        if not isinstance(disp, list):
            disp = [disp]

        awidth = 1.5E-2
        lacolor = ['k', 'r', 'g']
        # load
        # max disp
        for i, d in enumerate(disp):
            D = np.genfromtxt(d, usecols=(0, 1, 2, 3))
            max_disp = np.max(np.sqrt(np.sum(D[:, 2:4] ** 2, axis=1)))
            # if max_disp >0:
            #    D[:,3:5] = D[:,3:5] / max_disp
            q = ax.quiver(D[:, 0], D[:, 1], D[:, 2], D[:, 3], units='xy', scale=disp_scale, width=awidth,
                          color=lacolor[i])

    # force map limits to fit model
    plt.gca().set_xlim(xlim)
    plt.gca().set_ylim(ylim)

    # title
    title = ("%s slip (rate) min/max : %.1f/%.1f mm(/day)" % (title, z.min(), z.max()))
    if disp is not None:
        title = ("%s disp max : %.1f mm" % (title, max_disp))

    plt.title(title, fontsize=7)

    # save
    plot_name = ("%s/%s.pdf" % (outdir, shp.split('/')[-1][:-4]))
    plt.savefig(plot_name)

    # close
    plt.close('all')
###################################################################


