Tutorial static inversion
=========================

.. toctree::
   :maxdepth: 1


There are 3 steps for running an inversion:
* setup geodetic data (observations) files
* build the geometry
* calculate the Green's function
* run the inversion to get the models

So, we create four directories:

::

	mkdir data ; mkdir geometry ; mkdir green ; mkdir models

Input data
----------

You will need the following data sets:

* either a gmt format displacement/velocity file
* or a directory of time series in GAMIT/GLOBK PBO format
* a geometry of your fault provided as a gmt grd file

In the following, I'll assume that that you have a file data/mw67_coseismic.dat which looks like:

::

	# file mw67_coseismic.dat
	-78.62765   -1.26861  -0.20  -2.00   1.00   1.00   0.00 ABEC
 	-78.84741   -2.20216   4.44  -4.89   1.00   1.00   0.00 ALEC
 	-78.54539    0.91264   0.00  -3.05   1.00   1.00   0.00 ALTB
 	-78.16230   -0.46342   2.17  -4.00   1.00   1.00   0.00 ANTN
 	-78.17040   -0.49729  -0.40  -2.01   1.00   1.00   0.00 ANTS
 	-79.09820    0.09634  -4.42  -0.77   1.00   1.00   0.00 ARSH
	.
	. 


Building the fault geometry from a grd file
-------------------------------------------

This part assumes GMT 4.x to be properly installed on your system. Building the fault geometry can be done with the script 
::

	pyeq_parametrize_curve_surface_triangles.py

If not already done, you will first need to install the Polygon2 package, which is available like `here <https://pypi.python.org/pypi/Polygon2/>`_. 
::

	pip install Polygon2-2.0.6.tar.gz

You then need to link the grd file to you local directory 
::

	cd geometry
	ln -s /usr/local/geodesy/maps/geophysical_data/sam_slab1.0_clip_shifted.grd .

Finally, you run the script pyeq_parametrize_curve_surface_triangles.py. This routine makes an icosahedron partition of the earth surface, that is using equilateral triangles. n controls the numer of times the triangles will be divided.
The important parameters are:

* -g GRD            Netcdf grid
* -b BOUNDS         Bounds: /min_lon/max_lon/min_lat/max_lat
* -n N_SUBDIVISION  number of subdivisions of the global icosahedron. Triangles edge lengths are
* n=6  57.66 km
* n=7  28.83 km
* n=8  14.415 km
* n=9  7.208 km
* n=11  4.870 km
* -d DEPTH_RANGE*   Depth range min_depth/max_depth in km

Depending on n and the size, your selected region and your computer speed, this process might take some time.

::

	pyeq_parametrize_curve_surface_triangles.py -g sam_slab1.0_clip_shifted.grd -b /-81.3/-78.9/0.2/0.7 -n 11 -d 15/40 --verbose -e aftershocks_18_may_2016



You should see various files created. Among them, aftershocks_18_may_2016_geometry.dat is the most informative and should look like this:

::

	#     rdis_long    rdis_lat rdis_depth rdis_length rdis_width  rdis_area ratio_rdis_tdis    strike        dip centroid_long centroid_lat centroid_depth tdis_long1  tdis_lat1 tdis_depth1 tdis_long2  tdis_lat2 tdis_depth2 tdis_long3  tdis_lat3 tdis_depth3  tdis_area
	0000  -79.78222    0.46351      -28.49        2.46       2.46      6.03            1.01     18.77      20.79     -79.76887      0.47066         -28.48  -79.75219    0.48067      -29.02  -79.77126    0.45063      -28.66  -79.78317    0.48067      -27.78      6.07 
	0001  -79.77983    0.48354      -28.31        2.46       2.46      6.03            1.01     18.76      20.79     -79.76648      0.49068         -28.31  -79.76410    0.51071      -28.14  -79.75219    0.48067      -29.02  -79.78317    0.48067      -27.78      6.07 
	.
	.

This file includes the geometry information then required for calculating the Green's function either for rectangle dislocations, source points or triangles dislocation elements.
You can display the geometry in QGIS by creating a shapefile.

::

	pyacs_qgis_geometry2polygon.py -dat aftershocks_18_may_2016_geometry.dat --t

.. figure:: ../../images/qgis_geometry.png
    :width: 800px
    :align: center
    :height: 500px
    :alt: alternate text
    :figclass: align-center

    geometry file seen with QGIS 


Making the Green's functions
----------------------------

The next step is to calculate the Green's functions relating the unit slip to the observed displacements/velocities.

Then, assuming  run

::

	cd green
	pyeq_make_green.py -gps_h ../data/mw67_coseismic.dat -g ../geometry/aftershocks_18_may_2016_geometry.npy --tde --verbose -e green_aftershocks_18_may_2016

In addition to the Green's functions, so additional files are created. The most important is green_aftershocks_18_may_2016_input.npy . 
This files includes all the information required both calculating regularization matrices and building the linear system for the inversion.

Running the inversion
---------------------

We simply run a static inversion in the models directory by running:

::

	cd models
	pyeq_static_inversion.py -input_npy ../green/green_aftershocks_18_may_2016_input.npy -dc 20.0 -sigma 2000. -m0 0.0 -rake_type 'Euler' -rake_value /-179.551/82.858/0.4265/inverse -rake_constraint 0  --verbose --debug  -e m67	

Some explanation about the various options are explained in the help

::

	pyeq_static_inversion.py -h
	
::

	Required parameters are:
    - sigma      : regularization value. 1/sigma**2 is the weight given to regularization (Cm-1) (mm or mm/yr)
    - max_slip   : 0 no constraint, Euler calculated from Euler pole, or any float value
    - m0         : a priori model; any real number [0-1]; 0=null a priori model; 1= a priori max_slip
	optional arguments:
	  -h, --help            show this help message and exit
	  -input_npy INPUT_NPY  input npy file including all information
	  -dc DC                Correlation length in km
	  -sigma SIGMA          sigma in mm or mm/yr ; weight of the regularization
	  -m0 M0                a priori model
	  --c C                 constraints on a specific patch n/min_value/max_value
	  -rake_type RAKE_TYPE  rake_type: Euler or constant
	  -rake_value RAKE_OR_POLE
                        Euler pole /long/lat/w/style (style among inverse, normal, leftlateral,rightlateral) or rake value in degrees
	  -rake_constraint RAKE_CONSTRAINT
                        constraint on rake; 0 means fixed, any positive value means sigma on the complement of the principal rake
	  -max_slip MAX_SLIP    constraints on the maximum slip. 0 means maximum from Euler pole, any negative number means no bounds
	  --verbose, -v         verbose mode
	  --debug               debug mode
	  --save                save G,d,Cd, lower, upper matrices/vectors and stop
	  -e EXPERIMENT         experiment name

Getting the results
-------------------

All results files are put in a directory with automatic naming convention. The above example will create a directory m67_sigma_2000_dc_020_m0_000 with files:

::

	m67_sigma_2000_dc_020_m0_000_model.dat # model displacements/velocities GMT psvelo format
	m67_sigma_2000_dc_020_m0_000_residuals.dat # obs - model displacements/velocities GMT psvelo format
	m67_sigma_2000_dc_020_m0_000_slip_dir.dat # slip direction, projected onto the surface, GMT psvelo format
	m67_sigma_2000_dc_020_m0_000_sol_coupling.dat # if bounded inversion, this file is the coupling, GMT psxy format
	m67_sigma_2000_dc_020_m0_000_sol_slip.dat # inverted slip file, GMT psxy format
	m67_sigma_2000_dc_020_m0_000_sum.dat # summary file
	
Again, the results can be seen with converting the file to shapefiles

::

	pyacs_qgis_psvelo2shapefile.py -gmt m67_sigma_2000_dc_020_m0_000_model.dat 
	pyacs_qgis_psvelo2shapefile.py -gmt ../../data/mw67_coseismic.dat 
	pyacs_qgis_model2polygon.py -dat m67_sigma_2000_dc_020_m0_000_sol_slip.dat -g ../../geometry/aftershocks_18_may_2016_geometry.dat 
	
	
.. figure:: ../../images/qgis_model.png
    :width: 800px
    :align: center
    :height: 500px
    :alt: alternate text
    :figclass: align-center

    Model file seen with QGIS. Blue arrows are observed displacements. Red arrows are modeled displacements. 
 
