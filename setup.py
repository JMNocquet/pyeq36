from setuptools import setup

###############################################################################
# VERSION HISTORY
# 0.55.2 on 20230330 : added Green's function computation using cutde (https://github.com/tbenthompson/cutde)
# 0.55.1 on 20230328 : can now read 0-360 longitude files in pyeq_parametrize_curve_surface_triangle.py/get_depth
# 0.55.0 on 20230227 : added import pyacs in pyaks.py. Was causing an error when --cross_validation build with mpck option
# 0.54.9 on 20221106 : adding tensile slip for nikkhoo in pyaq_make_green.py and pyeq.make.green for nikkhoo_tde/rde
# 0.54.8 on 20221010 : New version avec Github update
# 0.54.7 on 20220904 : K-fold cross-validation operational
# 0.54.6 on 20220904 : start implementation of K-fold cross-validation approach for regularization parameters choice
# 0.54.5 on 20220329 : implemented slip resolved on vertices. model.geometry_type = 'TDV added'
#                    : added backward, ts_zero options
# 0.54.4 on 20211229 : new mesh using jigsaw. added script pyeq_mesh_triangles_jigsawpy.py
# 0.54.3 on 20210905 : change in setup.py.
# 0.54.2 on 20210617 : change in plot_model_shp to avoid shapely warning for centroid in geographical coordinates
# 0.54.1 on 20210602 : test to best choose sparse matrix choice in DLO regularization
# 0.54.0 on 20210602 : removing temporal smoothing at start seems to provide mitigate results depending on the case.
#                    : 0.54.0 puts back temporal smoothing at start.
# 0.53.9 on 20210520 : modify temporal smoothing and stf smoothing for the first 2-days (post-seismic)
# 0.53.8 on 20210519 : added stf smoothing
# 0.53.7 on 20210509 : renormalized cumulated slip constraint by 1/nstep
# 0.53.6 on 20210505 : renormalized cumulated slip constraint by 1/nstep**2
# 0.53.5 on 20210430 : added spatial laplacian constraints on the cumulated slip
# 0.53.4 on 20210426 : attempt to add resolution dependent laplacian constraints
# 0.53.3 on 20210422 : bug corrected on Green tensor dimension in nikkhoo_tde
# 0.53.2 on 20210409 : variable rake implemented + interpolation on it
# 0.53.1 on 20210322 : adding rake info in info
# 0.53.0 on 20210115 : first release after major refactoring - in test
# 0.52.2 on 20210111 : refactor conf and message
# 0.52.1 on 20210110 : new laplace regularization in pyeq.regularization.laplace in test
# 0.52.0 on 20201225 : refactoring coulomb package (stress/strain etc)
# 0.51.9 on 20201213 : refactoring initiated. green package, handle strain/stress calculation for meade tde.
# 0.51.8 on 20201208 : add time variable lambda_spatial_smoothing
# 0.51.7 on 20201117 : change in pyeq_plot_kinematics_shp and model.pck save
# 0.51.6 on 20201115 : added L1 inversion through CVXOPT - not operational
# 0.51.5 on 20201112 : new discrete laplacian operator implemented
# 0.51.4 on 20201104 : work on new laplacian_like regularization
# 0.51.3 on 20201102 : added pyeq_insert_model_into_model.py
# 0.51.2 on 20201019 : fix small bugs indicated by sphinx. Added sphinx doc built in setup.py
# 0.51.1 on 20200921 : work on bug plotting on some computers
# 0.51.0 on 20200629 : development version - adding back coulomb calculation
# 0.50.8 on 20200306 : development version for refactoring
# 0.50.7 on 20200218 : correct a long standing bug for Green tensor using tde for the up component. Should only impact pyeq_kinematics where up has not been used so far.
# 0.50.6 on 20190212 : added obs_tensor2sgts in pyeq.lib.obstensor
# 0.50.5 on 20191129 : added pyeq.lib.log.geometry2shp_gmt.py and use of it in pyeq_parametrize_curve_surface_triangles.py
# 0.50.4 on 20191129 : development version for build3 corresponding to the paper Nocquet, PYEQ, in prep. 
# 0.50.3 on 20191010 : development version for new build corresponding to the paper Nocquet, PYEQ, in prep. 
# 0.50.2 on 20190729 : work on pyeq.lib.okada_rde.
# 0.50.1 on 20190430 : corrected bug in pyeq_static_inversion.py when m0 is not 0. 
# 0.50.0 on 20190401 : pyeq moved as a separate project. See pyacs for history in changes
###############################################################################

from sphinx.setup_command import BuildDoc
cmdclass = {'build_sphinx': BuildDoc}

name = 'pyeq'
version = '0.55'
release = '0.55.2'

setup(name = name,
      version = release,
      description='PYEQ: PYACS modeling tools for fault slip inversion',
      long_description='PYACS modeling tools for static and kinematic slip inversion',
 
      classifiers=[
        'Development Status :: 3 - Alpha',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.6',
        'Topic :: Geodesy :: Geophysics',
      ],
      keywords='geodesy GPS earthquake elastic dislocation time series',
      url='',
      author='Jean-Mathieu Nocquet (Geoazur, IRD, CNRS, OCA, Cote d Azur University, France)',
      author_email='nocquet@geoazur.unice.fr',
      license='NA',
      packages=['pyeq',
                'pyeq.green',
                'pyeq.green.make',
                'pyeq.green.meade_tde',
                'pyeq.green.nikkhoo_tde',
                'pyeq.green.edcmp_rde',
                'pyeq.green.okada_rde',
                'pyeq.message',
                'pyeq.coulomb',
                'pyeq.regularization',
                'pyeq.regularization.laplace',
                'pyeq.regularization.damping',
                'pyeq.lib',
                'pyeq.conf',
                'pyeq.date',
                'pyeq.elastic_tensor',
                'pyeq.forward_model',
                'pyeq.lib.geometry',
                'pyeq.gps_time_series',
                'pyeq.log',
#                'pyeq.lib.meade_tde',
                'pyeq.optimization',
                'pyeq.optimization.nnls',
                'pyeq.optimization.wrapper',
                'pyeq.optimization.rake',
                'pyeq.lib.objects',
                'pyeq.obs_tensor',
#                'pyeq.lib.okada_rde',
                'pyeq.plot',
                'pyeq.lib.regularization',
                'pyeq.lib.inversion'],
       scripts=[
           #                # QGIS TRANSLATERS
                'pyeq/scripts/pyeq_static_model_to_qgis.py',
           #                # PYEQ INVERSION
                'pyeq/scripts/pyeq_kinematic_inversion.py',
           #                'pyeq/scripts/pyeq_plot_kinematic_results.py',\
                'pyeq/scripts/pyek_print_result_from_mpck.py',
           'pyeq/scripts/pyeq_plot_kinematics_shp.py',
           'pyeq/scripts/pyeq_static_inversion.py',
           'pyeq/scripts/pyeq_insert_model_into_model.py',
           #                # PYEQ GEOMETRY AND GREEN'S FUNCTIONS
           'pyeq/scripts/pyeq_parametrize_curve_surface_triangles.py',
           'pyeq/scripts/pyeq_mesh_triangles_jigsawpy.py',
           'pyeq/scripts/pyeq_make_rectangular_fault.py',
           'pyeq/scripts/pyeq_make_green.py',
           'pyeq/scripts/pyeq_model_to_disp.py',
           'pyeq/scripts/pyeq_interpolate_model.py',
           #                # PYEQ COULOMB
           'pyeq/scripts/pyeq_coulomb.py',
           #                # PYEQ MISC
                'pyeq/scripts/pyeq_magnitude.py',
           'pyeq/scripts/pyeq_scaling_laws.py',
           'pyeq/scripts/pyeq_kinematic_model_to_disp_time_series.py',
           # development version
           'pyeq/scripts/pyaks.py',

       ],
      # change for cargese summer school 2021
      install_requires=['pyacs>=0.65.5',
                        'progress==1.5',
                        'str2bool==1.1',
                        'termcolor==1.1.0',
                        'psutil==5.8.0',
                        'tqdm==4.59.0',
                        'netcdf4==1.5.7'],
          #'Polygon3',  \
#                        'hdf5',  \
#                        'netCDF4'],

#       install_requires=['ansicolors',
#                         'matplotlib>=3',
#                         'geopandas>=0.7.0',
#                         'pyacs>=0.62.0',
#                         'numpy>=1.17.2',
#                         'progress',
#                         'descartes' ,
#                         'str2bool',
#                         'termcolor'],
#           #'Polygon3',  \
# #                        'hdf5',  \
# #                        'netCDF4'],
      zip_safe=False,

      test_suite='nose.collector',
      tests_require=['nose'],

# for building documentation using sphinx
      # these are optional and override conf.py settings
      command_options={
          'build_sphinx': {
              'project': ('setup.py', name),
              'version': ('setup.py', release),
              'release': ('setup.py', release),
              'source_dir': ('setup.py', 'documentation/source'),
              'build_dir': ('setup.py', 'documentation/build')}}

)
