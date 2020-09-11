from setuptools import setup

###############################################################################
# VERSION HISTORY
# 0.51.0 on 20200629 : development version - adding back coulomb calculation
# 0.50.8 on 20200306 : development version for refactoring
# 0.50.7 on 20200218 : correct a long standing bug for Green tensor using tde for the up component. Should only impact pyeq_kinematics where up has not been used so far.
# 0.50.6 on 20190212 : aaded obs_ensor2sgts in pyeq.lib.obstensor
# 0.50.5 on 20191129 : added pyeq.lib.log.geometry2shp_gmt.py and use of it in pyeq_parametrize_curve_surface_triangles.py
# 0.50.4 on 20191129 : development version for build3 corresponding to the paper Nocquet, PYEQ, in prep. 
# 0.50.3 on 20191010 : development version for new build corresponding to the paper Nocquet, PYEQ, in prep. 
# 0.50.2 on 20190729 : work on pyeq.lib.okada. 
# 0.50.1 on 20190430 : corrected bug in pyeq_static_inversion.py when m0 is not 0. 
# 0.50.0 on 20190401 : pyeq moved as a separate project. See pyacs for history in changes
###############################################################################


setup(name='pyeq',
      version='0.51.0',
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
      packages=['pyeq',\
                'pyeq.lib',\
                'pyeq.lib.conf',\
                'pyeq.lib.date',\
                'pyeq.lib.elastic_tensor',\
                'pyeq.lib.forward_model',\
                'pyeq.lib.geometry', \
                'pyeq.lib.gps_time_series', \
                'pyeq.lib.log', \
                'pyeq.lib.meade_tde',\
                'pyeq.lib.nnls',\
                'pyeq.lib.objects',\
                'pyeq.lib.obs_tensor', \
                'pyeq.lib.okada', \
                'pyeq.lib.plot', \
                'pyeq.lib.regularization'],
       scripts=[\
#                # QGIS TRANSLATERS
                'pyeq/scripts/pyeq_static_model_to_qgis.py',\
#                # PYEQ INVERSION
                'pyeq/scripts/pyeq_kinematic_inversion.py',\
                'pyeq/scripts/pyeq_plot_kinematic_results.py',\
                'pyeq/scripts/pyeq_plot_kinematics_shp.py',\
                'pyeq/scripts/pyeq_static_inversion.py',\
#                # PYEQ GEOMETRY AND GREEN'S FUNCTIONS
                'pyeq/scripts/pyeq_parametrize_curve_surface_triangles.py',\
                'pyeq/scripts/pyeq_make_rectangular_fault.py',\
                'pyeq/scripts/pyeq_make_green.py',\
                'pyeq/scripts/pyeq_model_to_disp.py',\
                'pyeq/scripts/pyeq_interpolate_model.py',\
#                # PYEQ COULOMB
           'pyeq/scripts/pyeq_coulomb.py', \
#                # PYEQ MISC
                'pyeq/scripts/pyeq_magnitude.py',\
                'pyeq/scripts/pyeq_scaling_laws.py',\
                'pyeq/scripts/pyeq_kinematic_model_to_disp_time_series.py',\
                ],
      install_requires=['ansicolors', \
                        'pyacs>=0.62.0',\
                        'numpy>=1.17.2',\
                        'progress',\
                        'geopandas', \
                        'descartes' ,\
                        'str2bool'],
          #'Polygon3',  \
#                        'hdf5',  \
#                        'netCDF4'],
      zip_safe=False,

      test_suite='nose.collector',
      tests_require=['nose'],
)
