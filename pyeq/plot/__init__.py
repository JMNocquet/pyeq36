###############################################################################
# CUSTOM IMPORT
###############################################################################

__custom_import__ = ['plot_stf', 'plot_time_series', 'plot_model', 'plot_model_shp', 'plot_model_interseismic',
                     'model2shp_gmt', 'model_period2dat', 'make_plot', 'interpolate_model','plot_settings']

for mod in __custom_import__:
    exec('from .' + mod + ' import *')
