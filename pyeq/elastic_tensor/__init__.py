###############################################################################
# CUSTOM IMPORT
###############################################################################

__custom_import__ = ['exclude_dislocation','npz','shrink','check_obs_vs_green','make_green']

for mod in __custom_import__:
    exec('from .' + mod + ' import *')
