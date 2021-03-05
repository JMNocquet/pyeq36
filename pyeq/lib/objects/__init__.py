###############################################################################
# CUSTOM IMPORT
###############################################################################

__custom_import__ = ['pyeq_model','topology_trimesh']

for mod in __custom_import__:
    exec('from .' + mod + ' import *')

