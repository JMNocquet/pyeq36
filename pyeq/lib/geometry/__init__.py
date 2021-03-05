###############################################################################
# CUSTOM IMPORT
###############################################################################

__custom_import__ = ['make_topology_trimesh','is_triangle_up']

for mod in __custom_import__:
    exec('from .' + mod + ' import *')
