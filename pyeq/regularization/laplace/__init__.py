"""
Functions for spatial and temporal regularization through Discrete Laplace Operators (DLO)
"""

###############################################################################
# CUSTOM IMPORT
###############################################################################

__custom_import__ = ['add_laplace_cons',
                     'make_dlo_temporal',
                     'make_dlo_trimesh',
                     'make_dlo_trimesh_node',
                     'make_dlo_stf',
                     'make_normal_dlo_space_time',
                     'normalize_dlo_temporal',
                     'normalize_dlo_trimesh_resolution',
                     'normalize_dlo_trimesh']

for mod in __custom_import__:
    exec('from .' + mod + ' import *')
