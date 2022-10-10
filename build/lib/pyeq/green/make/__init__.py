###############################################################################
# CUSTOM IMPORT
###############################################################################

__custom_import__ = ['meade_tde','nikkhoo_tde','nikkhoo_rde','edcmp_rde','okada_rde']

for mod in __custom_import__:
    exec('from .' + mod + ' import *')

