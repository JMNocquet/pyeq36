###############################################################################
# CUSTOM IMPORT
###############################################################################

__custom_import__ = ['cvxopt']

for mod in __custom_import__:
    exec('from .' + mod + ' import *')
