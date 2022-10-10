###############################################################################
# CUSTOM IMPORT
###############################################################################

__custom_import__ = ['message',
                     'warning',
                     'verbose_message',
                     'debug_message',
                     'error',]

for mod in __custom_import__:
    exec('from .' + mod + ' import *')

