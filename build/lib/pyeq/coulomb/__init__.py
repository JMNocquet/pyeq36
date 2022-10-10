###############################################################################
# CUSTOM IMPORT
###############################################################################

__custom_import__ = ['unit_normal_vector',
                     'unit_dip_vector',
                     'unit_strike_vector',
                     'burger',
                     'strain2stress',
                     'stress2coulomb',
                     'stress2normal',
                     'stress2shear',
                     'stress2traction'
                     ]

for mod in __custom_import__:
    exec('from .' + mod + ' import *')
