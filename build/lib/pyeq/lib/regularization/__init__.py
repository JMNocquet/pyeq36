###############################################################################
# CUSTOM IMPORT
###############################################################################

__custom_import__ = ['make_model_spatial_correlation_matrix',\
                     'normalization_factor',\
                     'make_model_covariance_01',\
                     'make_model_covariance_02',\
                     'make_model_covariance_03',\
                     'make_regularization_laplacian',\
                     'time_to_2D_diff_time_in_days',\
                     'make_sqrt_inverse_corr_exponential_spatial',\
                     'make_regularization_valette_01',\
                     'add_spatial_constraint_average_operator',\
                     'get_inter_subfault_distance',\
                     'make_spatial_average_operator',\
                     'make_temporal_average_operator',\
                     'add_damping_constraint',\
                     'add_spatial_constraint_average_operator',\
                     'add_temporal_constraint_average_operator',\
                     'make_regularization_laplacian_like',
                     'make_discrete_laplace_trimesh',
                     'get_regularization_operators',
                     'normalize_lambda'
                     ]

for mod in __custom_import__:
    exec('from .' + mod + ' import *')

