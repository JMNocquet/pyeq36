###############################################################################
# CUSTOM IMPORT
###############################################################################

__custom_import__ = ['make_model_spatial_correlation_matrix',\
                     'normalization_factor',\
                     'make_model_covariance_01',\
                     'make_model_covariance_02',\
                     'make_model_covariance_03',\
                     'make_regularization_laplacian',\
                     'decipher_sigma_arg', \
                     'time_to_2D_diff_time_in_days',\
                     'make_sqrt_inverse_corr_exponential_spatial',\
                     'make_regularization_valette_01']

for mod in __custom_import__:
    exec('from .' + mod + ' import *')

