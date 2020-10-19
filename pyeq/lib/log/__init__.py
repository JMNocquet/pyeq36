###############################################################################
# CUSTOM IMPORT
###############################################################################

__custom_import__ = [
    'display_array',
    'get_current_memory_usage',
    'get_process_memory_usage',
    'get_user_cpu_time',
    'print_displacement_fields',
    'print_model_attributes',
    'print_model_tensors_shape',
    'print_modeled_time_series',
    'print_observed_time_series',
    'print_residual_time_series',
    'print_geometry',
    'print_dates',
    'print_offset',
    'print_results',
    'print_slip_time_series',
    'print_slip_model',
    'print_sol_to_slip',
    'print_stats',
    'print_stf',
    'print_conf',
    'print_summary']

for mod in __custom_import__:
    exec('from .' + mod + ' import *')
