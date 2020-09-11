"""
Exception class for pyeq
"""
from colors import red

###############################################################################
class PyeqError(Exception):
###############################################################################
    pass

###############################################################################
class PyeqNpzReadError(PyeqError):
###############################################################################

    def __init__(self , method_name , path_npz):
        self.method_name = method_name    
        self.path_npz = path_npz    
    
    def __str__(self):
        return_str = red( "[PYACS ERROR] Error in %s. Could not np.load npz file created by pyeq_make_green.py : %s" % 
                          ( self.method_name, self.path_npz ) ) 
        return( return_str )

###############################################################################
class PyeqNpzFormatError(PyeqError):
###############################################################################

    def __init__(self , method_name , path_npz):
        self.method_name = method_name    
        self.path_npz = path_npz    
    
    def __str__(self):
        return_str = red( "[PYACS ERROR] Format error in %s. Could not properly decode npz file created by pyeq_make_green.py : %s" % 
                          ( self.method_name, self.path_npz ) ) 
        return( return_str )
