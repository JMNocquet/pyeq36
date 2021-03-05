
def get_current_memory_usage( ):
    """
    return the current memory usage
    """
    import psutil
    
    byte_to_Gb = 1024 * 1024 * 1024
    
    return psutil.virtual_memory().used / byte_to_Gb