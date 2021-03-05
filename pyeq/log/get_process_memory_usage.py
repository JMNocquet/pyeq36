def get_process_memory_usage( ):
    """
    return the memory usage
    note: accounts for different operating systems
    """
    import resource
    import platform
    
    if platform.system().lower() == 'darwin':
        byte_to_Gb = 1024 * 1024 * 1024
    
    if platform.system().lower() == 'linux':
        byte_to_Gb = 1024 * 1024 

    if platform.system().lower() == 'windows':
        byte_to_Gb = 1024 * 1024 
    
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / byte_to_Gb