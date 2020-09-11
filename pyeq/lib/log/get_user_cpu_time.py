
def get_user_cpu_time( ):
    """
    return the user cpu time
    """
    import psutil
    
    
    return psutil.cpu_times().user