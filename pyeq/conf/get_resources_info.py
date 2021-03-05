def get_resources_info( model ):


    import pkg_resources
    import platform
    import getpass       
    import os                                                                                                                                                                                     
    import psutil
    import sys
    
    ###################################################################
    # RESOURCES & SYSTEM INFO
    ###################################################################

    import resource

    model.username = getpass.getuser() 
    model.os_version = platform.system()
    model.hostname = platform.node()
    model.wdir = os.getcwd()
    model.n_cpu = psutil.cpu_count( logical=False )   
    model.n_thread = psutil.cpu_count( logical=True )
    model.memory = psutil.virtual_memory().total / 1024 / 1024 / 1024
    model.python_version =  sys.version.split()[0]
    model.pyacs_version  = pkg_resources.get_distribution("pyacs").version
    model.pyeq_version  = pkg_resources.get_distribution("pyeq").version
    model.cmd_line  = " ".join(sys.argv)
    
    return model 
    
    