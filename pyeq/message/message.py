def message( str , level=0 ):
    """
    print message
    """

    from termcolor import colored

    if level > 0:
        banner = '###############################################################################'
        str = str.upper()
        print(colored(banner,attrs=['bold']))
        print(colored("[PYEQ] %s" % (str),attrs=['bold']))
        print(colored(banner,attrs=['bold']))

    else:
        print("[PYEQ] %s" % (str))

