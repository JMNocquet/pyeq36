def error( str , exit= False ):
    """
    print error message and optionnaly exit
    """

    from colors import red


    if exit:
        print(red("[PYEQ ERROR] %s. Exiting" % (str)))
        import sys
        sys.exit()

    else:
        print(red("[PYEQ ERROR] %s" % (str)))


