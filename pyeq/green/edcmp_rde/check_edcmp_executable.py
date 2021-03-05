def check_edcmp_executable():
    """Check whether edcmp is in PATH."""

    from distutils.spawn import find_executable

    return find_executable('edcmp') is not None
