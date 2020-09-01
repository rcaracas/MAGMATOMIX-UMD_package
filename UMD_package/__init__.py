"""
Created on Fri Aug 28 13:54:45 2020

@author: Laurent Gilquin
"""
#__all__ = [s for s in dir() if not s.startswith('_')]
#
# We first need to detect if we're being called as part of the UMD_package
# setup procedure itself in a reliable manner.
#try:
#    __UMD_SETUP__
#except NameError:
#    __UMD_SETUP__ = False
#
#
#if __UMD_SETUP__:
#    import sys as _sys
#    _sys.stderr.write('Running from UMD_package source directory.\n')
#    del _sys
#else:
#    msg = """Error importing UMD_package: you cannot import the package while
#    being in the source directory; please exit the source
#    tree first and relaunch your Python interpreter."""
#    raise ImportError(msg)
