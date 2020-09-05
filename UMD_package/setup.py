"""
Created on Wed Aug 26 10:33:25 2020

@author: Laurent Gilquin
"""

from os.path import join

def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration, get_numpy_include_dirs

    config = Configuration('UMD_package', parent_package, top_path)

    # gofr wrapper
    config.add_extension('c_gofr',
                         sources = [join('src', 'gofr.c')],
                         depends = [join('src', 'gofr.h')],
                         include_dirs = [get_numpy_include_dirs()],
                         define_macros = [('NPY_NO_DEPRECATED_API',
                                           'NPY_1_7_API_VERSION')]
                         )

    config.make_config_py()
    return config


if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
