#!/usr/bin/env python

"""
Created on Wed Aug 26 10:33:25 2020

@author: Laurent Gilquin
"""

import os
import sys
import textwrap
import warnings
import sysconfig
import builtins
from distutils.command.sdist import sdist
from numpy.distutils.command.build_clib import build_clib
from numpy.distutils.command.build_ext import build_ext


if sys.version_info[:2] < (3, 0):
    raise RuntimeError("Python version >= 3.0 required.")

CLASSIFIERS = """\
Development Status :: 1 - Production
Intended Audience :: Science/Research
License :: OSI Approved :: GNU-CPL v3 and CeCILL v2 licences
Programming Language :: C
Programming Language :: Python :: 3
Topic :: Software Development
Topic :: Scientific/Engineering
Operating System :: Unix
"""

MAJOR = 1
MINOR = 0
MICRO = 0
ISRELEASED = False
VERSION = '%d.%d.%d' % (MAJOR, MINOR, MICRO)


# BEFORE importing setuptools, remove MANIFEST. Otherwise it may not be
# properly updated when the contents of directories change (true for distutils,
# not sure about setuptools).
if os.path.exists('MANIFEST'):
    os.remove('MANIFEST')

# This is a bit hackish: we are setting a global variable so that the main
# umd_package __init__ can detect if it is being loaded by the setup routine, to
# avoid attempting to load components that aren't built yet.  While ugly, it's
# a lot more robust than what was previously being used.
#builtins.__UMD_SETUP__ = True

def get_build_ext_override():
    """
    Custom build_ext command to tweak extension building.
    """
    from numpy.distutils.command.build_ext import build_ext as old_build_ext

    class build_ext(old_build_ext):
        def finalize_options(self):
            super().finalize_options()

            # Disable distutils parallel build, due to race conditions
            # in numpy.distutils (Numpy issue gh-15957)
            if self.parallel:
                print("NOTE: -j build option not supportd. Set NPY_NUM_BUILD_JOBS=4 "
                      "for parallel build.")
            self.parallel = None

        def build_extension(self, ext):
            # When compiling with GNU compilers, use a version script to
            # hide symbols during linking.
            if self.__is_using_gnu_linker(ext):
                export_symbols = self.get_export_symbols(ext)
                text = '{global: %s; local: *; };' % (';'.join(export_symbols),)

                script_fn = os.path.join(self.build_temp, 'link-version-{}.map'.format(ext.name))
                with open(script_fn, 'w') as f:
                    f.write(text)
                    # line below fixes gh-8680
                    ext.extra_link_args = [arg for arg in ext.extra_link_args if not "version-script" in arg]
                    ext.extra_link_args.append('-Wl,--version-script=' + script_fn)

            # Allow late configuration
            hooks = getattr(ext, '_pre_build_hook', ())
            _run_pre_build_hooks(hooks, (self, ext))

            old_build_ext.build_extension(self, ext)

        def __is_using_gnu_linker(self, ext):
            if not sys.platform.startswith('linux'):
                return False

            # Fortran compilation with gfortran uses it also for
            # linking. For the C compiler, we detect gcc in a similar
            # way as distutils does it in
            # UnixCCompiler.runtime_library_dir_option
            if ext.language == 'f90':
                is_gcc = (self._f90_compiler.compiler_type in ('gnu', 'gnu95'))
            elif ext.language == 'f77':
                is_gcc = (self._f77_compiler.compiler_type in ('gnu', 'gnu95'))
            else:
                is_gcc = False
                if self.compiler.compiler_type == 'unix':
                    cc = sysconfig.get_config_var("CC")
                    if not cc:
                        cc = ""
                    compiler_name = os.path.basename(cc.split(" ")[0])
                    is_gcc = "gcc" in compiler_name or "g++" in compiler_name
            return is_gcc and sysconfig.get_config_var('GNULD') == 'yes'

    return build_ext


def get_build_clib_override():
    """
    Custom build_clib command to tweak library building.
    """
    from numpy.distutils.command.build_clib import build_clib as old_build_clib

    class build_clib(old_build_clib):
        def finalize_options(self):
            super().finalize_options()

            # Disable parallelization (see build_ext above)
            self.parallel = None

        def build_a_library(self, build_info, lib_name, libraries):
            # Allow late configuration
            hooks = build_info.get('_pre_build_hook', ())
            _run_pre_build_hooks(hooks, (self, build_info))
            old_build_clib.build_a_library(self, build_info, lib_name, libraries)

    return build_clib


def _run_pre_build_hooks(hooks, args):
    """Call a sequence of pre-build hooks, if any"""
    if hooks is None:
        hooks = ()
    elif not hasattr(hooks, '__iter__'):
        hooks = (hooks,)
    for hook in hooks:
        hook(*args)



def parse_setuppy_commands():
    """Check the commands and respond appropriately.  Disable broken commands.
    Return a boolean value for whether or not to run the build or not (avoid
    parsing Cython and template files if False).
    """
    args = sys.argv[1:]

    if not args:
        # User forgot to give an argument probably, let setuptools handle that.
        return True

    info_commands = ['--help-commands', '--name', '--version', '-V',
                     '--fullname', '--author', '--author-email',
                     '--maintainer', '--maintainer-email', '--contact',
                     '--contact-email', '--url', '--license', '--description',
                     '--platforms', '--classifiers', '--requires']

    for command in info_commands:
        if command in args:
            return False

    # Note that 'alias', 'saveopts' and 'setopt' commands also seem to work
    # fine as they are, but are usually used together with one of the commands
    # below and not standalone.  Hence they're not added to good_commands.
    good_commands = ('develop', 'sdist', 'build', 'build_ext', 'build_py',
                     'build_clib', 'build_scripts', 'bdist_wheel', 'bdist_rpm',
                     'bdist_wininst', 'bdist_msi', 'bdist_mpkg',
                     'build_sphinx')

    for command in good_commands:
        if command in args:
            return True

    # The following commands are supported, but we need to show more
    # useful messages to the user
    if 'install' in args:
        return True

    if '--help' in args or '-h' in sys.argv[1]:
        return False


    # The following commands aren't supported.  They can only be executed when
    # the user explicitly adds a --force command-line argument.
    bad_commands = dict(
        test="""`setup.py test` is not supported. """,
        upload="""
            `setup.py upload` is not supported, because it's insecure.""",
        upload_docs="`setup.py upload_docs` is not supported.",
        easy_install="`setup.py easy_install` is not supported.",
        clean="""`setup.py clean` is not supported.""",
        check="`setup.py check` is not supported.",
        register="`setup.py register` is not supported.",
        bdist_dumb="`setup.py bdist_dumb` is not supported.",
        bdist="`setup.py bdist` is not supported.",
        flake8="`setup.py flake8` is not supported.",
        )
    bad_commands['nosetests'] = bad_commands['test']
    for command in ('upload_docs', 'easy_install', 'bdist', 'bdist_dumb',
                     'register', 'check', 'install_data', 'install_headers',
                     'install_lib', 'install_scripts', ):
        bad_commands[command] = "`setup.py %s` is not supported" % command

    for command in bad_commands.keys():
        if command in args:
            print(textwrap.dedent(bad_commands[command]))
            sys.exit(1)

    # Commands that do more than print info, but also don't need Cython and
    # template parsing.
    other_commands = ['egg_info', 'install_egg_info', 'rotate']
    for command in other_commands:
        if command in args:
            return False

    # If we got here, we didn't detect what setup.py command was given
    warnings.warn("Unrecognized setuptools command ('{}').".format(' '.join(sys.argv[1:])))
    return True


def configuration(parent_package='', top_path=None): 
    from numpy.distutils.misc_util import Configuration

    config = Configuration(None, parent_package, top_path)
    config.set_options(ignore_setup_xxx_py=True,
                       assume_default_configuration=True,
                       delegate_options_to_subpackages=True,
                       quiet=True)

    config.add_subpackage('UMD_package')
    config.add_data_files(('UMD_package', 'LICENSE.txt'))

    return config


def setup_package():

    cmdclass = {'sdist': sdist}

    # Figure out whether to add ``*_requires = ['numpy']``.
    # We don't want to do that unconditionally, because we risk updating
    # an installed numpy which fails too often.  Just if it's not installed, we
    # may give it a try.  See gh-3379.
    build_requires = []
    try:
        import numpy
    except ImportError:  # We do not have numpy installed
        build_requires += ['numpy>=1.16.4']
    install_requires = build_requires
    try:
        import scipy
    except ImportError:  # We do not have scipy installed
        install_requires += ['scipy>=1.3.0']
    try:
        import pyopencl
    except ImportError:  # We do not have scipy installed
        install_requires += ['pyopencl>=2020.2']
    try:
        import matplotlib
    except ImportError:
        install_requires += ['matplotlib>=3.1.0']

#    if len(install_requires) != 0:
#        msg = list(map(lambda x: "Error importing " + x.split(">", 1)[0] +
#                       " requires version >" + x.split(">", 1)[1],
#                       install_requires))
#        
#        raise ImportError(". \n".join(msg))

    metadata = dict(
        name = 'UMD_package',
        version = '1.0',
        author = 'See AUTHORS',
        author_email = 'See AUTHORS',
        description = '''
        post-processing package that performs analysis of structural,
        transport, and thermodynamic properties from ab initio molecular dynamics
        simulations.
        ''',
        license = 'GNU GPL v3, CeCILLv2',
        maintainer = 'Razvan Caracas',
        maintainer_email = 'razvan.caracas@gmail.com',
        url = 'https://github.com/rcaracas',
        cmdclass = cmdclass,
        classifiers = [_f for _f in CLASSIFIERS.split('\n') if _f],
        platforms=["Linux"],
        requires = install_requires,
        python_requires= '>=3.0',
    )

    ## define entry points
    scripts_l = ["analyze_1gofr", "analyze_gofr_semi_automatic", "analyze_msd",
                 "averages", "check_overlap", "fullaverages", "gofrs_umd",
                 "insert_umd_xred", "insert_umd", "msd_all_umd", "msd_clusters_umd",
                 "msd_umd", "QBoxParser", "speciation_umd", "stat_concentrate", "umd2poscar",
                 "umd2xyz", "VaspParser", "vibr_spectrum_umd", "viscosity_umd"]
    entry_points = dict(
            console_scripts = ["{}=UMD_package.{}:main".format(*[scr]*2) for scr in scripts_l],
            )

    # Raise errors for unsupported commands, improve help output, etc.
    run_build = parse_setuppy_commands()

    # This import is here because it needs to be done before importing setup()
    # from numpy.distutils, but after the MANIFEST removing and sdist import
    # higher up in this file.
    from setuptools import setup

    if run_build:
        from numpy.distutils.core import setup
        # Customize extension building
        cmdclass['build_ext'] = get_build_ext_override()
        cmdclass['build_clib'] = get_build_clib_override()
        metadata['configuration'] = configuration
        metadata['zip_safe'] = False

    setup(**{**metadata, **{'entry_points':entry_points}})


if __name__ == '__main__':
    setup_package()
