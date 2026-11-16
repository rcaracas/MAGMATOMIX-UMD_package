#!/usr/bin/env python3
###
#AUTHORS: RAZVAN CARACAS
###
#Compiles the four C libraries of the package (shared libraries loaded with ctypes) into the folder of this script, with the extension of
#the operating system : .dylib (macOS), .so (Linux), .dll (Windows, with the gcc of MinGW-w64). Run it once after downloading the package :
#   python3 compile_libraries.py          (python on Windows)
#On Windows without gcc, take the compiled .dll of the GitHub release instead (umd-c-libraries-windows.zip), unzipped into this folder.

import sys, os, platform, subprocess, shutil

HERE = os.path.dirname(os.path.abspath(__file__))
LIBRARIES = (('c_UMDprocess',['-O2']),('c_bonds_full',[]),('c_clusters',[]),('c_gofr',['-O2']))   #name, optimization (as in README.txt)


def main(argv):
    compiler = os.environ.get('CC','gcc')
    if '-h' in argv:
        print('compile_libraries.py : compiles c_UMDprocess, c_bonds_full, c_clusters, c_gofr into the folder of the scripts (src),')
        print('as .dylib (macOS), .so (Linux) or .dll (Windows). The compiler is gcc, or the one of the environment variable CC.')
        print('Windows : gcc of MinGW-w64 (e.g. from https://www.msys2.org, or winget install BrechtSanders.WinLibs.POSIX.UCRT) ; without a')
        print('compiler, unzip the compiled libraries of the GitHub release (umd-c-libraries-windows.zip) into the src folder.')
        sys.exit()
    if shutil.which(compiler) is None:
        print('No C compiler',compiler,'found.')
        if platform.system() == 'Windows':
            print('Install MinGW-w64 gcc, or unzip the compiled libraries of the GitHub release (umd-c-libraries-windows.zip) into',HERE)
        elif platform.system() == 'Darwin':
            print('Install the command line tools of Xcode : xcode-select --install')
        else:
            print('Install gcc, e.g. sudo apt install gcc')
        sys.exit(1)
    system = platform.system()
    if system == 'Darwin':
        ext, flags = '.dylib', ['-dynamiclib']
    elif system == 'Windows':
        ext, flags = '.dll', ['-shared','-static-libgcc']               #no dependence on the DLLs of MinGW : the libraries can be copied
    else:
        ext, flags = '.so', ['-shared','-fPIC']
    failed = []
    for name, opt in LIBRARIES:
        target = os.path.join(HERE,name+ext)
        command = [compiler]+opt+flags+['-o',target,os.path.join(HERE,name+'.c')]
        print(' '.join(command))
        if subprocess.run(command).returncode != 0 or not os.path.isfile(target):
            failed.append(name)
    if failed:
        print('Compilation failed :',', '.join(failed))
        sys.exit(1)
    print('The four libraries are compiled in',HERE)


if __name__ == '__main__':
    main(sys.argv[1:])
