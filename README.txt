Several scripts use C functions compiled as shared libraries (loaded with ctypes):
    c_UMDprocess   reading of UMD and bonding files (umd_processes_fast.py, used by almost all scripts)
    c_gofr         pair distribution functions (gofr_umd.py)
    c_bonds_full   bonds (bonding_umd.py)
    c_clusters     chemical species (speciation_and_angles.py)

Compile them before running the scripts, in the src folder, where the python scripts are: the libraries must be in the same folder.
The compiled libraries are not distributed with the package.


-If you're on MacOS : .dylib libraries

      gcc -O2 -dynamiclib -o c_UMDprocess.dylib c_UMDprocess.c
      gcc -dynamiclib -o c_bonds_full.dylib c_bonds_full.c
      gcc -dynamiclib -o c_clusters.dylib c_clusters.c
      gcc -O2 -dynamiclib -o c_gofr.dylib c_gofr.c


-If you're on Linux : .so libraries

      gcc -O2 -fPIC -c c_UMDprocess.c
      gcc -fPIC -c c_bonds_full.c
      gcc -fPIC -c c_clusters.c
      gcc -O2 -fPIC -c c_gofr.c
      gcc -shared -o c_UMDprocess.so c_UMDprocess.o
      gcc -shared -o c_bonds_full.so c_bonds_full.o
      gcc -shared -o c_clusters.so c_clusters.o
      gcc -shared -o c_gofr.so c_gofr.o


-If you're on Windows : .dll libraries (e.g. with MinGW gcc), the same with .dll instead of .so


The commands are also in compilC-mac.txt (MacOS and Linux) and compilC-windows.txt.
