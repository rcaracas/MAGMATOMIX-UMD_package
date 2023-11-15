A lot of these scripts use c functions under the form of shared libraries. To use them, please follow the following instructions : 

In any case, make sure that the shared library file is in the same folder than the python script that uses it.

To get the shared library :


-If you're on Linux : 
  You'll need shared libraries under the .so form. You can either download them from this repository
  if they're available or create them from the .c file they pertain to.

  To do the latter, type at the command line the following instructions (with 'c_file.c' the name of your c script):
      
      gcc -fPIC -c c_file.c -o c_file.o
        Then
      gcc -shared -o c_file.so c_file.o



-If you're on MacOS :
  You'll need shared libraries under the .dylib form. Again, you can download them from this repository
  if they're available, or use the terminal to create them from the c file :

      gcc -dynamiclib -o c_UMDprocess.dylib c_UMDprocess.c




- Another option is to run:

source compilC.txt


