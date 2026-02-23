The cpp file binary_evol.cpp shows a simple example of the usage of SEVN as a library.
The code create a binary and evolve it until a given time is reached.

To compile the code with CMake

- mkdir build
- cd build
- cmake ..
- make

NOTICE: if the cmake command return an error due to the fact that the SEVN package was not found
	- Check if you installed the sevn libraries
	- If you have installed the sevn libraries is a custom directory add this directory (the one set with -DCMAKE_INSTALL_PREFIX)
	  to the cmake search path with:
	 	 cmake .. -DCMAKE_PREFIX_PATH <path_to_the_SEVN_install_dir>

To run the compiled code just type
./binary_evol.exe

you can add runtime parameters as in the usual SEVN run