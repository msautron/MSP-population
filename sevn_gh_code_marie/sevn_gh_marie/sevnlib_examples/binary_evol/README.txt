The cpp file binary_evol.cpp shows a simple example of the usage of SEVN as a library.
The code create a binary and evolve it until a given time is reached.

The script compile_binary_evolve.sh shows how to compile the file.
It assumes that the seven library has been properly installed (with make and make install)
and that the sevenv script is inside the PATH.
If the sevnenv is not inside the bin path you have to add the full path in the script.
In order to use directly the script make it executable (chmod u+x binary_evol.cpp)

To run the compiled code just type
./binary_evol.exe

you can add runtime parameters as in the usual SEVN run