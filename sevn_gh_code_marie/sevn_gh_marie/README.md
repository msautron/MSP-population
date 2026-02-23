![image](logo.png)

SEVN (Stellar EVolution for 𝑁-body) is a rapid binary population synthesis code.
It gets as input the initial conditions of stars or binaries (masses, spin, semi-major axis, eccentricity)
and evolve them.
SEVN calculates stellar evolution by interpolating pre-computed sets of stellar tracks.
Binary evolution is implemented by means of analytic and semi-analytic prescriptions.
The main advantage of this strategy is that it makes the implementation more general and flexible:
the stellar evolution models adopted in sevn can easily be changed or updated just by loading a new set of look-up tables.
SEVN allows to choose the stellar tables at runtime, without modifying the internal structure of the code or even recompiling it.

SEVN is written entirely in C++ (without external dependencies) following the object-oriented programming paradigm.
SEVN exploits the CPU-parallelisation through OpenMP.
The repository contains also SEVNpy, a Python companion  module can be used to easily  access the SEVN backend in Python. 

Additional information on the technical details of  SEVN can be found in the presentation paper ([Iorio et al., 2022](https://ui.adsabs.harvard.edu/abs/2022arXiv221111774I/abstract), see also [Spera et al., 2019](https://ui.adsabs.harvard.edu/abs/2019MNRAS.485..889S/abstract))
and  in the [user guide](resources/SEVN_userguide.pdf).

# News

**[16-04-24] Version 2.8.0 released**

The main change with respect to the other minor versions is that we changed the  Binary kick processes so that 
the SN kick (vkick in class Star) remains always with the same units (km/s). Before this versione they were changed 
to Rsun/yr after a SN explosion in a binary.


**[24-03-24] Bugs fixed!**

The current version contains fixes on minor bugs.
The only one that deserve a mention is a bug on the estimate of the common envelope lambda with the fitting equation
by [Klencki+21](https://ui.adsabs.harvard.edu/abs/2021A%26A...645A..54K/abstract), i.e. when the option _star_lambda -4_ 
is used. **All the other options are not impacted by this bug.**

Anyway, during the bug fixing, we found that when massive (>30 Msun) and metal-rich  (Z/Zsun>0.4) stars 
contract toward the end of their life, the lambda estimate (with option  _star_lambda -4_ and _star_lambda -41_) by
enters in a regime of extrapolations
with respect the domain studied in  [Klencki+21](https://ui.adsabs.harvard.edu/abs/2021A%26A...645A..54K/abstract).
As a consequence, the lambda becomes very small producing  (likely) unphysical very large binding energies (>1E52 ergs).
We are working on finding the best solutions (see [Issue#4](https://gitlab.com/sevncodes/sevn/-/issues/4)), 
**meanwhile, we suggest to not use 
the options _star_lambda -4_ and _star_lambda -41_ for production runs.**


**[29-08-23] New development line released**

The new SEVN development line *zelaous_redgiant* has been publicly released.
The old development line *humble_yellowdwarf* is still available at this [link](https://gitlab.com/sevncodes/sevn/-/tree/SEVN_humble_yellowdwarf?ref_type=heads).

The most import novelties are:

- The release of the first version of the Python companion module SEVNpy. It is included in the SEVN release. 
  SEVNpy contains both utilities to analyse the outputs of the SEVN runs and class and methods to directly
  access the SEVN backend and evolve star and binary with the same performance of the SEVN C++ executables, but with 
  all the flexibility of interpreted code such as Python. See the  SEVNpy section in the 
  userguide [user guide](resources/SEVN_userguide.pdf) and the online documentation at [http://sevn.rtfd.io/](http://sevn.rtfd.io/).
- SEVN Docker image. We added a Docker container including a working version of SEVN in the [Dockerhub](https://hub.docker.com).
  The are two images: [sevndocker](https://hub.docker.com/r/iogiul/sevndocker) including only a compiled version of SEVN  in a Linux environment,
  and [sevnpydocker](https://hub.docker.com/r/iogiul/sevnpydocker) including also a basic Python3 installation including SEVNpy
  and a running Jupyter notebook that can be accessed by the local browser. The Docker images can be run in essentially 
   any operative systems without caring about installation and dependencies (you need just to have a working version of Docker).
  See the Docker section in the [user guide](resources/SEVN_userguide.pdf) for additional details and instructions.
- Added a complete SEVN versioning policy. The version is specified in the format X.Y.Z, where X is the major version
   and will change only after a substantial refactoring of the code or after a significant change of the main SEVN algorithms.
   Y is the minor version and will change after changes that can break the compatibility with analysis tools or change of the 
   input, addition of new input parameters or moderate refactoring. The Z is the patch version and change almost after each commit.
  Moreover, we include also the concept of "development_line", a development line can be considered a development 
  branch that start to deviate from the main public branch for a variety of reason. When the variation are  too divergent 
  with respect to the current main branch, we create a new development_line. Not all the development line will reach 
  the status of new public branch.  In order to get the info about the SEVN version, just run the SEVN executables
  without runtime parameters. This call will return the list of SEVN parameters with their default value and the
  info about the SEVN version.

# Quickstart


## Requirements 

In order to install and use SEVN 
 - C++ compiler (std C++14, but it is compatible with older compilers)
   - SEVN can be compiled with any GNU C++ compiler, but a version >4.8 is warmly  suggested.
   - SEVN works  with the intel C++ compilers from the [Intel One api toolkit](https://www.intel.com/content/www/us/en/developer/tools/oneapi/toolkits.html#gs.ifxtu2)
   - SEVN can be compiled with clang (the default compiler in macOs systems). Notice that the default clang compiler 
   on Macs does not support OpenMP, therefore it is not possible to exploit the SEVN parallelisation.
 - Cmake  version > 2.8  (>3.2 is suggested)

SEVN uses OpenMP to parallelise the runs, if the Cmake compiler will not find a working OpenMP version (check the Cmake 
log), SEVN will run only in the serial mode. 

### Conda environments

We provide two conda environments that already fulfill the SEVN requirements. They are both located in the folder _resources/:

- [conda_sevn_env](https://gitlab.com/sevncodes/sevn/-/wikis/uploads/1985eb47f3bdc1a984cd2daa1ecf8ffb/conda_sevn_env.txt), simple conda environment containing a C++ compiler and Cmake 
- [conda_sevn_python_env_V2.txt](https://gitlab.com/sevncodes/sevn/-/wikis/uploads/531b124d3d696f30d5d5b932fc18aa94/conda_sevn_python_env.txt), this environment contains also a working Python installation (V 3.9) including some basic packges such as numpy, scipy, astropy, matplotlib, pands. Install this environment if you want to use the SEVNpy package.


In order to use them:
 - Install conda (if not alredy done): [https://conda.io/projects/conda/en/latest/user-guide](https://conda.io/projects/conda/en/latest/user-guide) 
 - Create the new env: `conda create —-name <env_name> —-file conda_sevn_env.txt` or `conda create —-name <env_name> —-file conda_sevn_python_env_V2.txt`
 - Active the new env: `conda activate <env_name>`
 - Install and use SEVN as usual

See the [user guide](resources/SEVN_userguide.pdf) (in _resources_) for additional information

## Installation
Cmake is used to compile SEVN, there are two options, using directly Cmake or using a compile script

### Compile script
   - Locate the [compile script](compile.sh)
   - Open it and edit the line `SEVN="<Insert absolute SEVNpath>"` and replace the string with the absolute path to the SEVN folder
   -  Make the script executable `chmod u+x compile.sh`
   -  Execute the script `./compile.sh`
   -  It is possibile to compile SEVN in parallel using `./compile.sh -j`.  **Be Careful**: this will speed up the compilation at the cost of a significant memory usage. If you machine has less than 4 GB of RAM  or if you are running other memory consuming processes do not use this option.
   -  The executables (_sevnB.x_ and _sevn.x_) will be located in _build/exe_
   -  The script has other runtime options, use ./compile.sh -h to list them

### Cmake 
   - Create a build folder (removing and already present build folder) `rm -r build; mkdir build`  
   - Enter in the build folder and execute cmake `cd build; cmake ..`
   - run make `make`
   - It is possibile to compile SEVN in parallel using `make -j`. **Be Careful**: this will speed up the compilation at the cost of a significant memory usage. If you machine has less than 4 GB of RAM  or if you are running other memory consuming processes do not use this option.
   - The executables (_sevnB.x_ and _sevn.x_) will be located in _build/exe_

## Run

The SEVN compilation produces two executable _sevn.x_ to simulate the single stellar evolution of a list of star, and
_sevnB.x_ to simulate the binary evolution of a list of binaries. 

The two executable can be run with a list of runtime parameters. All of them are optional except for the parameter 
-list containing the path of the input file storing the systems to evolve, e.g. 

`./sevnB.x -list ../../listBin.dat`

The available runtime parameters are described in the [user guide](resources/SEVN_userguide.pdf)

### Run scripts 

Writing a long list of runtime parameters can be tedious. For this reason we provide two run script 
template that can be used to quickly set a SEVN run. 
The two scripts are located in the [run_scripts folder](run_scripts):

   - [run_sse.sh](run_scripts/run_sse.sh) to simulate the single stellar evolution of a list of star
   - [run.sh](run_scripts/run.sh) to simulate the binary evolution of a list of binaries

In order to use the scripts:

   - Open the script and edit the line `SEVN="<Insert absolute SEVNpath>"`  replacing the string with the absolute path to the SEVN folder
   - Edit the line `LISTBIN="${SEVN}/run_scripts/listBin.dat"`  replacing the string  "${SEVN}/run_scripts/listBin.dat" with the path to the initial condition file
   - The script already contains all the runtime options. All of the are set to their default values. 
   - Make the script executable `chmod u+x compile.sh`
   - Run the script `./run.sh` or `./run_sse.sh`


# Documentation

The SEVN  [user guide](resources/SEVN_userguide.pdf) contains a general overview of SEVN, 
information on compilation and running, information on the runtime parameters and stellar tables, 
decription of the SEVN outputs and  some hints and tips to analyse them. 

A more [technical documentation](doc/doxyfiles/html/index.html)  on the code  and its components has been generated with Doxygen and can be found
in the [doc/doxyfiles](doc/doxyfiles) folder.

# Support
If you have problem in running SEVN or you find some weird behaviour/bugs
please open an issue on the Gitlab repository. 
Please use the Gitlab issue tracker also  to ask  for 
requiring additional features for the next SEVN versions. 

For any other questions/doubts/contributions send an email 
to: [giuliano.iorio.astro@gmail.com](mailto:giuliano.iorio.astro@gmail.com) 
and/or  [sevnpeople@gmail.com](mailto:sevnpeople@gmail.com)

# Contributing
SEVN is a general-purpose population-synthesis code and our idea is to create a network of users/developers. 
So, any contributions is highly encouraged. 
If you think it could be useful to add new processes/features in SEVN or you need it for a scientific 
project, please open an issue on Gitlab or contact the main SEVN team (see [Support section](#support)).

We are working to create a SEVN-guide for developers that will help to directly modify and extend the code.
The ideal way to contribute to the code is to fork it in your gitlab account and then send a merge-request when 
your version is ready to be included in the main SEVN repository.
Please, contact us (see [Support section](#support)) for any doubts. 


# Authors and acknowledgment

The original version of SEVN was developed by  Mario Spera, Michela Mapelli and Alessandro Alberto Trani. 

The current updated version is developed and maintained by Giuliano Iorio (main developer).
The SEVN core team includes  Guglielmo Costa, Gaston Escobar, Erika Korb, Michela Mapelli, Mario Spera, Cecilia Sgalletta. 

The developers thank all the people in the DEMOBLACK group for all the valuable comments and suggestions during the code development. 


# License
MIT License

Copyright (c) 2022 Giuliano Iorio, Michela Mapelli, Mario Spera, Guglielmo Costa

See file LICENSE in the repository
# Binary-pulsars-evolution
# Binary-pulsars-evolution
