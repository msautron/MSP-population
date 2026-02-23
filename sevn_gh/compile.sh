#!/usr/bin/env bash


SEVN="<Insert absolute SEVNpath>" #Complete path to the SEVN folder



# Usage info
show_help()
{
echo "
        Usage: [-p SEVNPATH] [-j] [-m] [-c] [-u] [-e]

        -p SEVNPATH     Path to the SEVNfolder containing the main CMakeLists.txt.
                        Alternatively the path can be set:
                          - Internally in the script setting the variable SEVN
                          - As environment variable SEVNPATH
                        If more than one of these variable are set the priority is:
                        -p argument, SEVN variable set in the script, environment variable
                        At least one of these variable has to be set.
        -j              Allow parallel compilation
        -m              Compile only, without running Cmake.
        -c              Clean only, run Make clean
        -u              Update the repository calling git pull before the compilation

        For Example: ./compile.sh -p Downloads/SEVN or ./compile.sh -p Downloads/SEVN -m (if the Cmake command has been already run)

        -H              Help
"
}

#Get options
UPDATE=false
ONLYMAKE=false
CLEAN=false
PARALLELMAKE=false
while getopts p:jumch flag
do
    case "${flag}" in
        p) SEVN=${OPTARG};;
        u) UPDATE=true;;
        m) ONLYMAKE=true;;
        c) CLEAN=true;;
        j) PARALLELMAKE=true;;
        h)
          show_help
          exit
          ;;
        \?)
            echo "Invalid option: -$OPTARG. Use -h flag for help."
            exit
            ;;
    esac
done



#Exit the script if some command fail
#set -e
RED="\033[1;31m"
GREEN="\033[1;32m"
BLUE="\033[1;34m"
NC="\033[0m" # No Color


#CHeck
if [ "$SEVN" = "<Insert absolute SEVNpath>" -a "$SEVNPATH" != "" ];
then
    SEVN=$SEVNPATH #Complete path to the SEVN folder
elif [ "$SEVN" = "<Insert absolute SEVNpath>" ];
then
    echo -e "\n${RED} ERROR: SEVN path is not set in the script and SEVNPATH is empty. You have to define one of the two${NC}"
    exit 1
fi


#Check pull update
if [ "$UPDATE" = true ]; then
  echo -e "\n${BLUE}-- Updating repository${NC}"
  sleep 2
  cd $SEVN/
  git pull || return
fi


#Cmake and Make
#--------------------------------------------------------------
#Possible Cmake options (defaults are marked with *)
CMAKE="cmake"

if [ "$PARALLELMAKE" != true ]; then
  MAKE="make"
else
  MAKE="make -j"
fi


#Use Cmake to make the Makefile
if [ "$ONLYMAKE" != true -a "$CLEAN" != true ]; then

  #CMAKE OPTIONS
  DBIN="ON" #Compile executable with binary stellar evolution - [ON*] [OFF]
  DSIN="ON" #Compile executable with single stellar evolution routine only - [ON*] [OFF]
  DH5="OFF" #Enable HDF5 format for output files - [ON] [OFF*]
  DEBUG="OFF" #Enable (verbose) debug information - [ON] [OFF*]
  TEST="OFF" #Compile unit tests - [ON] [OFF*]
  DOCS="OFF" #Build documentation - [ON] [OFF*]


  PSTRING="$PSTRING -Dbin=$DBIN"
  PSTRING="$PSTRING -Dsin=$DSIN"
  PSTRING="$PSTRING -Dh5=$DH5"
  PSTRING="$PSTRING -Ddebug=$DEBUG"
  PSTRING="$PSTRING -Dtest=$TEST"
  PSTRING="$PSTRING -Ddoc=$DOCS"



  echo -e "\n${BLUE}-- Generating Makefile for the SEVN code${NC}"
  sleep 2
  rm -rf $SEVN/build
  mkdir $SEVN/build
  cd $SEVN/build || return

  CMAKE_LAUNCH="$CMAKE .. $PSTRING"
  eval $CMAKE_LAUNCH
  retVal=$?

  if [ $retVal -eq 0 ]; then
      echo -e "${GREEN}-- Cmake execution completed\n\n${NC}"
  else
      echo "${RED}-- Cmake execution ERROR${NC}"
      exit 1
  fi
  #--------------------------------------------------------------


  sleep 2


  echo -e "${BLUE}Compiling SEVN ${NC}"
  sleep 2

fi


#Clean
if [ "$CLEAN" = true ]; then
   MAKE_LAUNCH="$MAKE clean"
else
   MAKE_LAUNCH="$MAKE"
fi

#Make or clean
cd $SEVN/build || return
eval $MAKE_LAUNCH
retVal=$?

if [ $retVal -eq 0 ]; then
    echo -e "${GREEN}[100%] SEVN compilation/cleaning was successful${NC}"
else
    echo "${RED}[XXX] SEVN compilation/cleaning has FAILED${NC}"
    exit 1
fi
cd -
#--------------------------------------------------------------



if [ "$CLEAN" != true ]; then
echo -e "${BLUE}\n**************************************************************************************************************"
echo -e "SEVN has been successfully compiled"
echo -e "We recommend to set all the relevant parameters in one the run script inside the run_scripts folder and to run \nSEVN through the same script (e.g. ./run.sh)."
echo -e "Have fun!"
echo -e "**************************************************************************************************************${NC}"
fi




