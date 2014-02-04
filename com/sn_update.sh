#!/bin/csh

echo " "
echo "This routine is NOT to be used to run a new TIME step SN model."
echo " "
echo "It is used to update SN data files when re-ruunning a SN model"
echo "based on updated calculations at the previous time steop."
echo " "
echo "Use drad_cpmod to start a new SN model. "
echo " "
echo "This routine should be used as follows"
echo " "
echo "     cpmod ts_10   ts_10_updated"
echo "     $cmfdist/com/sn_update.sh  ts_9_updated ts_10_updated"


# Test that we will not inadvertantly copy to the same directory
# We also test that the second directory does not already contain batch.sh,
#   possibly indicating model alredy exists in that directory.
 
if($1 == $2)then
  echo "Warning: directories are the same"
  echo "Aborting copy command"
  goto  SKIPSN
endif

if(-e $2/batch.sh)then
else
 echo "Warning: batch.sh should exist in the new model directory"
 echo "Read header informtion for this command"
 echo "Aborting copy command"
 exit
 goto  SKIPSN
endif

if(-e $1/JH_AT_CURRENT_TIME)then
  echo " "
else
  echo " "
  echo "Problem: SN files at old time step not available"
  echo " "
  goto  SKIPSN
endif


if(-e $1/JH_AT_CURRENT_TIME)then
  cp $1/JH_AT_CURRENT_TIME $2/JH_AT_OLD_TIME
  echo "    Copyed JH_AT_CURRENT_TIME to JH_AT_OLD_TIME"
else
  echo "    *** Unable to copy JH_AT_CURRENT_TIME"
endif

if(-e $1/JH_AT_CURRENT_TIME_INFO)then
  cp $1/JH_AT_CURRENT_TIME_INFO $2/JH_AT_OLD_TIME_INFO
  echo "    Copyed JH_AT_OLD_TIME_INFO to JDH_AT_OLD_TIME_INFO"
else
  echo "    *** Unable to copy JH_AT_OLD_TIME_INFO"
endif

if(-e $1/CUR_MODEL_DATA)then
  cp $1/CUR_MODEL_DATA $2/OLD_MODEL_DATA
  echo "    Copyed OLD_MODEL_DATA"
else
  echo "    *** Unable to copy OLD_MODEL_DATA"
endif

if(-e $1/SN_HYDRO_FOR_NEXT_MODEL)then
  cp $1/SN_HYDRO_FOR_NEXT_MODEL $2/SN_HYDRO_DATA
  echo "    Copyed SN_HYDRO_DATA"
else
  echo "    *** Unable to copy SN_HYDRO_DATA"
endif


if(-e $1/arnaud_rothenflug.dat)then
  cp $1/arnaud_rothenflug.dat         $2/
  echo " "
  echo " Copyed arnaud_rothenflug.dat"
endif


SKIPSN:

exit

