#!/bin/csh

#
# Use this procedure to cp files required to generate
# a SN model with a time-dependent radiation field.
#

echo " This procdure is used to copy the files needed to generate a NEW SN mode"
echo " The routine will not copy time depndent SN files if the input and output directories"
echo " are the same."
echo " "

# Test that we will not inadvertantly copy to the same directory
# We also test that the second directory does not already contain batch.sh,
#   possibly indicating model alredy exists in that directory.
 
if($1 == $2)then
  echo "Warning: directories are the same"
  echo "If you wish to continue enter yes"
  set answer=$<
  switch ($answer)
  case [yY][eE][sS]:
    echo "Continuing with the copy command"
    echo "Moving GREY_SCL_FAC_IN to GREY_SCL_FAC_SAVE"
    echo "You may need to rename this file back to GREY_SCL_FAC_IN"
    echo "Continuing with the copy command"
    cp $2/GREY_SCL_FAC_IN $2/GREY_SCL_FAC_SAVE
    breaksw
  default:
    echo "Aborting copy command"
    exit
  endsw
else
  if(-e $2/batch.sh)then
    echo "Warning: batch.sh exists in the new model directory"
    echo "Enter yes if you wish to continue the copy"
    set answer=$<
    switch ($answer)
      case [yY][eE][sS]:
      echo "Continuing with copy command"
      breaksw
    default:
      echo "Aborting copy command"
      exit
    endsw
  endif
endif

cp $1/batch.sh           $2/
cp $1/IN_ITS             $2/
cp $1/VADAT              $2/
cp $1/*OUT               $2/
cp $1/GAMMAS             $2/GAMMAS_IN
cp $1/MODEL_SPEC         $2/

# Copy files required for SN models with a time dependent radiation field.

if($1 == $2)then
 echo "WARNING"
 echo "Directories are the same"
 echo "Will not copy TIME dependent SN data'
else
  cp $1/JH_AT_CURRENT_TIME             $2/JH_AT_OLD_TIME
  cp $1/JH_AT_CURRENT_TIME_INFO        $2/JH_AT_OLD_TIME_INFO
  cp $1/NUC_DECAY_DATA                 $2/
  cp $1/SN_HYDRO_FOR_NEXT_MODEL        $2/SN_HYDRO_DATA
  cp $1/CUR_MODEL_DATA                 $2/OLD_MODEL_DATA
endif

#
# Copy files that do not necessarily exist for all models.
# 

if(-e $1/ADJUST_R_DEFAULTS)then
  cp $1/ADJUST_R_DEFAULTS              $2/
  echo " "
  echo " Copyed ADJUST_R_DEFAULTS: This file wll need to be edited if you wish the R grid"
  echo " to be revised during the run."
endif

if(-e $1/arnaud_rothenflug.dat)then
  cp $1/arnaud_rothenflug.dat         $2/
  echo " "
  echo " Copyed arnaud_rothenflug.dat"
endif

if(-e $1/IT_SPECIFIER)then
  cp $1/IT_SPECIFIER  $2/
  echo " "
  echo " Copyed IT_SPECIFIER"
endif

if(-e $1/RDINR)then
    echo " "
    echo " Warning: RDINR exist in old directory"
    echo " Enter y if you wish to copy this file"
    set answer=$<
    switch ($answer)
      case [yY]:
      cp $1/RDINR              $2/
      echo " Copyed RDINR"
      breaksw
    default:
      breaksw
    endsw
    echo " "
endif
# Change to the new model directory.

cd $2

# Rename the *OUT files to *IN

out2in

echo "You may need to delete the GREY_SCL_FAC_IN file since grey scaling may not work."

rm -f GREY_SCL_FAC_IN
