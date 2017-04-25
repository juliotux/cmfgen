#!/bin/csh

#
# Use the procedure to copy files required to generate a new model
# For models of normal stars (i.e., not SN).
#

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

echo " "
echo " Copyed startup files"

if(-e $1/RVSIG_COL)then
  cp $1/RVSIG_COL          $2/
  echo " "
  echo " Copyed RVSIG_COL: This file may need to be edited for a new model"
endif

if(-e $1/HYDRO_DEFAULTS)then
  cp $1/HYDRO_DEFAULTS     $2/
  echo " "
  echo " Copyed HYDRO_DEFAULTS: needed for calculations of the hydro structure"
  echo " Remember to edit the number of iterations"
endif

if(-e $1/ROSSELAND_LTE_TAB)then
  cp $1/ROSSELAND_LTE_TAB  $2/
  echo " "
  echo " Copyed ROSSELAND_LTE_TAB: needed for calculations of the hydro structure"
endif

if(-e $1/RDINR)then
  cp $1/RDINR $2/
  echo "    Copyed RDINR"
endif

if(-e $1/ADJUST_R_DEFAULTS)then
  cp $1/ADJUST_R_DEFAULTS $2/
  echo "    Copyed ADJUST_R_DEFAULTS"
endif

#If not a supernovae model, we exit.

if(-e $1/SN_HYDRO_DATA || -e $1/JH_AT_OLD_TIME)then
  echo " "
  echo " SN data files available "
  echo " "
  echo "This routine is NOT to be used to run a new TIME step SN model."
  echo "It can be used to copy SN output from a model so that the same"
  echo "  time step can be run in another directory."
  echo " "
  echo "Use drad_cpmod to start a new SN model. "
  echo " "
else
  goto  SKIPSN
endif

if(-e $1/NUC_DECAY_DATA)then
  cp $1/NUC_DECAY_DATA $2/
  echo " Copying SN data"
  echo "    Copyed NU_DECAY_DATA"
else
  echo "    *** Unable to copy NU_DECAY DATA"
endif

if(-e $1/JH_AT_OLD_TIME)then
  cp $1/JH_AT_OLD_TIME $2/
  echo "    Copyed JH_AT_OLD_TIME"
else
  echo "    *** Unable to copy JH_AT_OLD_TIME"
endif

if(-e $1/JH_AT_OLD_TIME_INFO)then
  cp $1/JH_AT_OLD_TIME_INFO $2/
  echo "    Copyed JH_AT_OLD_TIME_INFO"
else
  echo "    *** Unable to copy JH_AT_OLD_TIME_INFO"
endif

if(-e $1/OLD_MODEL_DATA)then
  cp $1/OLD_MODEL_DATA $2/
  echo "    Copyed OLD_MODEL_DATA"
else
  echo "    *** Unable to copy OLD_MODEL_DATA"
endif

if(-e $1/SN_HYDRO_DATA)then
  cp $1/SN_HYDRO_DATA $2/SN_HYDRO_DATA
  echo "    Copyed SN_HYDRO_DATA"
else
  echo "    *** Unable to copy SN_HYDRO_DATA"
endif

if(-e $1/arnaud_rothenflug.dat)then
  cp $1/arnaud_rothenflug.dat  $2/
  echo " "
  echo " Copyed arnaud_rothenflug.dat"
endif

if(-e $1/IT_SPECIFIER)then
  cp $1/IT_SPECIFIER  $2/
  echo " "
  echo " Copyed IT_SPECIFIER"
endif

SKIPSN:

cd $2

# Now rename the *OUT files to *IN

out2in
echo " "
echo " Renamed *OUT files to *_IN"
echo " "

