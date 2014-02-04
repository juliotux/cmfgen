#!/bin/csh

#
# Simple script to cycle through all the sub-directories
# and compare Makefile, *.f, and *.INC files in two
# different CMF distribution lists.
#

echo " "
echo "This program compares all Makefiles, Fortran files, and *.INC files"
echo "in 2 CMFGEN directory structures. If one directory arguments is supplied,"
echo "the comparison directory is taken as the pwd. Two directory arguments may"
echo "also be supplied. Output is to Diff_sum. Diff_output is corrupted."
echo " "
pwd
echo " "

if ($1 == "")then
  echo "Need to supply at least one directory  argument"
  exit
endif

rm -f Diff_sum

if ($2 == "")then
#  set X1=$dirstack
  set X1="."
  set X2=$1
else
  set X1=$1
  set X2=$2
endif

echo " " > Diff_sum
echo "Current directory is:" >> Diff_sum
pwd >> Diff_sum
echo " " >> Diff_sum

$cmfdist/com/com_diff.sh $X1/com $X2/com
cat Diff_output >> Diff_sum
 
$cmfdist/com/main_diff.sh $X1/blas $X2/blas
cat Diff_output >> Diff_sum

$cmfdist/com/main_diff.sh $X1/disp $X2/disp
cat Diff_output >> Diff_sum

$cmfdist/com/main_diff.sh $X1/disp/subs $X2/disp/subs
cat Diff_output >> Diff_sum

$cmfdist/com/main_diff.sh $X1/lpack $X2/lpack
cat Diff_output >> Diff_sum

$cmfdist/com/main_diff.sh $X1/lte_hydro $X2/lte_hydro
cat Diff_output >> Diff_sum

$cmfdist/com/main_diff.sh $X1/main $X2/main
cat Diff_output >> Diff_sum

$cmfdist/com/main_diff.sh $X1/misc $X2/misc
cat Diff_output >> Diff_sum

$cmfdist/com/main_diff.sh $X1/new_main $X2/new_main
cat Diff_output >> Diff_sum

$cmfdist/com/main_diff.sh $X1/new_main/mod_subs $X2/new_main/mod_subs
cat Diff_output >> Diff_sum

$cmfdist/com/main_diff.sh $X1/new_main/subs $X2/new_main/subs
cat Diff_output >> Diff_sum

$cmfdist/com/main_diff.sh $X1/new_main/subs/auto $X2/new_main/subs/auto
cat Diff_output >> Diff_sum

$cmfdist/com/main_diff.sh $X1/new_main/subs/chg $X2/new_main/subs/chg
cat Diff_output >> Diff_sum

$cmfdist/com/main_diff.sh $X1/new_main/subs/two $X2/new_main/subs/two
cat Diff_output >> Diff_sum

$cmfdist/com/main_diff.sh $X1/new_main/subs/non_therm $X2/new_main/subs/non_therm
cat Diff_output >> Diff_sum

$cmfdist/com/main_diff.sh $X1/newsubs $X2/newsubs
cat Diff_output >> Diff_sum

$cmfdist/com/main_diff.sh $X1/obs $X2/obs
cat Diff_output >> Diff_sum

$cmfdist/com/main_diff.sh $X1/pgplt $X2/pgplt
cat Diff_output >> Diff_sum

$cmfdist/com/main_diff.sh $X1/plane $X2/plane
cat Diff_output >> Diff_sum

$cmfdist/com/main_diff.sh $X1/spec_plt $X2/spec_plt
cat Diff_output >> Diff_sum

$cmfdist/com/main_diff.sh $X1/spec_plt/subs $X2/spec_plt/subs
cat Diff_output >> Diff_sum

$cmfdist/com/main_diff.sh $X1/subs $X2/subs
cat Diff_output >> Diff_sum

$cmfdist/com/main_diff.sh $X1/subs/chg $X2/subs/chg
cat Diff_output >> Diff_sum

$cmfdist/com/main_diff.sh $X1/subs/two $X2/subs/two
cat Diff_output >> Diff_sum

$cmfdist/com/main_diff.sh $X1/stark $X2/stark
cat Diff_output >> Diff_sum

$cmfdist/com/main_diff.sh $X1/txt_files $X2/txt_files
cat Diff_output >> Diff_sum

$cmfdist/com/main_diff.sh $X1/tools $X2/tools
cat Diff_output >> Diff_sum

$cmfdist/com/main_diff.sh $X1/unix  $X2/unix
cat Diff_output >> Diff_sum
