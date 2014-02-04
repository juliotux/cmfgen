#!/bin/csh
#
echo " "
echo "This program compares Makefiles, fortran files, *.INC in two directories"
echo "If one directory arguments is supplied, the comparison directory is taken"
echo "as the pwd. Two directory arguments may also be supplied"
echo "Output is to Diff_output"

rm -f Diff_output

if ($1 == "")then
  echo "Need to supply at least one directory  argument"
  exit
endif

#
# Get the directories
#
if ($2 == "")then
  set sm_main=$dirstack
  set main="Makefile *.f *.f90  *.INC *.txt"
  set second = $1
else
  set sm_main=$1
  set main = "$1/Makefile $1/*.f $1/*.f90 $1/*.INC $1/*.txt"
  set second = $2
endif

echo "Comparing files in" $sm_main "with those in " $second 
#echo "Comparing files " $main
echo " "

echo " " > Diff_output
echo "Comparing files in" $sm_main "with those in " $second  >> Diff_output
echo " " >> Diff_output


foreach i ($main) 
    echo $i >> Diff_output 
    sdiff -s -w200 $i $second >> Diff_output ;            #rename file on disk
end
