#!/bin/csh
#
echo " "
echo "This program compares .sh files in two directories"
echo "If one directory arguments is supplied, the comparison directory is taken as the pwd"
echo "Two directory arguments may also be supplied"
echo "Output is to Diff_output"

if ($2 =="")then
 echo " "
 echo "Present directory is "
 pwd
 echo " "
endif

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
  set main=*.sh
  set second = $1
else
  set sm_main=$1
  set main = $1/*.sh
  set second = $2
endif
echo " " > Diff_output
echo "Comparing files in" $sm_main "with those in " $second  >> Diff_output
echo " " >> Diff_output


foreach i ($main) 
    echo $i >> Diff_output 
    sdiff -s -w200 $i $second >> Diff_output ;            #rename file on disk
end
#
