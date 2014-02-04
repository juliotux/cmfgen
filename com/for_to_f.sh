#!/usr/bin/perl
#
# changes all local files named *OUT to *_in.dat
@outlist = split(" ",`ls *.for`);
foreach $outfile (@outlist) {
    $infile = $outfile;
    $infile =~ s/($infile)/\L$1\E/;   #convert string to lowercase
    $infile =~ s/\.for/\.f/;          #change *.for to *.f
    `mv $outfile $infile`;            #rename file on disk
}
@outlist = split(" ",`ls *.for~`);
foreach $outfile (@outlist) {
    $infile = $outfile;
    $infile =~ s/($infile)/\L$1\E/;   #convert string to lowercase
    $infile =~ s/\.for/\.f/;          #change *.for~ to *.f~
    `mv $outfile $infile`;            #rename file on disk
}
#
