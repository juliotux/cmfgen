#!/usr/bin/perl
#
# changes all local files named *OUT to *_in.dat
#
@outlist = split(" ",`ls *OUT`);
foreach $outfile (@outlist) {
    $infile = $outfile;
    $infile =~ s/OUT/_IN/;      #change *OUT to *_IN
    `mv $outfile $infile`;            #rename file on disk
}
