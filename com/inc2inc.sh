#!/usr/bin/perl
#
# changes all local files named *.inc to *.INC
#
@outlist = split(" ",`ls *\.inc`);
foreach $outfile (@outlist) {
    $infile = $outfile;
    $infile =~ s/($infile)/\U$1\E/;   #convert string to lowercase
    `mv $outfile $infile`;            #rename file on disk
}
