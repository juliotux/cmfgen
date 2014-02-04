#!/usr/bin/perl
#
# changes all local files named *_in.dat to *_IN
#
@outlist = split(" ",`ls *_in.dat`);
foreach $outfile (@outlist) {
    $infile = $outfile;
    $infile =~ s/($infile)/\U$1\E/;   #convert string to uppercase
    $infile =~ s/\.DAT//;             #chop off .DAT
    `mv $outfile $infile`;            #rename file on disk
}
