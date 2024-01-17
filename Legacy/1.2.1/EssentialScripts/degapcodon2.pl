#!/usr/bin/env perl

use strict;
use warnings;

our $header;
our $newseq;
our $line;
while ( $line = <> ) {
    chomp($line);
    my $seqlength = length ($line);

    if ($line =~  m/^>/ ) {
        $header = $line;

    }
    else{
        while ( $line =~ /(...)/g ) {
            if ( $1 !~ /NNN/ ) {
               $newseq .= $1;
            }
        }
    my $NUCcount = $newseq =~ tr/ACTG//;
    my $ratio = $NUCcount/$seqlength;
 	
    if ($ratio > 0.3) {
        print "$header\n";
        print "$newseq\n";}
        $newseq="";
    }
}

