#!/usr/bin/perl
use strict;
use warnings;

open (FAI, $ARGV[0]) or die "Can not read file ($ARGV[0])!\n";

while (my $FAIline=<FAI>) {
	chomp $FAIline; 
	my @FAIfields = split /\t/, $FAIline; 
		my $count = 0;
		for (1..$FAIfields[1]) {
		$count++;
		print "$FAIfields[0]\t$count\n";
		}
}
;

