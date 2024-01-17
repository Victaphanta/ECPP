#!/usr/bin/perl
use strict;
use warnings;

$/ = '>';
while (<>) {
    chomp;
    s/(.+?\n)(.+)/my $x = $2; $x =~ s|\s+||g; $_ = $x/se or next;
    print ">$1$_\n";
}
