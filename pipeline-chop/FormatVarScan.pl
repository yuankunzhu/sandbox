#!/usr/bin/perl -w 
use strict;

open IN, $ARGV[0] or die;
while (<IN>){
	if (/^#/) {
		print; next;
	} else {
		my @tmp = split;
		$tmp[3] =~ s/\//\,/;
		$tmp[4] =~ s/\//\,/;
		my $out = join ("\t", @tmp);
		print "$out\n";
	}
}close IN;