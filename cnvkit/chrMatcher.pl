#!/usr/bin/perl
open IN, $ARGV[0];
while (<IN>) {
	chomp;
	@a=split;
	$h{$a[1]}=$a[0];
}
close IN;

open IN2, $ARGV[1];
while (<IN2>) {
	@a=split;
	next if $a[2] eq "NO";
	s/$a[2]/$h{$a[2]}/;
	print;
}
close IN2;