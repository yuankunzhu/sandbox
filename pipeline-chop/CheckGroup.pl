#!/usr/bin/perl -w 
use strict;
my $CkpFile = $ARGV[0];
my $CkpInfo = $ARGV[1];
my $CkpNum  = $ARGV[2];

open FILE_HANDLE, $CkpFile or die $!;
while (<FILE_HANDLE>) {
	my @read_field=split;
	next unless $read_field[0] =~ /^$CkpInfo\./;
	exit 1 unless $read_field[1] == $CkpNum;
} close FILE_HANDLE;

exit 0;
