#!/usr/bin/perl -w 
use strict;
open FILE_HANDLE_CKP, $ARGV[0] or die $!;
while (<FILE_HANDLE_CKP>) {
	chomp;
	my @read_field=split;
	s/\t\d\t/\t$ARGV[3]\t/ if ( $read_field[0] eq $ARGV[1] and $read_field[2] eq $ARGV[2] );
	print "$_\n";
	undef @read_field;
}
close FILE_HANDLE_CKP;