#!/usr/bin/perl -w 
use strict;
my $CkpFileOne = $ARGV[0];
my $CkpFileTwo = $ARGV[1];
my $CheckInfo = $ARGV[2];
my $CheckNum  = $ARGV[3];

&_check_ckpFile ($CkpFileOne, $CheckInfo);
&_check_ckpFile ($CkpFileTwo, $CheckInfo);

exit 0;

sub _check_ckpFile {
	my $file = shift;
	my $info = shift;
	open FILE_HANDLE, $file or die $!;
	while (<FILE_HANDLE>) {
		my @read_field=split;
		next unless $read_field[0] eq $info;
		exit 1 unless $read_field[1] == $CheckNum;
	} close FILE_HANDLE;
}