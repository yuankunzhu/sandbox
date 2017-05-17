#!/usr/bin/perl -w
#usage: UpdateLog.pl log_file "running info" 0/1
use strict;
use Switch;

my $now = localtime;
open FILE_HANDLE_LOG, $ARGV[0] or die $!;
if ( -z FILE_HANDLE_LOG ){
	print "$ARGV[1] running @ $now\n" if $ARGV[2] == 1;
} else {
	while(<FILE_HANDLE_LOG>){
		switch ($ARGV[2]) {
			case 1 {
				if ($. == 1){
					print "$ARGV[1] running @ $now\n$_";
					next;
				} else {
					print;next;
				}
			} 
			case 0 {
				if (/^$ARGV[1]/) {
					chomp;print "\#$_ | finished @ $now\n";
					next;
				} else {
					print;next;
				}
			} 
			case 3 {
				if (/^$ARGV[1]/) {
					chomp;print "$_ | failed @ $now\n";
					next;
				} else {
					print;next;
				}
			}
			case 4 {
				if (/^#/) {
					print;next;
				} else {
					print"\#$_";next;
				}
			}
		}
	} close FILE_HANDLE_LOG;
}
