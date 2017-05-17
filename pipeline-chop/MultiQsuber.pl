#!/usr/bin/perl -w
use strict;

foreach(@ARGV){
	my $count += 0;
	system "/cm/shared/apps/sge/univa/8.1.3/bin/lx-amd64/qsub $_";
	#sleep 1 unless $#ARGV-$count==0;
}
