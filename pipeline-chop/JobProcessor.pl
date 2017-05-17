#!/usr/bin/perl -w 
# usage: JobProcessor.pl $OUTDIR/$Sample_ID $Analysis_Ckp Check1 Check2
use strict;
use Switch;
use FindBin '$Bin';
use File::Basename;
use lib $Bin;
my $Work_Dir = $ARGV[0];
my $Sample_ID = basename ($Work_Dir);
my ( %Ckp_Number, %Ckp_Type );
my ( $Method_MAP, $Method_SNV, $Method_IDV, $Method_CNV, $Method_SV );
my ( @RunChr_Group );
&_readCheckPoint;

switch ( $ARGV[2] ){
	case 'run'   { &_CheckPointIs_run($ARGV[3]) }
	case 'lane'  { &_CheckPointIs_lane($ARGV[3]) }
	case 'lib'   { &_CheckPointIs_lib($ARGV[3]) }
	case 'map'   { &_CheckPointIs_map($ARGV[3]) }
	case /^bqsr\.G\d+$/  { &_CheckPointIs_bqsr($ARGV[2]) }
	case /^snv\.G\d+$/  { &_CheckPointIs_snv($ARGV[2]) }
	case /^idv\.G\d+$/  { &_CheckPointIs_idv($ARGV[2]) }
	case /^cnv\.G\d+$/  { &_CheckPointIs_cnv($ARGV[2]) }
	case 'sv'    { &_CheckPointIs_sv($ARGV[3]) }
}

close FILE_HANDLE_CKP;

sub _readCheckPoint {
	open FILE_HANDLE_CKP, $ARGV[1] or die "[ERROR] from $ARGV[1]: $!";
	while ( <FILE_HANDLE_CKP> ) {
		chomp;
		my @read_field = split;
		$Ckp_Number{$read_field[0]}{$read_field[2]} = $read_field[1];
		$Ckp_Type{$read_field[0]} = 1;
		$Method_MAP = $read_field[2] if $read_field[0] eq 'map';
		$Method_SNV = $read_field[2] if $read_field[0] eq 'snv';
		$Method_IDV = $read_field[2] if $read_field[0] eq 'idv';
		$Method_CNV = $read_field[2] if $read_field[0] eq 'cnv';
		$Method_SV  = $read_field[2] if $read_field[0] eq 'sv';
		if ( $read_field[0] =~ s/^bqsr\.// ) {
			push @RunChr_Group, $read_field[0];
		} 
		undef @read_field;
	} 
	close FILE_HANDLE_CKP;
}

sub _CheckPointIs_run {
	my $sub_argv_run = shift;
	if ( $sub_argv_run == 1 ) { # $ARGV[3] == 1 eq new run 
		foreach ( keys %{$Ckp_Number{lane}} ){
			my $check_str = "fq2sai $_";
			my $qsub_bash = "step_1-1_fq2sai_$_.sh";
			$qsub_bash = "step_1-2_fq2sai_$_.sh" if /\_2$/;
			&_Check_To_Run ($check_str, $qsub_bash);
		}
	} else { # else is re-run
		foreach my $ckp_type (keys %Ckp_Type) {
			switch ($ckp_type) {
				case 'lane' { map {&_CheckPointIs_lane($_) if s/\_1//} (keys %{$Ckp_Number{lane}}) }
				case 'lib'  { map {&_CheckPointIs_lib($_)} (keys %{$Ckp_Number{lib}}) }
				case 'map'  { &_CheckPointIs_map($Method_MAP) }
				case 'sv'   { &_CheckPointIs_sv($Method_SV) }
				case /^spt\.G\d+$/ { &_CheckPointIs_spt($ckp_type) }
				case /^snv\.G\d+$/ { &_CheckPointIs_snv($ckp_type) }
				case /^idv\.G\d+$/ { &_CheckPointIs_idv($ckp_type) }
				case /^cnv\.G\d+$/ { &_CheckPointIs_cnv($ckp_type) }
			}
		}

	}
	undef $sub_argv_run;
}

sub _CheckPointIs_lane {
	my $sub_argv_lane = shift;
	my $read1 = "$sub_argv_lane\_1";
	my $read2 = "$sub_argv_lane\_2";
	switch ( $Ckp_Number{lane}{$read1} ) {
		case 1 {
			my $check_str = "fq2sai $read1";
			my $qsub_bash = "step_1-1_fq2sai_$read1.sh";
			&_Check_To_Run ($check_str, $qsub_bash);
			#my $sub_return = &_Check_To_Run ($check_str, $qsub_bash);
		}
		case 2 {
			if ( $Ckp_Number{lane}{$read2} == 1 ) {
				my $check_str = "fq2sai $read2";
				my $qsub_bash = "step_1-2_fq2sai_$read2\.sh";
				&_Check_To_Run ($check_str, $qsub_bash);
			} elsif ( $Ckp_Number{lane}{$read2} == 2 ) {
				my $check_str = "sai2sam $sub_argv_lane";
				my $qsub_bash = "step_1-3_sai2sam_$sub_argv_lane.sh";
				&_Check_To_Run ($check_str, $qsub_bash);
			}
		}
	}
	undef $sub_argv_lane;
}

sub _CheckPointIs_lib {
	my $sub_argv_lib = shift;
	my $not_ready_mark = 0;
	switch ( $Ckp_Number{lib}{$sub_argv_lib} ){
		case 1 {
			foreach my $lane_id ( keys %{$Ckp_Number{lane}} ) {
				next unless $lane_id =~ /$sub_argv_lib/;
				if ( $Ckp_Number{lane}{$lane_id} != 4 ) {
					$lane_id =~ s/\_[12]//;
					&_CheckPointIs_lane($lane_id);
					$not_ready_mark += 1;
					next;
				}
			}
			return unless $not_ready_mark == 0;

			my $check_str = "mergeLib $sub_argv_lib";
			my $qsub_bash = "step_1-5_merge_$sub_argv_lib.sh";
			&_Check_To_Run ($check_str, $qsub_bash);
		}
	}
	undef $sub_argv_lib;
	undef $not_ready_mark;
}

sub _CheckPointIs_map {
	my $sub_argv_map = shift;
	my $not_ready_mark = 0;
	switch ( $Ckp_Number{map}{$sub_argv_map} ){
		case 1 {
			foreach my $lib_id ( keys %{$Ckp_Number{lib}} ) {
				if ( $Ckp_Number{lib}{$lib_id} != 2 ) {
					&_CheckPointIs_lib($lib_id);
					$not_ready_mark += 1;
					next;
				}
			}
			return unless $not_ready_mark == 0;
			my $check_str = "mergeSample $Sample_ID";
			my $qsub_bash = "step_1-6_merge_$Sample_ID.sh";
			&_Check_To_Run ($check_str, $qsub_bash);
		}
		case 2 {
			my $check_str = "merge recal bam $Method_MAP";
			my $qsub_bash = "step_1-8_merge_recal_$Sample_ID.sh";
			&_Check_To_Run ($check_str, $qsub_bash);
		}
		case 998 {
			my $check_str = "IndelRealigner $Sample_ID";
			my $qsub_bash = "step_1-7_IndelRealigner_$Sample_ID.sh";
			&_Check_To_Run ($check_str, $qsub_bash);
		}
		case 999 {
			my $check_str = "BaseRecalibrator $Sample_ID";
			my $qsub_bash = "step_1-8_BaseRecalibrator_$Sample_ID.sh";
			&_Check_To_Run ($check_str, $qsub_bash);
		}
	}
	undef $sub_argv_map;
	undef $not_ready_mark;
}

sub _CheckPointIs_bqsr {
	my $sub_argv_grp = shift;
	$sub_argv_grp =~ s/^bqsr\.//;
	switch ( $Ckp_Number{"bqsr.$sub_argv_grp"}{$Method_MAP} ){
		case 1 {
			return unless $Ckp_Number{map}{$Method_MAP} == 2;
			my $check_str = "RealignerTargetCreator $sub_argv_grp";
			my $qsub_bash = "step_1-7-1_RealignerTargetCreator_$Sample_ID.$sub_argv_grp.sh";
			&_Check_To_Run ($check_str, $qsub_bash);
		}
		case 2 {
			my $check_str = "IndelRealigner $sub_argv_grp";
			my $qsub_bash = "step_1-7-2_IndelRealigner_$Sample_ID.$sub_argv_grp.sh";
			&_Check_To_Run ($check_str, $qsub_bash);
		}
		case 3 {
			my $check_str = "BaseRecalibrator $sub_argv_grp";
			my $qsub_bash = "step_1-7-3_BaseRecalibrator_$Sample_ID.$sub_argv_grp.sh";
			&_Check_To_Run ($check_str, $qsub_bash);
		}
		case 4 {
			my $check_str = "PrintReads $sub_argv_grp";
			my $qsub_bash = "step_1-7-4_PrintReads_$Sample_ID.$sub_argv_grp.sh";
			&_Check_To_Run ($check_str, $qsub_bash);
		}
	}
}

sub _CheckPointIs_snv {
	my $sub_argv_grp = shift;
	$sub_argv_grp =~ s/^snv\.//;
	switch ( $Ckp_Number{"snv.$sub_argv_grp"}{$Method_SNV} ){
		case 1 {
			return unless $Ckp_Number{"bqsr.$sub_argv_grp"}{$Method_MAP} == 5;
			my $check_str = "call snv $sub_argv_grp";
			my $qsub_bash = "step_2-1-1_call_snv_$Method_SNV.$sub_argv_grp.sh";
			&_Check_To_Run ($check_str, $qsub_bash);
		}
		case 2 {
			my $check_str = "vqsr snv $sub_argv_grp";
			my $qsub_bash = "step_2-1-2_vqsr_snv_$Method_SNV.$sub_argv_grp.sh";
			&_Check_To_Run ($check_str, $qsub_bash);
		}
		case 3 {
			my $check_str = "somatic snv $sub_argv_grp";
			my $qsub_bash = "step_2-1-3_somatic_snv_$Method_SNV.$sub_argv_grp.sh";
			&_Check_To_Run ($check_str, $qsub_bash);
		}
		case 4 {
			my $check_str = "ann snv $sub_argv_grp";
			my $qsub_bash = "step_2-1-4_ann_snv_$Method_SNV.$sub_argv_grp.sh";
			&_Check_To_Run ($check_str, $qsub_bash);
		}
	}
	undef $sub_argv_grp;
}

sub _CheckPointIs_idv {
	my $sub_argv_grp = shift;
	$sub_argv_grp =~ s/^idv\.//;
	switch ( $Ckp_Number{"idv.$sub_argv_grp"}{$Method_IDV} ){
		case 1 {
			return unless $Ckp_Number{"bqsr.$sub_argv_grp"}{$Method_MAP} == 5;
			my $check_str = "call idv $sub_argv_grp";
			my $qsub_bash = "step_2-2-1_call_idv_$Method_IDV.$sub_argv_grp.sh";
			&_Check_To_Run ($check_str, $qsub_bash);
		}
		case 2 {
			my $check_str = "vqsr idv $sub_argv_grp";
			my $qsub_bash = "step_2-2-2_vqsr_idv_$Method_IDV.$sub_argv_grp.sh";
			&_Check_To_Run ($check_str, $qsub_bash);
		}
		case 3 {
			my $check_str = "somatic idv $sub_argv_grp";
			my $qsub_bash = "step_2-2-3_somatic_idv_$Method_IDV.$sub_argv_grp.sh";
			&_Check_To_Run ($check_str, $qsub_bash);
		}
		case 4 {
			my $check_str = "ann idv $sub_argv_grp";
			my $qsub_bash = "step_2-2-4_ann_idv_$Method_IDV.$sub_argv_grp.sh";
			&_Check_To_Run ($check_str, $qsub_bash);
		}
	}
	undef $sub_argv_grp;
}

sub _CheckPointIs_cnv {
	my $sub_argv_grp = shift;
	$sub_argv_grp =~ s/^cnv\.//;
	switch ( $Ckp_Number{"cnv.$sub_argv_grp"}{$Method_CNV} ){
		case 1 {
			return unless $Ckp_Number{"bqsr.$sub_argv_grp"}{$Method_MAP} == 5;
			my $check_str = "call cnv tree $sub_argv_grp";
			my $qsub_bash = "step_2-3-1_call_cnv_$Method_CNV.$sub_argv_grp\_tree.sh";
			&_Check_To_Run ($check_str, $qsub_bash);
		}
		case 2 {
			my $check_str = "call cnv histogram $sub_argv_grp";
			my $qsub_bash = "step_2-3-2_call_cnv_$Method_CNV.$sub_argv_grp\_histogram.sh";
			&_Check_To_Run ($check_str, $qsub_bash);
		}
		case 3 {
			my $check_str = "call cnv stat $sub_argv_grp";
			my $qsub_bash = "step_2-3-3_call_cnv_$Method_CNV.$sub_argv_grp\_stat.sh";
			&_Check_To_Run ($check_str, $qsub_bash);
		}
		case 4 {
			my $check_str = "call cnv partition $sub_argv_grp";
			my $qsub_bash = "step_2-3-4_call_cnv_$Method_CNV.$sub_argv_grp\_partition.sh";
			&_Check_To_Run ($check_str, $qsub_bash);
		}
		case 5 {
			my $check_str = "call cnv call $sub_argv_grp";
			my $qsub_bash = "step_2-3-5_call_cnv_$Method_CNV.$sub_argv_grp\_call.sh";
			&_Check_To_Run ($check_str, $qsub_bash);
		}
	}
	undef $sub_argv_grp;
}

sub _CheckPointIs_sv {
	my $sub_argv_idv = shift;
	switch ( $Ckp_Number{sv}{$sub_argv_idv} ){
		case 1 {
		return unless $Ckp_Number{map}{$Method_MAP} == 3;
			foreach my $type ('DEL', 'DUP', 'INV', 'TRA') {				
				my $check_str = "call sv $type";
				my $qsub_bash = "step_2-4-1_call_sv_$Method_SV\_$type.sh";
				&_Check_To_Run ($check_str, $qsub_bash);
			}
		}
	}
	undef $sub_argv_idv;
}

sub _Check_To_Run {
	my $check_info = shift;
	my $qsub_shell = shift;
	my $run_marker = 0;
	my $log_file = "$Work_Dir/process/analysis.log";
	
	open FILE_HANDLE_LOG, $log_file or die $!;
	while (<FILE_HANDLE_LOG>) { $run_marker = 1 if (/^$check_info/) }
	close FILE_HANDLE_LOG;
	return 1 unless $run_marker == 0;
	
	while ( -e "$log_file.tmp" ){
		sleep 1;
	}
	
	system "while [ -e $log_file.tmp ]; do sleep 1; done; perl $Bin/UpdateLog.pl $log_file \"$check_info\" 1 > $log_file.tmp; mv $log_file.tmp $log_file";
	chdir "$Work_Dir/process/shells/";
	system "perl $Bin/MultiQsuber.pl $qsub_shell >> $Work_Dir/process/qsub.log; sleep 1";
	return 0;
}