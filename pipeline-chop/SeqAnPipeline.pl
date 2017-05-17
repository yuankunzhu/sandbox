#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Basename;
use File::Path;
use FindBin '$Bin';
use Cwd 'abs_path';
use lib $Bin;
use Switch;

#------------------------------------------------------------------
# get/check with parameters
#------------------------------------------------------------------
my ( $MANIFEST, $CONFIG, $CANCERPAIR );
my $OUTDIR = './sap';
my $UPDATE = 0;

my $REMOVES = 1;
my $Rm_Mark = "rm";
GetOptions(
	"m=s" => \$MANIFEST,
	"c=s" => \$CONFIG,
	"o=s" => \$OUTDIR,
	"p=s" => \$CANCERPAIR,
	"u=i" => \$UPDATE,
	"d=i" => \$REMOVES
);
if ( ! defined $MANIFEST ) {
	print "usage: $0 -m <manifest> -c <config> -o <output>\n";
	print "-m <manifest> file is required\n";
	exit 1;
}
if ( ! defined $CONFIG ) {
	print "usage: $0 -m <manifest> -c <config> -o <output>\n";
	print "-c <config> file is required\n";
	exit 1;
}
$MANIFEST = abs_path ( $MANIFEST );
$CONFIG = abs_path ( $CONFIG );
$OUTDIR = abs_path ( $OUTDIR );

# remove mediate file or not
$Rm_Mark = "# rm" if $REMOVES == 0;

#------------------------------------------------------------------
# read from config
#------------------------------------------------------------------
my ( $Method_REF, $Method_MAP, $Method_SNV, $Method_IDV );
my ( $Method_CNV, $Method_SV, $Method_ANN, $Method_RUN, $Method_VQRS);
my ( %ReadCfg, %Group_Content, @RunChr_Group, $Species );
open FILE_HANDLE_CONFIG, $CONFIG or die "[ERROR] from $CONFIG: $!";
while ( <FILE_HANDLE_CONFIG> ) {
	next if ( /^\s*$/ or /^\s*\#/ );
	chomp;
	my @read_field = split/\=/;
	if ( $read_field[0] eq 'AnaysisMethod') {
		&_extractMethod ( $read_field[1] );
	} else {
		$ReadCfg{$read_field[0]} = "$read_field[1]";
	}
	switch ( $Method_REF ) {
		case 'GRCh38' { $Species = "human" }
		case 'GRCh37' { $Species = "human" }
		case 'GRCh36' { $Species = "human" }
		case 'hg38'   { $Species = "human" }
		case 'hg19'   { $Species = "human" }
		case 'hg18'   { $Species = "human" }
		case 'GRCm37' { $Species = "mouse" }
		case 'GRCm38' { $Species = "mouse" }
		case 'mm10'   { $Species = "mouse" }
		case 'mm9'    { $Species = "mouse" }
		else		  { die "bad species/ref setting!"}
	}
	
	undef @read_field;
} close FILE_HANDLE_CONFIG;
&_defineSplicChrGroup;


#------------------------------------------------------------------
# read from fq manifest
#------------------------------------------------------------------
my ( @SampleID, %SampleHash, %Read1, %Read2, %ReadLen, %LaneInLib, %LibNumber, %InsertSize );
open FILE_HANDLE_MANIFEST, $MANIFEST or die "[ERROR] from $MANIFEST: $!";
while ( <FILE_HANDLE_MANIFEST> ) {
	my @read_field = split;
	my $cohort_id = $read_field[0];
	my $sample_id = $read_field[1];
	my $lib_id = $read_field[2];
	my $lane_id = $read_field[5];
	
	#$Cohort{$sample_id} = $cohort_id;
	$SampleHash{$sample_id} = 1;
	if($read_field[6]==1){
		$Read1{$lane_id} = "$read_field[7]/$sample_id/$lane_id\_$lib_id\_$read_field[6].fq.gz";
	}else{
		$Read2{$lane_id} = "$read_field[7]/$sample_id/$lane_id\_$lib_id\_$read_field[6].fq.gz";
	}
	
	#$ReadLen{$lane_id} = $read_field[3];
	$ReadLen{$sample_id} = $read_field[3];
	$InsertSize{$lane_id} = $read_field[4];

	$LaneInLib{$lib_id}{$lane_id} = 1;
	$LibNumber{$sample_id}{$lib_id} += 1;
	
	undef @read_field;
	undef $sample_id;
	undef $lane_id;
} close FILE_HANDLE_MANIFEST;
map { push @SampleID, $_ } ( keys %SampleHash );
undef %SampleHash;


#------------------------------------------------------------------
# read from cancer pair list
#------------------------------------------------------------------
my ( %NormalSample, %TumorSample, %PairLink, %PairMarker );
my $CancerMode = 0;
if ( defined $CANCERPAIR ) {
	$CancerMode = 1;
	open FILE_HANDLE_CANCERPAIR, $CANCERPAIR or die $!;
	while ( <FILE_HANDLE_CANCERPAIR> ) {
		chomp;
		my @read_field = split;
		die "bad cancer pair design, some sample repeat" if ( defined $PairMarker{$read_field[0]} or defined $PairMarker{$read_field[1]} );
		$PairLink{$read_field[0]} = $read_field[1];
		$PairLink{$read_field[1]} = $read_field[0];
		$PairMarker{$read_field[0]} = "$read_field[0]\_vs_$read_field[1]";
		$PairMarker{$read_field[1]} = "$read_field[0]\_vs_$read_field[1]";
		$NormalSample{$PairMarker{$read_field[0]}} = $read_field[0];
		$TumorSample{$PairMarker{$read_field[1]}} = $read_field[1];
		undef @read_field
	} close FILE_HANDLE_CANCERPAIR;
}

#------------------------------------------------------------------
# pre-file reading dong, now start real work
#------------------------------------------------------------------
my ( %Process_Path, %Results_Path );
my ( %Analysis_Log, %Analysis_Ckp );
my ( %Check_newRun, %Ckp_Type );
my $Now_Time = time;

#------------------------------------------------------------------
# define key structure 
#------------------------------------------------------------------
map { &_defineWorkSpace($_) } (@SampleID);

foreach my $sample_id ( @SampleID ) {
	my $check_newRun;
	#--------------------------------------------------------------
	# check output/job dir status
	#--------------------------------------------------------------
	$check_newRun = &_checkNewRun ( $sample_id );
	#	&_mapScripts ( $sample_id ) if defined $Ckp_Type{map};
	#	&_snvScripts ( $sample_id ) if defined $Ckp_Type{snv};
	#	&_idvScripts ( $sample_id ) if defined $Ckp_Type{idv};
	#	&_cnvScripts ( $sample_id ) if defined $Ckp_Type{cnv};
	#	&_svScripts  ( $sample_id ) if defined $Ckp_Type{sv};
	#--------------------------------------------------------------
	# create file structure for new run
	#--------------------------------------------------------------
	if ( $check_newRun == 1 ){
		&_creatWorkSpace ( $sample_id ) if $UPDATE == 0;
		&_mapScripts ( $sample_id ) if defined $Ckp_Type{map};
		&_snvScripts ( $sample_id ) if defined $Ckp_Type{snv};
		&_idvScripts ( $sample_id ) if defined $Ckp_Type{idv};
		&_cnvScripts ( $sample_id ) if defined $Ckp_Type{cnv};
		&_svScripts  ( $sample_id ) if defined $Ckp_Type{sv};
	} else { 
		$Now_Time = time;
		system "cp $Analysis_Log{$sample_id} $Analysis_Log{$sample_id}.$Now_Time";
	}
	system "perl $Bin/UpdateLog.pl $Analysis_Log{$sample_id} $sample_id 4 > $Analysis_Log{$sample_id}.tmp && mv $Analysis_Log{$sample_id}.tmp $Analysis_Log{$sample_id}";
	
	#--------------------------------------------------------------
	# qsub un-run jobs
	#--------------------------------------------------------------
	system "perl $Bin/JobProcessor.pl $OUTDIR/$sample_id $Analysis_Ckp{$sample_id} run $check_newRun";

	undef $sample_id;
}


#------------------------------------------------------------------
# generate run command shell 
#------------------------------------------------------------------
if ( $UPDATE == 0 ) {
	open FILE_HANDEL_RUN, ">$OUTDIR/command.sh" or die $!;
	mkpath "$OUTDIR/bin";
	system "cp $MANIFEST $CONFIG $OUTDIR/bin";
	$MANIFEST = basename ($MANIFEST);
	$CONFIG   = basename ($CONFIG);
	print FILE_HANDEL_RUN "perl $Bin/SeqAnPipeline.pl -c bin/$CONFIG -m bin/$MANIFEST -d $REMOVES -u $UPDATE -o $OUTDIR";
	if (defined $CANCERPAIR) {
		system "cp $CANCERPAIR $OUTDIR/bin";
		$CANCERPAIR = basename ($CANCERPAIR);
		print FILE_HANDEL_RUN " -p bin/$CANCERPAIR\n";
	}else {
		print FILE_HANDEL_RUN "\n";
	}close FILE_HANDEL_RUN;
}


#------------------------------------------------------------------
# subroutines
#------------------------------------------------------------------
sub _extractMethod {
	my $read_field = shift;
	my @read_field = split/\s+/, $read_field;
	foreach my $tmp (0 .. $#read_field){
		next if $tmp%2 == 1; 
		switch ( $read_field[$tmp] ) {
			case '-ref' { $Method_REF = $read_field[$tmp+1] }
			case '-map' { $Method_MAP = $read_field[$tmp+1]; $Ckp_Type{map} = 1 }
			case '-snv' { $Method_SNV = $read_field[$tmp+1]; $Ckp_Type{snv} = 1 }
			case '-idv' { $Method_IDV = $read_field[$tmp+1]; $Ckp_Type{idv} = 1 }
			case '-cnv' { $Method_CNV = $read_field[$tmp+1]; $Ckp_Type{cnv} = 1 }
			case '-sv'  { $Method_SV  = $read_field[$tmp+1]; $Ckp_Type{sv}  = 1 }
			case '-ann' { $Method_ANN = $read_field[$tmp+1] }
			case '-vqrs' { $Method_VQRS = $read_field[$tmp+1] }
			case '-runAs' { $Method_RUN = $read_field[$tmp+1] } 
			
		}
	}
	undef @read_field;
	undef $read_field;
}

sub _defineSplicChrGroup {
	switch ( $Method_REF ) {
		case 'GRCh37' {
			foreach (1 .. 18) {
				my $tmp = sprintf"G%02d",$_;
				push @RunChr_Group, $tmp;
				$Group_Content{$tmp} = $_ if $_ < 14;
			}
			$Group_Content{G14} = ("14 15");
			$Group_Content{G15} = ("16 17");
			$Group_Content{G16} = ("18 19");
			$Group_Content{G17} = ("20 21 22");
			$Group_Content{G18} = ("X Y MT GL000191.1 GL000192.1 GL000193.1 GL000194.1 GL000195.1 GL000196.1 GL000197.1 GL000198.1 GL000199.1 GL000200.1 GL000201.1 GL000202.1 GL000203.1 GL000204.1 GL000205.1 GL000206.1 GL000207.1 GL000208.1 GL000209.1 GL000210.1 GL000211.1 GL000212.1 GL000213.1 GL000214.1 GL000215.1 GL000216.1 GL000217.1 GL000218.1 GL000219.1 GL000220.1 GL000221.1 GL000222.1 GL000223.1 GL000224.1 GL000225.1 GL000226.1 GL000227.1 GL000228.1 GL000229.1 GL000230.1 GL000231.1 GL000232.1 GL000233.1 GL000234.1 GL000235.1 GL000236.1 GL000237.1 GL000238.1 GL000239.1 GL000240.1 GL000241.1 GL000242.1 GL000243.1 GL000244.1 GL000245.1 GL000246.1 GL000247.1 GL000248.1 GL000249.1");
		}
		case 'GRCm38' {
			foreach (1 .. 18) {
				my $tmp = sprintf"G%02d",$_;
				push @RunChr_Group, $tmp;
			}
			$Group_Content{G01} = ("1");
			$Group_Content{G02} = ("2");
			$Group_Content{G03} = ("X");
			$Group_Content{G04} = ("3");
			$Group_Content{G05} = ("4");
			$Group_Content{G06} = ("5");
			$Group_Content{G07} = ("6");
			$Group_Content{G08} = ("7");
			$Group_Content{G09} = ("10");
			$Group_Content{G10} = ("8");
			$Group_Content{G11} = ("14");
			$Group_Content{G12} = ("9");
			$Group_Content{G13} = ("11");
			$Group_Content{G14} = ("13");
			$Group_Content{G15} = ("12");
			$Group_Content{G17} = ("16 17");
			$Group_Content{G18} = ("Y 18 19");
			$Group_Content{G16} = ("15 JH584299.1 GL456233.1 JH584301.1 GL456211.1 GL456350.1 JH584293.1 GL456221.1 JH584297.1 JH584296.1 GL456354.1 JH584294.1 JH584298.1 JH584300.1 GL456219.1 GL456210.1 JH584303.1 JH584302.1 GL456212.1 JH584304.1 GL456379.1 GL456216.1 GL456393.1 GL456366.1 GL456367.1 GL456239.1 GL456213.1 GL456383.1 GL456385.1 GL456360.1 GL456378.1 GL456389.1 GL456372.1 GL456370.1 GL456381.1 GL456387.1 GL456390.1 GL456394.1 GL456392.1 GL456382.1 GL456359.1 GL456396.1 GL456368.1 MT JH584292.1 JH584295.1");
		}
		case 'mm10' {
			foreach (1 .. 17) {
				my $tmp = sprintf"G%02d",$_;
				push @RunChr_Group, $tmp;
			}
			$Group_Content{G01} = ("chr1");
			$Group_Content{G02} = ("chr2");
			$Group_Content{G03} = ("chr3");
			$Group_Content{G04} = ("chr4");
			$Group_Content{G05} = ("chr5");
			$Group_Content{G06} = ("chr6");
			$Group_Content{G07} = ("chr7");
			$Group_Content{G08} = ("chr8");
			$Group_Content{G09} = ("chr9");
			$Group_Content{G10} = ("chr10");
			$Group_Content{G11} = ("chr11");
			$Group_Content{G12} = ("chr12");
			$Group_Content{G13} = ("chr13");
			$Group_Content{G14} = ("chr14");
			$Group_Content{G15} = ("chr15 chr16");
			$Group_Content{G16} = ("chr17 chr18 chr19");
			$Group_Content{G17} = ("chrX chrY");
		}
	}
}

sub _defineWorkSpace {
	my $sub_argv = shift;
	$Process_Path{$sub_argv} = "$OUTDIR/$sub_argv/process";
	$Results_Path{$sub_argv} = "$OUTDIR/$sub_argv/results";
	$Analysis_Log{$sub_argv} = "$OUTDIR/$sub_argv/process/analysis.log";
	$Analysis_Ckp{$sub_argv} = "$OUTDIR/$sub_argv/process/analysis.ckp";
	if ( $CancerMode == 1 ) {
		$Process_Path{$PairMarker{$sub_argv}} = "$OUTDIR/cancer.pair-$PairMarker{$sub_argv}/process";
		$Analysis_Ckp{$PairMarker{$sub_argv}} = "$OUTDIR/cancer.pair-$PairMarker{$sub_argv}/analysis.log";
		$Analysis_Log{$PairMarker{$sub_argv}} = "$OUTDIR/cancer.pair-$PairMarker{$sub_argv}/analysis.ckp";
	}
	undef $sub_argv;
}

sub _checkNewRun {
	my $sub_argv = shift;
	if ( -e $Process_Path{$sub_argv} ) {
		if ( -e $Results_Path{$sub_argv} ) {
			if ( -e $Analysis_Ckp{$sub_argv} ) {
				print "[STATUS] sample $sub_argv old run found, continue from last breakpoint, ^_^\n";
				return 0;
			} else {
				print "[ERROR] sample $sub_argv missing .ckp file, stop now. T_T\n";
				exit 1;
			}
		} else {
			print "[ERROR] sample $sub_argv missing results folder, stop now. T_T\n";
			exit 1;
		}
	} elsif ( -e $Results_Path{$sub_argv} ) {
		print "[ERROR] sample $sub_argv missing process folder, stop now. T_T\n";
		exit 1;
	} else {
		print "[STATUS] sample $sub_argv new run will start, ^_^\n";
		return 1;
	}
	undef $sub_argv;
}

sub _creatWorkSpace {
	my $sub_argv = shift;
	my $lib_info;
	my $lane_info;
	
	mkpath "$Process_Path{$sub_argv}/runs/step_1-1_map";
	mkpath "$Process_Path{$sub_argv}/runs/step_2-1_snv";
	mkpath "$Process_Path{$sub_argv}/runs/step_2-2_idv";
	mkpath "$Process_Path{$sub_argv}/runs/step_2-3_cnv";
	mkpath "$Process_Path{$sub_argv}/runs/step_2-4_sv";
	mkpath "$Process_Path{$sub_argv}/shells";
	mkpath "$Results_Path{$sub_argv}";
	`touch $Analysis_Ckp{$sub_argv}`;
	`touch $Analysis_Log{$sub_argv}`;
	 
	 if ( $CancerMode == 1 && !(-e $Process_Path{$PairMarker{$sub_argv}}) ) {
	 	mkpath "$Process_Path{$PairMarker{$sub_argv}}/step_2-1_snv";
		mkpath "$Process_Path{$PairMarker{$sub_argv}}/step_2-2_idv";
		mkpath "$Process_Path{$PairMarker{$sub_argv}}/step_2-3_cnv";
		mkpath "$Process_Path{$PairMarker{$sub_argv}}/step_2-4_sv";
		`touch $Analysis_Ckp{$PairMarker{$sub_argv}} $Analysis_Log{$PairMarker{$sub_argv}}`;
	 }
	
	#--------------------------------------------------------------
	# initial fq read group info && checkpoint file for sample
	#--------------------------------------------------------------
	open FILE_HANDLE_SAMPLE_RGP, ">$Process_Path{$sub_argv}/runs/step_1-1_map/$sub_argv.rgp" or die $!;
	foreach my $tmp1 ( keys %{$LibNumber{$sub_argv}} ) {
		open FILE_HANDLE_LIB_RGP, ">$Process_Path{$sub_argv}/runs/step_1-1_map/$tmp1.rgp" or die $!;
		$lib_info .= "lib\t1\t$tmp1\n";
		foreach my $tmp2 ( keys %{$LaneInLib{$tmp1}} ) {
			$lane_info .= "lane\t1\t$tmp2\_1\nlane\t1\t$tmp2\_2\n";
			my $rgp_info = "\@RG\tID:$tmp2\tSM:$sub_argv\tLB:$tmp1\tPL:Illumina\n";
			print FILE_HANDLE_LIB_RGP $rgp_info;
			print FILE_HANDLE_SAMPLE_RGP $rgp_info;
		}
		close FILE_HANDLE_LIB_RGP;
	}
	close FILE_HANDLE_SAMPLE_RGP;
	my $initial_ckp = $lib_info.$lane_info;
	foreach my $tmp ('map', 'snv', 'idv', 'cnv', 'sv') {
		switch ($tmp) {
			case 'map' { next unless defined $Ckp_Type{$tmp}; $initial_ckp .= "map\t1\t$Method_MAP\n" }
			case 'snv' { next unless defined $Ckp_Type{$tmp}; $initial_ckp .= "snv\t1\t$Method_SNV\n" }
			case 'idv' { next unless defined $Ckp_Type{$tmp}; $initial_ckp .= "idv\t1\t$Method_IDV\n" }
			case 'cnv' { next unless defined $Ckp_Type{$tmp}; $initial_ckp .= "cnv\t1\t$Method_CNV\n" }
			case 'sv'  { next unless defined $Ckp_Type{$tmp}; $initial_ckp .= "sv\t1\t$Method_SV\n" }
		}
	}
	
	#--------------------------------------------------------------
	# initial checkpoint and log file for chr groups
	#--------------------------------------------------------------
	foreach my $tmp ('map', 'snv', 'idv', 'cnv', 'sv') {
		foreach my $chr_group ( @RunChr_Group ) {
			switch ($tmp) {
				case 'map' { next unless defined $Ckp_Type{$tmp}; $initial_ckp .= "bqsr.$chr_group\t1\t$Method_MAP\n" }
				case 'snv' { next unless defined $Ckp_Type{$tmp}; $initial_ckp .= "snv.$chr_group\t1\t$Method_SNV\n" }
				case 'idv' { next unless defined $Ckp_Type{$tmp}; $initial_ckp .= "idv.$chr_group\t1\t$Method_IDV\n" }
				case 'cnv' { next unless defined $Ckp_Type{$tmp}; $initial_ckp .= "cnv.$chr_group\t1\t$Method_CNV\n" }
			}
		}
	}
	open FILE_HANDLE_CKP, ">$Analysis_Ckp{$sub_argv}" or die "[ERROR] from $Analysis_Ckp{$sub_argv}: $!";
	print FILE_HANDLE_CKP $initial_ckp;
	close FILE_HANDLE_CKP;
	
	#system "egrep '^snv|^idv|^cnv|^sv' $Analysis_Ckp{$sub_argv} > $Analysis_Ckp{$PairMarker{$sub_argv}}";

	undef $initial_ckp;
	undef $sub_argv;
	undef $lib_info;
	undef $lane_info;
	
}

sub _mapScripts {
	my $sub_argv = shift;

	#--------------------------------------------------------------
	# initial mapping scripts
	#--------------------------------------------------------------
	my $all_lib_bam;
	my $lib_count = 0;
	foreach my $lib_id ( keys %{$LibNumber{$sub_argv}} ) {
		$lib_count += 1;
		my $all_lane_bam;
		my $lane_count = 0;
		foreach my $lane_id ( keys %{$LaneInLib{$lib_id}} ) {
			$lane_count += 1;
			switch ( $Method_MAP ) {
				case ('BWA') {
					#----------------------------------------------
					# bwa fq2sai read_1
					#----------------------------------------------
					open FILE_HANDLE_LIB_1, ">$Process_Path{$sub_argv}/shells/step_1-1_fq2sai_$lane_id\_1.sh" or die $!;
					my $initial_map_lib1 = <<"SCRIPTS";
#!/bin/sh
## bwa fq2sai for $lane_id read_1
#\$ -cwd -l virtual_free=45,h_vmem=75G -pe smp 36

{ $ReadCfg{BWA} aln $ReadCfg{BwaAlnOpt} $ReadCfg{$Method_REF} $Read1{$lane_id} > $Process_Path{$sub_argv}/runs/step_1-1_map/$lane_id\_1.sai && echo "$lane_id read1 fq2sai OK"; } || { echo "$lane_id read1 fq2sai failed"; while [ -e $Analysis_Log{$sub_argv}.tmp ]; do sleep 1; done; perl $Bin/UpdateLog.pl $Analysis_Log{$sub_argv} "fq2sai $lane_id\_1" 3 > $Analysis_Log{$sub_argv}.tmp && mv $Analysis_Log{$sub_argv}.tmp $Analysis_Log{$sub_argv}; exit 1; }

while [ -e $Analysis_Log{$sub_argv}.tmp ]; do sleep 1; done; perl $Bin/UpdateLog.pl $Analysis_Log{$sub_argv} "fq2sai $lane_id\_1" 0 > $Analysis_Log{$sub_argv}.tmp && mv $Analysis_Log{$sub_argv}.tmp $Analysis_Log{$sub_argv}

while [ -e $Analysis_Ckp{$sub_argv}.tmp ]; do sleep 1; done; perl $Bin/UpdateCkp.pl $Analysis_Ckp{$sub_argv} lane $lane_id\_1 2 > $Analysis_Ckp{$sub_argv}.tmp && mv $Analysis_Ckp{$sub_argv}.tmp $Analysis_Ckp{$sub_argv}

perl $Bin/JobProcessor.pl $OUTDIR/$sub_argv $Analysis_Ckp{$sub_argv} lane $lane_id

SCRIPTS
					print FILE_HANDLE_LIB_1 $initial_map_lib1;
					close FILE_HANDLE_LIB_1;
					undef $initial_map_lib1; 


					#----------------------------------------------
					# bwa fq2sai read_2
					#----------------------------------------------
					open FILE_HANDLE_LIB_2, ">$Process_Path{$sub_argv}/shells/step_1-2_fq2sai_$lane_id\_2.sh" or die $!;
					my $initial_map_lib2 = <<"SCRIPTS";
#!/bin/sh
## bwa fq2sai for $lane_id read_2
#\$ -cwd -l virtual_free=45,h_vmem=75G -pe smp 36

{ $ReadCfg{BWA} aln $ReadCfg{BwaAlnOpt} $ReadCfg{$Method_REF} $Read2{$lane_id} > $Process_Path{$sub_argv}/runs/step_1-1_map/$lane_id\_2.sai && echo "$lane_id read2 fq2sai OK"; } || { echo "$lane_id read2 fq2sai failed"; while [ -e $Analysis_Log{$sub_argv}.tmp ]; do sleep 1; done; perl $Bin/UpdateLog.pl $Analysis_Log{$sub_argv} "fq2sai $lane_id\_2" 3 > $Analysis_Log{$sub_argv}.tmp && mv $Analysis_Log{$sub_argv}.tmp $Analysis_Log{$sub_argv}; exit 1; }

while [ -e $Analysis_Log{$sub_argv}.tmp ]; do sleep 1; done; perl $Bin/UpdateLog.pl $Analysis_Log{$sub_argv} "fq2sai $lane_id\_2" 0 > $Analysis_Log{$sub_argv}.tmp && mv $Analysis_Log{$sub_argv}.tmp $Analysis_Log{$sub_argv}

while [ -e $Analysis_Ckp{$sub_argv}.tmp ]; do sleep 1; done; perl $Bin/UpdateCkp.pl $Analysis_Ckp{$sub_argv} lane $lane_id\_2 2 > $Analysis_Ckp{$sub_argv}.tmp && mv $Analysis_Ckp{$sub_argv}.tmp $Analysis_Ckp{$sub_argv}

perl $Bin/JobProcessor.pl $OUTDIR/$sub_argv $Analysis_Ckp{$sub_argv} lane $lane_id
SCRIPTS
					print FILE_HANDLE_LIB_2 $initial_map_lib2;
					close FILE_HANDLE_LIB_2;
					undef $initial_map_lib2;


					#----------------------------------------------
					# bwa sai2sam
					#----------------------------------------------
					open FILE_HANDLE_SAI2SAM, ">$Process_Path{$sub_argv}/shells/step_1-3_sai2sam_$lane_id.sh" or die $!;
					my $initial_map_sai2sam = <<"SCRIPTS";
#!/bin/sh
## bwa sampe for $lane_id
#\$ -cwd -l virtual_free=10G,h_vmem=20G -pe smp 3

{ $ReadCfg{BWA} sampe $ReadCfg{BwaSamOpt} -r "\@RG\\tID:$lane_id\\tSM:$sub_argv\\tLB:$lib_id\\tPL:Illumina" $ReadCfg{$Method_REF} $Process_Path{$sub_argv}/runs/step_1-1_map/$lane_id\_1.sai $Process_Path{$sub_argv}/runs/step_1-1_map/$lane_id\_2.sai $Read1{$lane_id} $Read2{$lane_id} > $Process_Path{$sub_argv}/runs/step_1-1_map/$lane_id.sam && echo "$lane_id sai2sam OK"; } || { echo "$lane_id sai2sam failed"; while [ -e $Analysis_Log{$sub_argv}.tmp ]; do sleep 1; done; perl $Bin/UpdateLog.pl $Analysis_Log{$sub_argv} "sai2sam $lane_id" 3 > $Analysis_Log{$sub_argv}.tmp && mv $Analysis_Log{$sub_argv}.tmp $Analysis_Log{$sub_argv}; exit 1; }

while [ -e $Analysis_Log{$sub_argv}.tmp ]; do sleep 1; done; perl $Bin/UpdateLog.pl $Analysis_Log{$sub_argv} "sai2sam $lane_id" 0 > $Analysis_Log{$sub_argv}.tmp && mv $Analysis_Log{$sub_argv}.tmp $Analysis_Log{$sub_argv}
$Rm_Mark $Process_Path{$sub_argv}/runs/step_1-1_map/$lane_id\_1.sai $Process_Path{$sub_argv}/runs/step_1-1_map/$lane_id\_2.sai

while [ -e $Analysis_Ckp{$sub_argv}.tmp ]; do sleep 1; done; perl $Bin/UpdateCkp.pl $Analysis_Ckp{$sub_argv} lane $lane_id\_1 3 > $Analysis_Ckp{$sub_argv}.tmp && mv $Analysis_Ckp{$sub_argv}.tmp $Analysis_Ckp{$sub_argv}

while [ -e $Analysis_Ckp{$sub_argv}.tmp ]; do sleep 1; done; perl $Bin/UpdateCkp.pl $Analysis_Ckp{$sub_argv} lane $lane_id\_2 3 > $Analysis_Ckp{$sub_argv}.tmp && mv $Analysis_Ckp{$sub_argv}.tmp $Analysis_Ckp{$sub_argv}

while [ -e $Analysis_Log{$sub_argv}.tmp ]; do sleep 1; done; perl $Bin/UpdateLog.pl $Analysis_Log{$sub_argv} "sam2bam $lane_id" 1 > $Analysis_Log{$sub_argv}.tmp && mv $Analysis_Log{$sub_argv}.tmp $Analysis_Log{$sub_argv}

/cm/shared/apps/sge/univa/8.1.3/bin/lx-amd64/qsub $Process_Path{$sub_argv}/shells/step_1-4_sam2bam_$lane_id.sh
SCRIPTS
					print FILE_HANDLE_SAI2SAM $initial_map_sai2sam;
					close FILE_HANDLE_SAI2SAM;
					undef $initial_map_sai2sam;


					#----------------------------------------------
					# samtools sam2bam & sort
					#----------------------------------------------
					open FILE_HANDLE_SAM2BAM, ">$Process_Path{$sub_argv}/shells/step_1-4_sam2bam_$lane_id.sh" or die $!;
					my $initial_map_sam2bam = <<"SCRIPTS";
#!/bin/sh
## samtools sam2bam for $lane_id
#\$ -cwd -l virtual_free=2G,h_vmem=10G -pe smp 2

{ $ReadCfg{Samtools} view -Sb $Process_Path{$sub_argv}/runs/step_1-1_map/$lane_id.sam -t $ReadCfg{$Method_REF}\.fai | $ReadCfg{Samtools} sort -m 2000000000 - $Process_Path{$sub_argv}/runs/step_1-1_map/$lane_id.sort && echo "$lane_id sam2bam OK"; } || { echo "$lane_id sam2bam failed"; while [ -e $Analysis_Log{$sub_argv}.tmp ]; do sleep 1; done; perl $Bin/UpdateLog.pl $Analysis_Log{$sub_argv} "sam2bam $lane_id" 3 > $Analysis_Log{$sub_argv}.tmp && mv $Analysis_Log{$sub_argv}.tmp $Analysis_Log{$sub_argv} exit 1; }

while [ -e $Analysis_Log{$sub_argv}.tmp ]; do sleep 1; done; perl $Bin/UpdateLog.pl $Analysis_Log{$sub_argv} "sam2bam $lane_id" 0 > $Analysis_Log{$sub_argv}.tmp && mv $Analysis_Log{$sub_argv}.tmp $Analysis_Log{$sub_argv}
$Rm_Mark $Process_Path{$sub_argv}/runs/step_1-1_map/$lane_id.sam

while [ -e $Analysis_Ckp{$sub_argv}.tmp ]; do sleep 1; done; perl $Bin/UpdateCkp.pl $Analysis_Ckp{$sub_argv} lane $lane_id\_1 4 > $Analysis_Ckp{$sub_argv}.tmp && mv $Analysis_Ckp{$sub_argv}.tmp $Analysis_Ckp{$sub_argv}

while [ -e $Analysis_Ckp{$sub_argv}.tmp ]; do sleep 1; done; perl $Bin/UpdateCkp.pl $Analysis_Ckp{$sub_argv} lane $lane_id\_2 4 > $Analysis_Ckp{$sub_argv}.tmp && mv $Analysis_Ckp{$sub_argv}.tmp $Analysis_Ckp{$sub_argv}

perl $Bin/JobProcessor.pl $OUTDIR/$sub_argv $Analysis_Ckp{$sub_argv} lib $lib_id
SCRIPTS
					print FILE_HANDLE_SAM2BAM $initial_map_sam2bam;
					close FILE_HANDLE_SAM2BAM;
					undef $initial_map_sam2bam;

					$all_lane_bam .= " $Process_Path{$sub_argv}/runs/step_1-1_map/$lane_id.sort.bam";
				}
			}
				#case ('Novoalign')
		}


		#----------------------------------------------------------
		# samtools merge & rmdup by lib
		#----------------------------------------------------------
		open FILE_HANDLE_MERGE, ">$Process_Path{$sub_argv}/shells/step_1-5_merge_$lib_id.sh" or die $!;
		my $initial_map_merge;
		## check bam file number, see merge is taken or not
		if ( $lane_count == 1 ) {
			$initial_map_merge = <<"SCRIPTS";
#!/bin/sh
## samtools merge/rmdup for $lib_id
#\$ -cwd -l virtual_free=2G,h_vmem=10G -pe smp 2

mv $all_lane_bam $Process_Path{$sub_argv}/runs/step_1-1_map/$lib_id.bam

$ReadCfg{Samtools} rmdup $Process_Path{$sub_argv}/runs/step_1-1_map/$lib_id.bam $Process_Path{$sub_argv}/runs/step_1-1_map/$lib_id.rmdup.bam
$Rm_Mark $Process_Path{$sub_argv}/runs/step_1-1_map/$lib_id.bam

while [ -e $Analysis_Log{$sub_argv}.tmp ]; do sleep 1; done; perl $Bin/UpdateLog.pl $Analysis_Log{$sub_argv} "mergeLib $lib_id" 0 > $Analysis_Log{$sub_argv}.tmp && mv $Analysis_Log{$sub_argv}.tmp $Analysis_Log{$sub_argv}

while [ -e $Analysis_Ckp{$sub_argv}.tmp ]; do sleep 1; done; perl $Bin/UpdateCkp.pl $Analysis_Ckp{$sub_argv} lib $lib_id 2 > $Analysis_Ckp{$sub_argv}.tmp && mv $Analysis_Ckp{$sub_argv}.tmp $Analysis_Ckp{$sub_argv}

perl $Bin/JobProcessor.pl $OUTDIR/$sub_argv $Analysis_Ckp{$sub_argv} map $Method_MAP
SCRIPTS

		} else {
			$initial_map_merge = <<"SCRIPTS";
#!/bin/sh
## samtools merge/rmdup for $lib_id
#\$ -cwd -l virtual_free=2G,h_vmem=10G -pe smp 2

{ $ReadCfg{Samtools} merge -h $Process_Path{$sub_argv}/runs/step_1-1_map/$lib_id.rgp $Process_Path{$sub_argv}/runs/step_1-1_map/$lib_id.bam $all_lane_bam && echo "mergeLib $lib_id OK"; } || { echo "mergeLib $lib_id failed"; while [ -e $Analysis_Log{$sub_argv}.tmp ]; do sleep 1; done; while [ -e $Analysis_Log{$sub_argv}.tmp ]; do sleep 1; done; perl $Bin/UpdateLog.pl $Analysis_Log{$sub_argv} "mergeLib $lib_id" 3 > $Analysis_Log{$sub_argv}.tmp && mv $Analysis_Log{$sub_argv}.tmp $Analysis_Log{$sub_argv}; exit 1; }
$Rm_Mark $all_lane_bam

$ReadCfg{Samtools} rmdup $Process_Path{$sub_argv}/runs/step_1-1_map/$lib_id.bam $Process_Path{$sub_argv}/runs/step_1-1_map/$lib_id.rmdup.bam
$Rm_Mark $Process_Path{$sub_argv}/runs/step_1-1_map/$lib_id.bam

while [ -e $Analysis_Log{$sub_argv}.tmp ]; do sleep 1; done; perl $Bin/UpdateLog.pl $Analysis_Log{$sub_argv} "mergeLib $lib_id" 0 > $Analysis_Log{$sub_argv}.tmp && mv $Analysis_Log{$sub_argv}.tmp $Analysis_Log{$sub_argv}

while [ -e $Analysis_Ckp{$sub_argv}.tmp ]; do sleep 1; done; perl $Bin/UpdateCkp.pl $Analysis_Ckp{$sub_argv} lib $lib_id 2 > $Analysis_Ckp{$sub_argv}.tmp && mv $Analysis_Ckp{$sub_argv}.tmp $Analysis_Ckp{$sub_argv}

perl $Bin/JobProcessor.pl $OUTDIR/$sub_argv $Analysis_Ckp{$sub_argv} map $Method_MAP
SCRIPTS

		}

		print FILE_HANDLE_MERGE $initial_map_merge;
		close FILE_HANDLE_MERGE;
		undef $initial_map_merge;
		undef $lane_count;

		$all_lib_bam .= " $Process_Path{$sub_argv}/runs/step_1-1_map/$lib_id.rmdup.bam";
	}


	#--------------------------------------------------------------
	# samtools merge by sample
	#--------------------------------------------------------------
	open FILE_HANDLE_MERGE, ">$Process_Path{$sub_argv}/shells/step_1-6_merge_$sub_argv.sh" or die $!;
	my $initial_map_merge;
	## check bam file number, see merge is taken or not
	if ( $lib_count == 1 ) {
		$initial_map_merge = <<"SCRIPTS";
#!/bin/sh
## samtools merge/rmdup for $sub_argv
#\$ -cwd -l virtual_free=2G,h_vmem=10G -pe smp 2

mv $all_lib_bam $Process_Path{$sub_argv}/runs/step_1-1_map/$sub_argv.raw.bam
$ReadCfg{Samtools} index $Process_Path{$sub_argv}/runs/step_1-1_map/$sub_argv.raw.bam

while [ -e $Analysis_Log{$sub_argv}.tmp ]; do sleep 1; done; perl $Bin/UpdateCkp.pl $Analysis_Ckp{$sub_argv} map $Method_MAP 2 > $Analysis_Ckp{$sub_argv}.tmp && mv $Analysis_Ckp{$sub_argv}.tmp $Analysis_Ckp{$sub_argv}

SCRIPTS
	} else {
		$initial_map_merge = <<"SCRIPTS";
#!/bin/sh
## samtools merge/rmdup for $sub_argv
#\$ -cwd -l virtual_free=2G,h_vmem=10G -pe smp 2

{ $ReadCfg{Samtools} merge -h $Process_Path{$sub_argv}/runs/step_1-1_map/$sub_argv.rgp $Process_Path{$sub_argv}/runs/step_1-1_map/$sub_argv.raw.bam $all_lib_bam && echo "mergeSample $Method_MAP OK"; $ReadCfg{Samtools} index $Process_Path{$sub_argv}/runs/step_1-1_map/$sub_argv.raw.bam; } || { echo "mergeSample $Method_MAP failed"; while [ -e $Analysis_Log{$sub_argv}.tmp ]; do sleep 1; done; perl $Bin/UpdateLog.pl $Analysis_Log{$sub_argv} "mergeSample $Method_MAP" 3 > $Analysis_Log{$sub_argv}.tmp && mv $Analysis_Log{$sub_argv}.tmp $Analysis_Log{$sub_argv}; exit 1; }
$Rm_Mark $all_lib_bam

while [ -e $Analysis_Log{$sub_argv}.tmp ]; do sleep 1; done; perl $Bin/UpdateLog.pl $Analysis_Log{$sub_argv} "mergeSample $sub_argv" 0 > $Analysis_Log{$sub_argv}.tmp && mv $Analysis_Log{$sub_argv}.tmp $Analysis_Log{$sub_argv}

while [ -e $Analysis_Ckp{$sub_argv}.tmp ]; do sleep 1; done; perl $Bin/UpdateCkp.pl $Analysis_Ckp{$sub_argv} map $Method_MAP 2 > $Analysis_Ckp{$sub_argv}.tmp && mv $Analysis_Ckp{$sub_argv}.tmp $Analysis_Ckp{$sub_argv}

SCRIPTS
	}
	
	map { $initial_map_merge .= "perl $Bin/JobProcessor.pl $OUTDIR/$sub_argv $Analysis_Ckp{$sub_argv} bqsr.$_ $Method_MAP\n"} (@RunChr_Group);
	
	print FILE_HANDLE_MERGE $initial_map_merge;
	close FILE_HANDLE_MERGE;
	undef $initial_map_merge;
	undef $lib_count;

	
	#--------------------------------------------------------------
	# BQSR for each chr group
	#--------------------------------------------------------------
	foreach my $runchr_group ( @RunChr_Group ) {
	
	
		#----------------------------------------------------------
		# gatk indel re-aligner RealignerTargetCreator
		#----------------------------------------------------------
		open FILE_HANDLE_REALIGNER, ">$Process_Path{$sub_argv}/shells/step_1-7-1_RealignerTargetCreator_$sub_argv.$runchr_group.sh" or die $!;
		my $known_indel;
		switch ( $Species ){
			case 'human' { $known_indel = "--known $ReadCfg{'1000g.indel'} --known $ReadCfg{'mills.indel'}" }
			case 'mouse' { $known_indel = "--known $ReadCfg{'sanger.indel'}" }
		}
		my $initial_map_realigner = <<"SCRIPTS";
#!/bin/sh
## IndelRealigner for $sub_argv
#\$ -cwd -l virtual_free=5G,h_vmem=10G -pe smp 3

$ReadCfg{Java} -Xmx4G -Djava.io.tmpdir=/mnt/isilon/cbttc/tester/javaTmp/$sub_argv/RealignerTargetCreator -jar $ReadCfg{GATK} \\
	-T RealignerTargetCreator \\
	-R $ReadCfg{$Method_REF} \\
	-I $Process_Path{$sub_argv}/runs/step_1-1_map/$sub_argv.raw.bam \\
	-o $Process_Path{$sub_argv}/runs/step_1-1_map/$sub_argv.$runchr_group.IndelRealigner.intervals \\
	-nt 16 \\
	-L $ReadCfg{"$Method_REF.root"}/newDefine/$runchr_group.intervals \\
	$known_indel || { echo 'RealignerTargetCreator failed'; while [ -e $Analysis_Log{$sub_argv}.tmp ]; do sleep 1; done; perl $Bin/UpdateLog.pl $Analysis_Log{$sub_argv} "RealignerTargetCreator $runchr_group" 3 > $Analysis_Log{$sub_argv}.tmp && mv $Analysis_Log{$sub_argv}.tmp $Analysis_Log{$sub_argv}; exit 1; }

while [ -e $Analysis_Log{$sub_argv}.tmp ]; do sleep 1; done; perl $Bin/UpdateLog.pl $Analysis_Log{$sub_argv} "RealignerTargetCreator $runchr_group" 0 > $Analysis_Log{$sub_argv}.tmp && mv $Analysis_Log{$sub_argv}.tmp $Analysis_Log{$sub_argv}

while [ -e $Analysis_Ckp{$sub_argv}.tmp ]; do sleep 1; done; perl $Bin/UpdateCkp.pl $Analysis_Ckp{$sub_argv} bqsr.$runchr_group $Method_MAP 2 > $Analysis_Ckp{$sub_argv}.tmp && mv $Analysis_Ckp{$sub_argv}.tmp $Analysis_Ckp{$sub_argv}

perl $Bin/JobProcessor.pl $OUTDIR/$sub_argv $Analysis_Ckp{$sub_argv} bqsr.$runchr_group $Method_MAP
SCRIPTS
		print FILE_HANDLE_REALIGNER $initial_map_realigner;
		close FILE_HANDLE_REALIGNER;
		undef $initial_map_realigner;
		

		#----------------------------------------------------------
		# gatk indel re-aligner IndelRealigner
		#----------------------------------------------------------
		open FILE_HANDLE_REALIGNER, ">$Process_Path{$sub_argv}/shells/step_1-7-2_IndelRealigner_$sub_argv.$runchr_group.sh" or die $!;
		switch ( $Species ){
			case 'human' { $known_indel = "-known $ReadCfg{'1000g.indel'} -known $ReadCfg{'mills.indel'}" }
			case 'mouse' { $known_indel = "-known $ReadCfg{'sanger.indel'}" }
		}
		$initial_map_realigner = <<"SCRIPTS";
#!/bin/sh
## IndelRealigner for $sub_argv
#\$ -cwd -l virtual_free=5G,h_vmem=10G -pe smp 3

$ReadCfg{Java} -Xmx4G -Djava.io.tmpdir=/mnt/isilon/cbttc/tester/javaTmp/$sub_argv/IndelRealigner -jar $ReadCfg{GATK} \\
	-T IndelRealigner \\
	-R $ReadCfg{$Method_REF} \\
	-I $Process_Path{$sub_argv}/runs/step_1-1_map/$sub_argv.raw.bam \\
	-targetIntervals $Process_Path{$sub_argv}/runs/step_1-1_map/$sub_argv.$runchr_group.IndelRealigner.intervals \\
	-o $Process_Path{$sub_argv}/runs/step_1-1_map/$sub_argv.$runchr_group.realigned.bam \\
	-L $ReadCfg{"$Method_REF.root"}/newDefine/$runchr_group.intervals \\
	$known_indel || { echo 'IndelRealigner failed'; while [ -e $Analysis_Log{$sub_argv}.tmp ]; do sleep 1; done; perl $Bin/UpdateLog.pl $Analysis_Log{$sub_argv} "IndelRealigner $runchr_group" 3 > $Analysis_Log{$sub_argv}.tmp && mv $Analysis_Log{$sub_argv}.tmp $Analysis_Log{$sub_argv}; exit 1; }

while [ -e $Analysis_Log{$sub_argv}.tmp ]; do sleep 1; done; perl $Bin/UpdateLog.pl $Analysis_Log{$sub_argv} "IndelRealigner $runchr_group" 0 > $Analysis_Log{$sub_argv}.tmp && mv $Analysis_Log{$sub_argv}.tmp $Analysis_Log{$sub_argv}

while [ -e $Analysis_Ckp{$sub_argv}.tmp ]; do sleep 1; done; perl $Bin/UpdateCkp.pl $Analysis_Ckp{$sub_argv} bqsr.$runchr_group $Method_MAP 3 > $Analysis_Ckp{$sub_argv}.tmp && mv $Analysis_Ckp{$sub_argv}.tmp $Analysis_Ckp{$sub_argv}

perl $Bin/JobProcessor.pl $OUTDIR/$sub_argv $Analysis_Ckp{$sub_argv} bqsr.$runchr_group $Method_MAP
SCRIPTS
		print FILE_HANDLE_REALIGNER $initial_map_realigner;
		close FILE_HANDLE_REALIGNER;
		undef $initial_map_realigner;
	
	
		#----------------------------------------------------------
		# gatk base recalibration BaseRecalibrator
		#----------------------------------------------------------
		open FILE_HANDLE_BQSR, ">$Process_Path{$sub_argv}/shells/step_1-7-3_BaseRecalibrator_$sub_argv.$runchr_group.sh" or die $!;
		switch ( $Species ){
			case 'human' { $known_indel = "--knownSites $ReadCfg{'1000g.indel'} --knownSites $ReadCfg{'mills.indel'}" }
			case 'mouse' { $known_indel = "--knownSites $ReadCfg{'sanger.indel'}" }
		}
		my $initial_map_bqsr = <<"SCRIPTS";
#!/bin/sh
## BaseRecalibrator for $sub_argv
#\$ -cwd -l virtual_free=5G,h_vmem=10G -pe smp 3

$ReadCfg{Java} -Xmx4G -Djava.io.tmpdir=/mnt/isilon/cbttc/tester/javaTmp/$sub_argv/BaseRecalibrator -jar $ReadCfg{GATK} \\
	-T BaseRecalibrator \\
	-R $ReadCfg{$Method_REF} \\
	-I $Process_Path{$sub_argv}/runs/step_1-1_map/$sub_argv.$runchr_group.realigned.bam \\
	-o $Process_Path{$sub_argv}/runs/step_1-1_map/$sub_argv.$runchr_group.bqsr.tables \\
	-L $ReadCfg{"$Method_REF.root"}/newDefine/$runchr_group.intervals \\
	-nct 8 \\
	$known_indel ||  { echo 'BaseRecalibrator failed'; while [ -e $Analysis_Log{$sub_argv}.tmp ]; do sleep 1; done; perl $Bin/UpdateLog.pl $Analysis_Log{$sub_argv} "BaseRecalibrator $runchr_group" 3 > $Analysis_Log{$sub_argv}.tmp && mv $Analysis_Log{$sub_argv}.tmp $Analysis_Log{$sub_argv}; exit 1; }

while [ -e $Analysis_Log{$sub_argv}.tmp ]; do sleep 1; done; perl $Bin/UpdateLog.pl $Analysis_Log{$sub_argv} "BaseRecalibrator $runchr_group" 0 > $Analysis_Log{$sub_argv}.tmp && mv $Analysis_Log{$sub_argv}.tmp $Analysis_Log{$sub_argv}

while [ -e $Analysis_Ckp{$sub_argv}.tmp ]; do sleep 1; done; perl $Bin/UpdateCkp.pl $Analysis_Ckp{$sub_argv} bqsr.$runchr_group $Method_MAP 4 > $Analysis_Ckp{$sub_argv}.tmp && mv $Analysis_Ckp{$sub_argv}.tmp $Analysis_Ckp{$sub_argv}

perl $Bin/JobProcessor.pl $OUTDIR/$sub_argv $Analysis_Ckp{$sub_argv} bqsr.$runchr_group $Method_MAP
SCRIPTS
		print FILE_HANDLE_BQSR $initial_map_bqsr;
		close FILE_HANDLE_BQSR;
		undef $initial_map_bqsr;
		
		
		#----------------------------------------------------------
		# gatk base recalibration PrintReads
		#----------------------------------------------------------
		open FILE_HANDLE_BQSR, ">$Process_Path{$sub_argv}/shells/step_1-7-4_PrintReads_$sub_argv.$runchr_group.sh" or die $!;
		$initial_map_bqsr = <<"SCRIPTS";
#!/bin/sh
## BaseRecalibrator for $sub_argv
#\$ -cwd -l virtual_free=5G,h_vmem=10G -pe smp 3
   
$ReadCfg{Java} -Xmx4G -Djava.io.tmpdir=/mnt/isilon/cbttc/tester/javaTmp/$sub_argv/PrintReads -jar $ReadCfg{GATK} \\
	-T PrintReads \\
	-R $ReadCfg{$Method_REF} \\
	-I $Process_Path{$sub_argv}/runs/step_1-1_map/$sub_argv.$runchr_group.realigned.bam \\
	-BQSR $Process_Path{$sub_argv}/runs/step_1-1_map/$sub_argv.$runchr_group.bqsr.tables \\
	-o $Process_Path{$sub_argv}/runs/step_1-1_map/$sub_argv.$runchr_group.bam \\
	-L $ReadCfg{"$Method_REF.root"}/newDefine/$runchr_group.intervals \\
	-nct 8 ||  { echo 'PrintReads failed'; while [ -e $Analysis_Log{$sub_argv}.tmp ]; do sleep 1; done; perl $Bin/UpdateLog.pl $Analysis_Log{$sub_argv} "PrintReads $runchr_group" 3 > $Analysis_Log{$sub_argv}.tmp && mv $Analysis_Log{$sub_argv}.tmp $Analysis_Log{$sub_argv}; exit 1; }

$Rm_Mark $Process_Path{$sub_argv}/runs/step_1-1_map/$sub_argv.$runchr_group.realigned.bam $Process_Path{$sub_argv}/runs/step_1-1_map/$sub_argv.$runchr_group.realigned.bai $Process_Path{$sub_argv}/runs/step_1-1_map/$sub_argv.$runchr_group.bqsr.tables $Process_Path{$sub_argv}/runs/step_1-1_map/$sub_argv.$runchr_group.IndelRealigner.intervals

while [ -e $Analysis_Log{$sub_argv}.tmp ]; do sleep 1; done; perl $Bin/UpdateLog.pl $Analysis_Log{$sub_argv} "PrintReads $runchr_group" 0 > $Analysis_Log{$sub_argv}.tmp && mv $Analysis_Log{$sub_argv}.tmp $Analysis_Log{$sub_argv}

while [ -e $Analysis_Ckp{$sub_argv}.tmp ]; do sleep 1; done; perl $Bin/UpdateCkp.pl $Analysis_Ckp{$sub_argv} bqsr.$runchr_group $Method_MAP 5 > $Analysis_Ckp{$sub_argv}.tmp && mv $Analysis_Ckp{$sub_argv}.tmp $Analysis_Ckp{$sub_argv}

# perl $Bin/CheckGroup.pl $Analysis_Ckp{$sub_argv} bqsr 5 && { while [ -e $Analysis_Log{$sub_argv}.tmp ]; do sleep 1; done; perl $Bin/UpdateLog.pl $Analysis_Log{$sub_argv} "merge recal bam $Method_MAP" 1 > $Analysis_Log{$sub_argv}.tmp; mv $Analysis_Log{$sub_argv}.tmp $Analysis_Log{$sub_argv}; /cm/shared/apps/sge/univa/8.1.3/bin/lx-amd64/qsub $Process_Path{$sub_argv}/shells/step_1-8_merge_recal_$sub_argv.sh; }

perl $Bin/CheckGroup.pl $Analysis_Ckp{$sub_argv} bqsr 5 && perl $Bin/JobProcessor.pl $OUTDIR/$sub_argv $Analysis_Ckp{$sub_argv} map $Method_MAP


SCRIPTS
		foreach ('snv', 'idv', 'cnv', 'sv'){
			next unless defined $Ckp_Type{$_};
			$initial_map_bqsr .= "perl $Bin/JobProcessor.pl $OUTDIR/$sub_argv $Analysis_Ckp{$sub_argv} snv.$runchr_group $Method_SNV\n" if /snv/;
			$initial_map_bqsr .= "perl $Bin/JobProcessor.pl $OUTDIR/$sub_argv $Analysis_Ckp{$sub_argv} idv.$runchr_group $Method_IDV\n" if /idv/;
			$initial_map_bqsr .= "perl $Bin/JobProcessor.pl $OUTDIR/$sub_argv $Analysis_Ckp{$sub_argv} cnv.$runchr_group $Method_CNV\n" if /cnv/;
		}
		print FILE_HANDLE_BQSR $initial_map_bqsr;
		close FILE_HANDLE_BQSR;
		undef $initial_map_bqsr;
	}

	#------------------------------------------------------------------
	# final merge the recal bams
	#------------------------------------------------------------------
	open FILE_HANDLE_MERGE, ">$Process_Path{$sub_argv}/shells/step_1-8_merge_recal_$sub_argv.sh" or die $!;
	
	# make all the recal bams into str
	my $join_bam = join( ".bam $Process_Path{$sub_argv}/runs/step_1-1_map/$sub_argv.", @RunChr_Group);
	$initial_map_merge = <<"SCRIPTS";
#!/bin/sh
## samtools merge/rmdup for $sub_argv
#\$ -cwd -l virtual_free=2G,h_vmem=10G -pe smp 2

$Rm_Mark $Process_Path{$sub_argv}/runs/step_1-1_map/$sub_argv.raw.bam $Process_Path{$sub_argv}/runs/step_1-1_map/$sub_argv.raw.bam.bai

{ $ReadCfg{Samtools} merge -h $Process_Path{$sub_argv}/runs/step_1-1_map/$sub_argv.rgp $Process_Path{$sub_argv}/runs/step_1-1_map/$sub_argv.bam $Process_Path{$sub_argv}/runs/step_1-1_map/$sub_argv.$join_bam.bam && echo "mergeSample $Method_MAP OK"; $ReadCfg{Samtools} index $Process_Path{$sub_argv}/runs/step_1-1_map/$sub_argv.bam; } || { echo "mergeSample $Method_MAP failed"; while [ -e $Analysis_Log{$sub_argv}.tmp ]; do sleep 1; done; perl $Bin/UpdateLog.pl $Analysis_Log{$sub_argv} "merge recal bam $Method_MAP" 3 > $Analysis_Log{$sub_argv}.tmp && mv $Analysis_Log{$sub_argv}.tmp $Analysis_Log{$sub_argv}; exit 1; }

$ReadCfg{Samtools} index $Process_Path{$sub_argv}/runs/step_1-1_map/$sub_argv.bam

while [ -e $Analysis_Log{$sub_argv}.tmp ]; do sleep 1; done; perl $Bin/UpdateLog.pl $Analysis_Log{$sub_argv} "merge recal bam $Method_MAP" 0 > $Analysis_Log{$sub_argv}.tmp && mv $Analysis_Log{$sub_argv}.tmp $Analysis_Log{$sub_argv}

while [ -e $Analysis_Ckp{$sub_argv}.tmp ]; do sleep 1; done; perl $Bin/UpdateCkp.pl $Analysis_Ckp{$sub_argv} map $Method_MAP 3 > $Analysis_Ckp{$sub_argv}.tmp && mv $Analysis_Ckp{$sub_argv}.tmp $Analysis_Ckp{$sub_argv}



SCRIPTS
	$initial_map_merge .= "perl $Bin/JobProcessor.pl $OUTDIR/$sub_argv $Analysis_Ckp{$sub_argv} sv $Method_SV\n" if defined $Ckp_Type{sv};
	print FILE_HANDLE_MERGE $initial_map_merge;
	close FILE_HANDLE_MERGE;
	undef $initial_map_merge;
	
}

sub _snvScripts {
	my $sub_argv = shift;
	foreach my $runchr_group ( @RunChr_Group ) {
		#--------------------------------------------------------------
		# initial snv scripts
		#--------------------------------------------------------------
		switch ( $Method_SNV ) {
			case ('GATK') {
				#------------------------------------------------------
				# gatk HaplotypeCaller
				#------------------------------------------------------
				open FILE_HANDLE_CALL, ">$Process_Path{$sub_argv}/shells/step_2-1-1_call_snv_$Method_SNV.$runchr_group.sh" or die $!;
				my $initial_call = <<"SCRIPTS";
#!/bin/sh
## BaseRecalibrator for $sub_argv
#\$ -cwd -l virtual_free=30G,h_vmem=50G -pe smp 8

$ReadCfg{Java} -Xmx10G -Xms8G -Djava.io.tmpdir=/mnt/isilon/cbttc/tester/javaTmp/$sub_argv/HaplotypeCaller -jar $ReadCfg{GATK} \\
	-T HaplotypeCaller \\
	-R $ReadCfg{$Method_REF} \\
	-I $Process_Path{$sub_argv}/runs/step_1-1_map/$sub_argv.$runchr_group.bam \\
	-o $Process_Path{$sub_argv}/runs/step_2-1_snv/$sub_argv.$runchr_group.HaplotypeCaller.vcf \\
	--dbsnp $ReadCfg{"snp138.$Method_REF"} \\
	-L $ReadCfg{"$Method_REF.root"}/newDefine/$runchr_group.intervals \\
	-stand_call_conf 50.0 \\
	-stand_emit_conf 10.0 ||  { while [ -e $Analysis_Log{$sub_argv}.tmp ]; do sleep 1; done; perl $Bin/UpdateLog.pl $Analysis_Log{$sub_argv} "call snv $runchr_group" 3 > $Analysis_Log{$sub_argv}.tmp && mv $Analysis_Log{$sub_argv}.tmp $Analysis_Log{$sub_argv}; while [ -e $Analysis_Log{$sub_argv}.tmp ]; do sleep 1; done; perl $Bin/UpdateLog.pl $Analysis_Log{$sub_argv} "call idv $runchr_group" 3 > $Analysis_Log{$sub_argv}.tmp && mv $Analysis_Log{$sub_argv}.tmp $Analysis_Log{$sub_argv}; exit 1; }
   
while [ -e $Analysis_Log{$sub_argv}.tmp ]; do sleep 1; done; perl $Bin/UpdateLog.pl $Analysis_Log{$sub_argv} "call snv $runchr_group" 0 > $Analysis_Log{$sub_argv}.tmp && mv $Analysis_Log{$sub_argv}.tmp $Analysis_Log{$sub_argv}

while [ -e $Analysis_Ckp{$sub_argv}.tmp ]; do sleep 1; done; perl $Bin/UpdateCkp.pl $Analysis_Ckp{$sub_argv} snv.$runchr_group $Method_SNV 2 > $Analysis_Ckp{$sub_argv}.tmp && mv $Analysis_Ckp{$sub_argv}.tmp $Analysis_Ckp{$sub_argv}

while [ -e $Analysis_Log{$sub_argv}.tmp ]; do sleep 1; done; perl $Bin/UpdateLog.pl $Analysis_Log{$sub_argv} "call idv $runchr_group" 0 > $Analysis_Log{$sub_argv}.tmp && mv $Analysis_Log{$sub_argv}.tmp $Analysis_Log{$sub_argv}

while [ -e $Analysis_Ckp{$sub_argv}.tmp ]; do sleep 1; done; perl $Bin/UpdateCkp.pl $Analysis_Ckp{$sub_argv} idv.$runchr_group $Method_IDV 2 > $Analysis_Ckp{$sub_argv}.tmp && mv $Analysis_Ckp{$sub_argv}.tmp $Analysis_Ckp{$sub_argv}

perl $Bin/JobProcessor.pl $OUTDIR/$sub_argv $Analysis_Ckp{$sub_argv} snv.$runchr_group $Method_SNV
perl $Bin/JobProcessor.pl $OUTDIR/$sub_argv $Analysis_Ckp{$sub_argv} idv.$runchr_group $Method_SNV
SCRIPTS
				print FILE_HANDLE_CALL $initial_call;
				close FILE_HANDLE_CALL;
				undef $initial_call;


				#------------------------------------------------------
				# snv vqsr training
				#------------------------------------------------------
				open FILE_HANDLE_VQSR, ">$Process_Path{$sub_argv}/shells/step_2-1-2_vqsr_snv_$Method_SNV.$runchr_group.sh" or die $!;
				my $initial_vqsr_on  = <<"VQSR_ON";
#!/bin/sh
## BaseRecalibrator for snv $sub_argv
#\$ -cwd -l virtual_free=30G,h_vmem=50G -pe smp 8

$ReadCfg{Java} -Xms8G -Xmx10G -Djava.io.tmpdir=/mnt/isilon/cbttc/tester/javaTmp/$sub_argv/SelectVariants -jar $ReadCfg{GATK} \\
	-T SelectVariants \\
	-R $ReadCfg{$Method_REF} \\
	-V $Process_Path{$sub_argv}/runs/step_2-1_snv/$sub_argv.$runchr_group.HaplotypeCaller.vcf \\
	-L $ReadCfg{"$Method_REF.root"}/newDefine/$runchr_group.intervals \\
	-selectType SNP \\
	-o $Process_Path{$sub_argv}/runs/step_2-1_snv/$sub_argv.$runchr_group.snv.vcf || { while [ -e $Analysis_Log{$sub_argv}.tmp ]; do sleep 1; done; perl $Bin/UpdateLog.pl $Analysis_Log{$sub_argv} "vqsr snv $runchr_group" 3 > $Analysis_Log{$sub_argv}.tmp && mv $Analysis_Log{$sub_argv}.tmp $Analysis_Log{$sub_argv}; exit 1; }

$ReadCfg{Java} -Xmx10G -Xms8G -Djava.io.tmpdir=/mnt/isilon/cbttc/tester/javaTmp/$sub_argv/VariantRecalibrator -jar $ReadCfg{GATK} \\
	-T VariantRecalibrator \\
	-R $ReadCfg{$Method_REF} \\
	-input $Process_Path{$sub_argv}/runs/step_2-1_snv/$sub_argv.$runchr_group.snv.vcf \\
	-recalFile $Process_Path{$sub_argv}/runs/step_2-1_snv/$sub_argv.$runchr_group.snv.recal \\
	-tranchesFile $Process_Path{$sub_argv}/runs/step_2-1_snv/$sub_argv.$runchr_group.snv.tranches \\
	-L $ReadCfg{"$Method_REF.root"}/newDefine/$runchr_group.intervals \\
	-nt 16 \\
	-mode SNP \\
	$ReadCfg{'VariantRecalibrator.snp'} \\
	-resource:hapmap,known=false,training=true,truth=true,prior=15.0 $ReadCfg{'hapmap_3.3.b37'} \\
	-resource:omni,known=false,training=true,truth=true,prior=12.0 $ReadCfg{'1000G_omni2.5.b37'} \\
	-resource:1000G,known=false,training=true,truth=false,prior=10.0 $ReadCfg{'1000G_p1.snps.hq'} \\
	-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $ReadCfg{'snp138.GRCh37'} || { while [ -e $Analysis_Log{$sub_argv}.tmp ]; do sleep 1; done; perl $Bin/UpdateLog.pl $Analysis_Log{$sub_argv} "vqsr snv $runchr_group" 3 > $Analysis_Log{$sub_argv}.tmp && mv $Analysis_Log{$sub_argv}.tmp $Analysis_Log{$sub_argv}; exit 1; }

$ReadCfg{Java} -Xmx10G -Xms8G -Djava.io.tmpdir=/mnt/isilon/cbttc/tester/javaTmp/$sub_argv/ApplyRecalibration -jar $ReadCfg{GATK} \\
	-T ApplyRecalibration \\
	-R $ReadCfg{$Method_REF} \\
	-input $Process_Path{$sub_argv}/runs/step_2-1_snv/$sub_argv.$runchr_group.snv.vcf \\
	-recalFile $Process_Path{$sub_argv}/runs/step_2-1_snv/$sub_argv.$runchr_group.snv.recal \\
	-tranchesFile $Process_Path{$sub_argv}/runs/step_2-1_snv/$sub_argv.$runchr_group.snv.tranches \\
	-L $ReadCfg{"$Method_REF.root"}/newDefine/$runchr_group.intervals \\
	-o $Process_Path{$sub_argv}/runs/step_2-1_snv/$sub_argv.$runchr_group.recalibrated.snv.vcf \\
	--ts_filter_level 99.5 \\
	-mode SNP --excludeFiltered || { while [ -e $Analysis_Log{$sub_argv}.tmp ]; do sleep 1; done; perl $Bin/UpdateLog.pl $Analysis_Log{$sub_argv} "vqsr snv $runchr_group" 3 > $Analysis_Log{$sub_argv}.tmp && mv $Analysis_Log{$sub_argv}.tmp $Analysis_Log{$sub_argv}; exit 1; }

while [ -e $Analysis_Log{$sub_argv}.tmp ]; do sleep 1; done; perl $Bin/UpdateLog.pl $Analysis_Log{$sub_argv} "vqsr snv $runchr_group" 0 > $Analysis_Log{$sub_argv}.tmp && mv $Analysis_Log{$sub_argv}.tmp $Analysis_Log{$sub_argv}

while [ -e $Analysis_Ckp{$sub_argv}.tmp ]; do sleep 1; done; perl $Bin/UpdateCkp.pl $Analysis_Ckp{$sub_argv} snv.$runchr_group $Method_SNV 4 > $Analysis_Ckp{$sub_argv}.tmp && mv $Analysis_Ckp{$sub_argv}.tmp $Analysis_Ckp{$sub_argv}

perl $Bin/JobProcessor.pl $OUTDIR/$sub_argv $Analysis_Ckp{$sub_argv} snv.$runchr_group $Method_SNV
VQSR_ON
				my $initial_vqsr_off = <<"VQSR_OFF";
#!/bin/sh
## BaseRecalibrator for snv $sub_argv
#\$ -cwd -l virtual_free=15G,h_vmem=30G -pe smp 4

$ReadCfg{Java} -Xms8G -Xmx10G -Djava.io.tmpdir=/mnt/isilon/cbttc/tester/javaTmp/$sub_argv/SelectVariants -jar $ReadCfg{GATK} \\
	-T SelectVariants \\
	-R $ReadCfg{$Method_REF} \\
	-V $Process_Path{$sub_argv}/runs/step_2-1_snv/$sub_argv.$runchr_group.HaplotypeCaller.vcf \\
	-L $ReadCfg{"$Method_REF.root"}/newDefine/$runchr_group.intervals \\
	-selectType SNP \\
	-o $Process_Path{$sub_argv}/runs/step_2-1_snv/$sub_argv.$runchr_group.snv.vcf || { while [ -e $Analysis_Log{$sub_argv}.tmp ]; do sleep 1; done; perl $Bin/UpdateLog.pl $Analysis_Log{$sub_argv} "vqsr snv $runchr_group" 3 > $Analysis_Log{$sub_argv}.tmp && mv $Analysis_Log{$sub_argv}.tmp $Analysis_Log{$sub_argv}; exit 1; }

ln $Process_Path{$sub_argv}/runs/step_2-1_snv/$sub_argv.$runchr_group.snv.vcf $Process_Path{$sub_argv}/runs/step_2-1_snv/$sub_argv.$runchr_group.recalibrated.snv.vcf

while [ -e $Analysis_Log{$sub_argv}.tmp ]; do sleep 1; done; perl $Bin/UpdateLog.pl $Analysis_Log{$sub_argv} "vqsr snv $runchr_group" 0 > $Analysis_Log{$sub_argv}.tmp && mv $Analysis_Log{$sub_argv}.tmp $Analysis_Log{$sub_argv}

while [ -e $Analysis_Ckp{$sub_argv}.tmp ]; do sleep 1; done; perl $Bin/UpdateCkp.pl $Analysis_Ckp{$sub_argv} snv.$runchr_group $Method_SNV 4 > $Analysis_Ckp{$sub_argv}.tmp && mv $Analysis_Ckp{$sub_argv}.tmp $Analysis_Ckp{$sub_argv}

perl $Bin/JobProcessor.pl $OUTDIR/$sub_argv $Analysis_Ckp{$sub_argv} snv.$runchr_group $Method_SNV
VQSR_OFF
				if ( !defined $ReadCfg{'GATK.VQSR'} ) {$ReadCfg{'GATK.VQSR'} = 'off';}
				
				if ( $ReadCfg{'GATK.VQSR'} eq 'on' ){
					print FILE_HANDLE_VQSR $initial_vqsr_on;
				} elsif ( $ReadCfg{'GATK.VQSR'} eq 'off' ) {
					print FILE_HANDLE_VQSR $initial_vqsr_off;
				} else {
					die "bad vqsr setting in config file, only \'on\' & \'off\' allow to use, fix to run\n";
				}
				close FILE_HANDLE_VQSR;
				undef $initial_vqsr_off;
				undef $initial_vqsr_on;

				#------------------------------------------------------
				# snv annotation
				#------------------------------------------------------
				open FILE_HANDLE_ANN, ">$Process_Path{$sub_argv}/shells/step_2-1-4_ann_snv_$Method_SNV.$runchr_group.sh" or die $!;
				my $initial_ann;
				switch ( $Species ) {
					case 'human' {
						$initial_ann = <<"SCRIPTS";
#!/bin/sh
## annotation $Method_ANN for snv $sub_argv
#\$ -cwd -l virtual_free=15G,h_vmem=30G -pe smp 4

$ReadCfg{Java} -Xmx10G -Xms8G -Djava.io.tmpdir=/mnt/isilon/cbttc/tester/javaTmp/$sub_argv/SelectVariants -jar $ReadCfg{$Method_SNV} \\
	-T SelectVariants \\
	-R $ReadCfg{$Method_REF} \\
	--discordance $ReadCfg{'hapmap_3.3.b37'} \\
	-V $Process_Path{$sub_argv}/runs/step_2-1_snv/$sub_argv.$runchr_group.recalibrated.snv.vcf \\
	-L $ReadCfg{"$Method_REF.root"}/newDefine/$runchr_group.intervals \\
	-o $Process_Path{$sub_argv}/runs/step_2-1_snv/$sub_argv.$runchr_group.rmsnp1.snv.vcf || { while [ -e $Analysis_Log{$sub_argv}.tmp ]; do sleep 1; done; perl $Bin/UpdateLog.pl $Analysis_Log{$sub_argv} "ann snv $runchr_group" 3 > $Analysis_Log{$sub_argv}.tmp && mv $Analysis_Log{$sub_argv}.tmp $Analysis_Log{$sub_argv}; exit 1; }

$ReadCfg{Java} -Xmx10G -Xms8G -Djava.io.tmpdir=/mnt/isilon/cbttc/tester/javaTmp/$sub_argv/SelectVariants -jar $ReadCfg{$Method_SNV} \\
	-T SelectVariants \\
	-R $ReadCfg{$Method_REF} \\
	--discordance $ReadCfg{'1000G_omni2.5.b37'} \\
	-L $ReadCfg{"$Method_REF.root"}/newDefine/$runchr_group.intervals \\
	-V $Process_Path{$sub_argv}/runs/step_2-1_snv/$sub_argv.$runchr_group.rmsnp1.snv.vcf \\
	-o $Process_Path{$sub_argv}/runs/step_2-1_snv/$sub_argv.$runchr_group.rmsnp2.snv.vcf || { while [ -e $Analysis_Log{$sub_argv}.tmp ]; do sleep 1; done; perl $Bin/UpdateLog.pl $Analysis_Log{$sub_argv} "ann snv $runchr_group" 3 > $Analysis_Log{$sub_argv}.tmp && mv $Analysis_Log{$sub_argv}.tmp $Analysis_Log{$sub_argv}; exit 1; }

$ReadCfg{Java} -Xmx10G -Xms8G -Djava.io.tmpdir=/mnt/isilon/cbttc/tester/javaTmp/$sub_argv/SelectVariants -jar $ReadCfg{$Method_SNV} \\
	-T SelectVariants \\
	-R $ReadCfg{$Method_REF} \\
	--discordance $ReadCfg{'1000G_p1.snps.hq'} \\
	-L $ReadCfg{"$Method_REF.root"}/newDefine/$runchr_group.intervals \\
	-V $Process_Path{$sub_argv}/runs/step_2-1_snv/$sub_argv.$runchr_group.rmsnp2.snv.vcf \\
	-o $Process_Path{$sub_argv}/runs/step_2-1_snv/$sub_argv.$runchr_group.rmsnp.snv.vcf || { while [ -e $Analysis_Log{$sub_argv}.tmp ]; do sleep 1; done; perl $Bin/UpdateLog.pl $Analysis_Log{$sub_argv} "ann snv $runchr_group" 3 > $Analysis_Log{$sub_argv}.tmp && mv $Analysis_Log{$sub_argv}.tmp $Analysis_Log{$sub_argv}; exit 1; }

$ReadCfg{Java} -Xms4G -Xmx8G -Djava.io.tmpdir=/mnt/isilon/cbttc/tester/javaTmp/$sub_argv/snpeff -jar $ReadCfg{snpEff} \\
	-c $ReadCfg{snpEffCfg} \\
	$ReadCfg{'snpEff.human'} \\
	$Process_Path{$sub_argv}/runs/step_2-1_snv/$sub_argv.$runchr_group.rmsnp.snv.vcf > $Process_Path{$sub_argv}/runs/step_2-1_snv/$sub_argv.$runchr_group.snv.snpeff.vcf || { while [ -e $Analysis_Log{$sub_argv}.tmp ]; do sleep 1; done; perl $Bin/UpdateLog.pl $Analysis_Log{$sub_argv} "ann snv $runchr_group" 3 > $Analysis_Log{$sub_argv}.tmp && mv $Analysis_Log{$sub_argv}.tmp $Analysis_Log{$sub_argv}; exit 1; }

$ReadCfg{Java} -Xms4G -Xmx8G -Djava.io.tmpdir=/mnt/isilon/cbttc/tester/javaTmp/$sub_argv/snpSift -jar $ReadCfg{snpSift} \\
	dbnsfp -v $ReadCfg{dbNSFP} \\
$Process_Path{$sub_argv}/runs/step_2-1_snv/$sub_argv.$runchr_group.snv.snpeff.vcf > $Process_Path{$sub_argv}/runs/step_2-1_snv/$sub_argv.$runchr_group.dbnsfp.snv.vcf || { while [ -e $Analysis_Log{$sub_argv}.tmp ]; do sleep 1; done; perl $Bin/UpdateLog.pl $Analysis_Log{$sub_argv} "ann snv $runchr_group" 3 > $Analysis_Log{$sub_argv}.tmp && mv $Analysis_Log{$sub_argv}.tmp $Analysis_Log{$sub_argv}; exit 1; }

while [ -e $Analysis_Log{$sub_argv}.tmp ]; do sleep 1; done; perl $Bin/UpdateLog.pl $Analysis_Log{$sub_argv} "ann snv $runchr_group" 0 > $Analysis_Log{$sub_argv}.tmp && mv $Analysis_Log{$sub_argv}.tmp $Analysis_Log{$sub_argv}

while [ -e $Analysis_Ckp{$sub_argv}.tmp ]; do sleep 1; done; perl $Bin/UpdateCkp.pl $Analysis_Ckp{$sub_argv} snv.$runchr_group $Method_SNV 5 > $Analysis_Ckp{$sub_argv}.tmp && mv $Analysis_Ckp{$sub_argv}.tmp $Analysis_Ckp{$sub_argv}

perl $Bin/JobProcessor.pl $OUTDIR/$sub_argv $Analysis_Ckp{$sub_argv} snv.$runchr_group $Method_SNV
SCRIPTS
					}
					
					case 'mouse' {
						$initial_ann = <<"SCRIPTS";
#!/bin/sh
## annotation $Method_ANN for snv $sub_argv
#\$ -cwd -l virtual_free=15G,h_vmem=30G -pe smp 4

$ReadCfg{Java} -Xmx10G -Xms8G -Djava.io.tmpdir=/mnt/isilon/cbttc/tester/javaTmp/$sub_argv/SelectVariants -jar $ReadCfg{$Method_SNV} \\
	-T SelectVariants \\
	-R $ReadCfg{$Method_REF} \\
	--discordance $ReadCfg{'sanger.snp'} \\
	-V $Process_Path{$sub_argv}/runs/step_2-1_snv/$sub_argv.$runchr_group.recalibrated.snv.vcf \\
	-L $ReadCfg{"$Method_REF.root"}/newDefine/$runchr_group.intervals \\
	-o $Process_Path{$sub_argv}/runs/step_2-1_snv/$sub_argv.$runchr_group.rmsnp.snv.vcf || { while [ -e $Analysis_Log{$sub_argv}.tmp ]; do sleep 1; done; perl $Bin/UpdateLog.pl $Analysis_Log{$sub_argv} "ann snv $runchr_group" 3 > $Analysis_Log{$sub_argv}.tmp && mv $Analysis_Log{$sub_argv}.tmp $Analysis_Log{$sub_argv}; exit 1; }

$ReadCfg{Java} -Xms4G -Xmx8G -Djava.io.tmpdir=/mnt/isilon/cbttc/tester/javaTmp/$sub_argv/snpEff -jar $ReadCfg{snpEff} \\
	-c $ReadCfg{snpEffCfg} \\
	$ReadCfg{'snpEff.mouse'} \\
	$Process_Path{$sub_argv}/runs/step_2-1_snv/$sub_argv.$runchr_group.rmsnp.snv.vcf > $Process_Path{$sub_argv}/runs/step_2-1_snv/$sub_argv.$runchr_group.snv.snpeff.vcf || { while [ -e $Analysis_Log{$sub_argv}.tmp ]; do sleep 1; done; perl $Bin/UpdateLog.pl $Analysis_Log{$sub_argv} "ann snv $runchr_group" 3 > $Analysis_Log{$sub_argv}.tmp && mv $Analysis_Log{$sub_argv}.tmp $Analysis_Log{$sub_argv}; exit 1; }

while [ -e $Analysis_Log{$sub_argv}.tmp ]; do sleep 1; done; perl $Bin/UpdateLog.pl $Analysis_Log{$sub_argv} "ann snv $runchr_group" 0 > $Analysis_Log{$sub_argv}.tmp && mv $Analysis_Log{$sub_argv}.tmp $Analysis_Log{$sub_argv}

while [ -e $Analysis_Ckp{$sub_argv}.tmp ]; do sleep 1; done; perl $Bin/UpdateCkp.pl $Analysis_Ckp{$sub_argv} snv.$runchr_group $Method_SNV 5 > $Analysis_Ckp{$sub_argv}.tmp && mv $Analysis_Ckp{$sub_argv}.tmp $Analysis_Ckp{$sub_argv}

perl $Bin/JobProcessor.pl $OUTDIR/$sub_argv $Analysis_Ckp{$sub_argv} snv.$runchr_group $Method_SNV
SCRIPTS
					}
				}
				print FILE_HANDLE_ANN $initial_ann;
				close FILE_HANDLE_ANN;
				undef $initial_ann;
			}
			
			case ('muTect') {
				die "not cancer mode, bad snv method setting!" unless $CancerMode == 1;
				my $normal_sample = $NormalSample{$PairMarker{$sub_argv}};
				my $tumor_sample = $TumorSample{$PairMarker{$sub_argv}};
				#------------------------------------------------------
				# check cancer pair status 
				#------------------------------------------------------
				open FILE_HANDLE_CALL, ">$Process_Path{$sub_argv}/shells/step_2-1-1_call_snv_$Method_SNV.$runchr_group.sh" or die $!;
				my $initial_call = <<"SCRIPTS";
#!/bin/sh
## call somatic mutation $PairMarker{$sub_argv}
#\$ -cwd -l virtual_free=0.1G,h_vmem=0.1G

# check pair status
perl $Bin/CheckPair.pl $Analysis_Ckp{$normal_sample} $Analysis_Ckp{$tumor_sample} bqsr.$runchr_group 5 || { echo "Pair mate $PairLink{$sub_argv} is not ready" && exit 1; }

while [ -e $Analysis_Ckp{$normal_sample}.tmp ]; do sleep 1; done; perl $Bin/UpdateCkp.pl $Analysis_Ckp{$normal_sample} snv.$runchr_group $Method_SNV 3 > $Analysis_Ckp{$normal_sample}.tmp && mv $Analysis_Ckp{$normal_sample}.tmp $Analysis_Ckp{$normal_sample}

while [ -e $Analysis_Log{$normal_sample}.tmp ]; do sleep 1; done; perl $Bin/UpdateLog.pl $Analysis_Log{$normal_sample} "call snv $runchr_group" 0 > $Analysis_Log{$normal_sample}.tmp && mv $Analysis_Log{$normal_sample}.tmp $Analysis_Log{$normal_sample}
while [ -e $Analysis_Log{$tumor_sample}.tmp ]; do sleep 1; done; perl $Bin/UpdateLog.pl $Analysis_Log{$tumor_sample} "call snv $runchr_group" 0 > $Analysis_Log{$tumor_sample}.tmp && mv $Analysis_Log{$tumor_sample}.tmp $Analysis_Log{$tumor_sample}

# next step
perl $Bin/JobProcessor.pl $OUTDIR/$normal_sample $Analysis_Ckp{$normal_sample} snv.$runchr_group $Method_SNV

SCRIPTS
				print FILE_HANDLE_CALL $initial_call;
				close FILE_HANDLE_CALL;
				undef $initial_call;
				
				#------------------------------------------------------
				# muTect call somatic
				#------------------------------------------------------
				open FILE_HANDLE_SOMATIC, ">$Process_Path{$sub_argv}/shells/step_2-1-3_somatic_snv_$Method_SNV.$runchr_group.sh" or die $!;
				my $initial_somatic = <<"SCRIPTS";
#!/bin/sh
## call somatic mutation $PairMarker{$normal_sample}
#\$ -cwd -l virtual_free=15G,h_vmem=20G -pe smp 8

# check pair status
# perl $Bin/CheckPair.pl $Analysis_Ckp{$normal_sample} $Analysis_Ckp{$tumor_sample} bqsr.$runchr_group 5 || { echo "Pair mate $PairLink{$sub_argv} is not ready" && exit 1; }

$ReadCfg{'Java1.6'} -Xms8G -Xmx10G -Djava.io.tmpdir=/mnt/isilon/cbttc/tester/javaTmp/$sub_argv/muTect -jar $ReadCfg{muTect} \\
	--analysis_type MuTect \\
	--reference_sequence $ReadCfg{$Method_REF} \\
	--dbsnp $ReadCfg{"snp138.$Method_REF"} \\
	--cosmic $ReadCfg{cosmic} \\
	--input_file:normal $Process_Path{$normal_sample}/runs/step_1-1_map/$normal_sample.$runchr_group.bam \\
	--input_file:tumor $Process_Path{$tumor_sample}/runs/step_1-1_map/$tumor_sample.$runchr_group.bam \\
	--out $Process_Path{$PairMarker{$sub_argv}}/step_2-1_snv/$runchr_group.muTect.txt \\
	--vcf $Process_Path{$PairMarker{$sub_argv}}/step_2-1_snv/$runchr_group.muTect.raw.vcf \\
	-L $ReadCfg{"$Method_REF.root"}/newDefine/$runchr_group.intervals \\
	-nt 8 || { while [ -e $Analysis_Log{$normal_sample}.tmp ]; do sleep 1; done; perl $Bin/UpdateLog.pl $Analysis_Log{$normal_sample} "somatic snv $runchr_group" 3 > $Analysis_Log{$normal_sample}.tmp && mv $Analysis_Log{$normal_sample}.tmp $Analysis_Log{$normal_sample}; exit 1; }
# grep -v REJECT $Process_Path{$PairMarker{$sub_argv}}/step_2-1_snv/$runchr_group.muTect.raw.vcf > $Process_Path{$PairMarker{$sub_argv}}/step_2-1_snv/$runchr_group.muTect.pass.vcf

# mark pair done
# while [ -e $Analysis_Log{$PairMarker{$sub_argv}}.tmp ]; do sleep 1; done; perl $Bin/UpdateLog.pl $Analysis_Log{$PairMarker{$sub_argv}} "somatic snv $runchr_group" 0 > $Analysis_Log{$PairMarker{$sub_argv}}.tmp && mv $Analysis_Log{$PairMarker{$sub_argv}}.tmp $Analysis_Log{$PairMarker{$sub_argv}}
# while [ -e $Analysis_Ckp{$PairMarker{$sub_argv}}.tmp ]; do sleep 1; done; perl $Bin/UpdateCkp.pl $Analysis_Ckp{$PairMarker{$sub_argv}} snv.$runchr_group $Method_SNV 3 > $Analysis_Ckp{$PairMarker{$sub_argv}}.tmp && mv $Analysis_Ckp{$PairMarker{$sub_argv}}.tmp $Analysis_Ckp{$PairMarker{$sub_argv}}

# mark somatic done
while [ -e $Analysis_Log{$normal_sample}.tmp ]; do sleep 1; done; perl $Bin/UpdateLog.pl $Analysis_Log{$normal_sample} "somatic snv  $runchr_group" 0 > $Analysis_Log{$normal_sample}.tmp && mv $Analysis_Log{$normal_sample}.tmp $Analysis_Log{$normal_sample}

while [ -e $Analysis_Ckp{$normal_sample}.tmp ]; do sleep 1; done; perl $Bin/UpdateCkp.pl $Analysis_Ckp{$normal_sample} snv.$runchr_group $Method_SNV 4 > $Analysis_Ckp{$normal_sample}.tmp && mv $Analysis_Ckp{$normal_sample}.tmp $Analysis_Ckp{$normal_sample}

# next step
perl $Bin/JobProcessor.pl $OUTDIR/$normal_sample $Analysis_Ckp{$normal_sample} snv.$runchr_group $Method_SNV

SCRIPTS
				print FILE_HANDLE_SOMATIC $initial_somatic;
				close FILE_HANDLE_SOMATIC;
				undef $initial_somatic;


				#------------------------------------------------------
				# snv annotation
				#------------------------------------------------------
				open FILE_HANDLE_ANN, ">$Process_Path{$sub_argv}/shells/step_2-1-4_ann_snv_$Method_SNV.$runchr_group.sh" or die $!;
				my $initial_ann = <<"SCRIPTS";
#!/bin/sh
## annotation $Method_ANN for snv $sub_argv
#\$ -cwd -l virtual_free=10G,h_vmem=15G -pe smp 3

$ReadCfg{Java} -Xmx8G -Xms5G -Djava.io.tmpdir=/mnt/isilon/cbttc/tester/javaTmp/$sub_argv/SelectVariants -jar $ReadCfg{GATK} \\
	-T SelectVariants \\
	-R $ReadCfg{$Method_REF} \\
	--discordance $ReadCfg{'hapmap_3.3.b37'} \\
	-V $Process_Path{$PairMarker{$sub_argv}}/step_2-1_snv/$runchr_group.muTect.raw.vcf \\
	-L $ReadCfg{"$Method_REF.root"}/newDefine/$runchr_group.intervals \\
	-o $Process_Path{$PairMarker{$sub_argv}}/step_2-1_snv/$runchr_group.muTect.rmsnp1.snv.vcf || { while [ -e $Analysis_Log{$normal_sample}.tmp ]; do sleep 1; done; perl $Bin/UpdateLog.pl $Analysis_Log{$normal_sample} "ann snv $runchr_group" 3 > $Analysis_Log{$normal_sample}.tmp && mv $Analysis_Log{$normal_sample}.tmp $Analysis_Log{$normal_sample}; exit 1; }
	
$ReadCfg{Java} -Xmx8G -Xms5G -Djava.io.tmpdir=/mnt/isilon/cbttc/tester/javaTmp/$sub_argv/SelectVariants -jar $ReadCfg{GATK} \\
	-T SelectVariants \\
	-R $ReadCfg{$Method_REF} \\
	--discordance $ReadCfg{'1000G_omni2.5.b37'} \\
	-V $Process_Path{$PairMarker{$sub_argv}}/step_2-1_snv/$runchr_group.muTect.rmsnp1.snv.vcf\\
	-L $ReadCfg{"$Method_REF.root"}/newDefine/$runchr_group.intervals \\
	-o $Process_Path{$PairMarker{$sub_argv}}/step_2-1_snv/$runchr_group.muTect.rmsnp2.snv.vcf || { while [ -e $Analysis_Log{$normal_sample}.tmp ]; do sleep 1; done; perl $Bin/UpdateLog.pl $Analysis_Log{$normal_sample} "ann snv $runchr_group" 3 > $Analysis_Log{$normal_sample}.tmp && mv $Analysis_Log{$normal_sample}.tmp $Analysis_Log{$normal_sample}; exit 1; }
	
$ReadCfg{Java} -Xmx8G -Xms5G -Djava.io.tmpdir=/mnt/isilon/cbttc/tester/javaTmp/$sub_argv/SelectVariants -jar $ReadCfg{GATK} \\
	-T SelectVariants \\
	-R $ReadCfg{$Method_REF} \\
	--discordance $ReadCfg{'1000G_p1.snps.hq'} \\
	-V $Process_Path{$PairMarker{$sub_argv}}/step_2-1_snv/$runchr_group.muTect.rmsnp2.snv.vcf \\
	-L $ReadCfg{"$Method_REF.root"}/newDefine/$runchr_group.intervals \\
	-o $Process_Path{$PairMarker{$sub_argv}}/step_2-1_snv/$runchr_group.muTect.rmsnp.snv.vcf || { while [ -e $Analysis_Log{$normal_sample}.tmp ]; do sleep 1; done; perl $Bin/UpdateLog.pl $Analysis_Log{$normal_sample} "ann snv $runchr_group" 3 > $Analysis_Log{$normal_sample}.tmp && mv $Analysis_Log{$normal_sample}.tmp $Analysis_Log{$normal_sample}; exit 1; }

$ReadCfg{Java} -Xms4G -Xmx8G -Djava.io.tmpdir=/mnt/isilon/cbttc/tester/javaTmp/$sub_argv/snpEff -jar $ReadCfg{snpEff} \\
	-c $ReadCfg{snpEffCfg} \\
	$ReadCfg{'snpEff.human'} \\
	$Process_Path{$PairMarker{$sub_argv}}/step_2-1_snv/$runchr_group.muTect.rmsnp.snv.vcf > $Process_Path{$PairMarker{$sub_argv}}/step_2-1_snv/$runchr_group.muTect.rmsnp.snv.snpeff.vcf || { while [ -e $Analysis_Log{$normal_sample}.tmp ]; do sleep 1; done; perl $Bin/UpdateLog.pl $Analysis_Log{$normal_sample} "ann snv $runchr_group" 3 > $Analysis_Log{$normal_sample}.tmp && mv $Analysis_Log{$normal_sample}.tmp $Analysis_Log{$normal_sample}; exit 1; }

$ReadCfg{Java} -Xms4G -Xmx8G -Djava.io.tmpdir=/mnt/isilon/cbttc/tester/javaTmp/$sub_argv/snpSift -jar $ReadCfg{snpSift} \\
	dbnsfp -v $ReadCfg{dbNSFP} \\
$Process_Path{$PairMarker{$sub_argv}}/step_2-1_snv/$runchr_group.muTect.rmsnp.snv.snpeff.vcf >$Process_Path{$PairMarker{$sub_argv}}/step_2-1_snv/$runchr_group.muTect.dbnsfp.snv.vcf || { while [ -e $Analysis_Log{$normal_sample}.tmp ]; do sleep 1; done; perl $Bin/UpdateLog.pl $Analysis_Log{$normal_sample} "ann snv $runchr_group" 3 > $Analysis_Log{$normal_sample}.tmp && mv $Analysis_Log{$normal_sample}.tmp $Analysis_Log{$normal_sample}; exit 1; }

while [ -e $Analysis_Log{$normal_sample}.tmp ]; do sleep 1; done; perl $Bin/UpdateLog.pl $Analysis_Log{$normal_sample} "ann snv $runchr_group" 0 > $Analysis_Log{$normal_sample}.tmp && mv $Analysis_Log{$normal_sample}.tmp $Analysis_Log{$normal_sample}

while [ -e $Analysis_Ckp{$normal_sample}.tmp ]; do sleep 1; done; perl $Bin/UpdateCkp.pl $Analysis_Ckp{$normal_sample} snv.$runchr_group $Method_SNV 5 > $Analysis_Ckp{$normal_sample}.tmp && mv $Analysis_Ckp{$normal_sample}.tmp $Analysis_Ckp{$normal_sample}

# perl $Bin/JobProcessor.pl $OUTDIR/$sub_argv $Analysis_Ckp{$sub_argv} snv.$runchr_group $Method_SNV
SCRIPTS
				print FILE_HANDLE_ANN $initial_ann;
				close FILE_HANDLE_ANN;
				undef $initial_ann;
			}
			
			case ('VarScan') {
				die "not cancer mode, bad snv method setting!" unless $CancerMode == 1;
				my $normal_sample = $NormalSample{$PairMarker{$sub_argv}};
				my $tumor_sample = $TumorSample{$PairMarker{$sub_argv}};
				#------------------------------------------------------
				# mpileup && check cancer pair status 
				#------------------------------------------------------
				open FILE_HANDLE_CALL, ">$Process_Path{$sub_argv}/shells/step_2-1-1_call_snv_$Method_SNV.$runchr_group.sh" or die $!;
				my $initial_call = <<"SCRIPTS";
#!/bin/sh
## call somatic InDel $PairMarker{$sub_argv}
#\$ -cwd -l virtual_free=1G,h_vmem=1G

if [ ! -e $Process_Path{$sub_argv}/runs/step_1-1_map/$sub_argv.$runchr_group.pileup ]
then
	$ReadCfg{Samtools} mpileup -f $ReadCfg{$Method_REF} $Process_Path{$sub_argv}/runs/step_1-1_map/$sub_argv.$runchr_group.bam > $Process_Path{$sub_argv}/runs/step_1-1_map/$sub_argv.$runchr_group.pileup || { echo "mpileup $sub_argv failed"; while [ -e $Analysis_Log{$sub_argv}.tmp ]; do sleep 1; done; perl $Bin/UpdateLog.pl $Analysis_Log{$sub_argv} "call snv $runchr_group" 3 > $Analysis_Log{$sub_argv}.tmp && mv $Analysis_Log{$sub_argv}.tmp $Analysis_Log{$sub_argv}; exit 1; }
fi

while [ -e $Analysis_Ckp{$normal_sample}.tmp ]; do sleep 1; done; perl $Bin/UpdateCkp.pl $Analysis_Ckp{$normal_sample} snv.$runchr_group $Method_SNV 3 > $Analysis_Ckp{$normal_sample}.tmp && mv $Analysis_Ckp{$normal_sample}.tmp $Analysis_Ckp{$normal_sample}

# check pair status
perl $Bin/CheckPair.pl $Analysis_Ckp{$normal_sample} $Analysis_Ckp{$tumor_sample} bqsr.$runchr_group 5 || { echo "Pair mate $PairLink{$sub_argv} is not ready" && exit 1; }

while [ -e $Analysis_Log{$normal_sample}}.tmp ]; do sleep 1; done; perl $Bin/UpdateLog.pl $Analysis_Log{$normal_sample} "call snv $runchr_group" 0 > $Analysis_Log{$normal_sample}.tmp && mv $Analysis_Log{$normal_sample}.tmp $Analysis_Log{$normal_sample}
while [ -e $Analysis_Log{$tumor_sample}.tmp ]; do sleep 1; done; perl $Bin/UpdateLog.pl $Analysis_Log{$tumor_sample} "call snv $runchr_group" 0 > $Analysis_Log{$tumor_sample}.tmp && mv $Analysis_Log{$tumor_sample}.tmp $Analysis_Log{$tumor_sample}

# next step
perl $Bin/JobProcessor.pl $OUTDIR/$normal_sample $Analysis_Ckp{$normal_sample} snv.$runchr_group $Method_SNV

SCRIPTS
				print FILE_HANDLE_CALL $initial_call;
				close FILE_HANDLE_CALL;
				undef $initial_call;
				
				
				#------------------------------------------------------
				# mpileup for both samples
				#------------------------------------------------------
				open FILE_HANDLE_CALL, ">$Process_Path{$sub_argv}/shells/step_2-2-3_call_snv_$Method_SNV.$runchr_group.sh" or die $!;
				my $initial_somatic = <<"SCRIPTS";
#!/bin/sh
## call somatic InDel $PairMarker{$sub_argv}
#\$ -cwd -l virtual_free=20G,h_vmem=35G

$ReadCfg{Java} -Xms5G -Xmx15G -Djava.io.tmpdir=/mnt/isilon/cbttc/tester/javaTmp/$sub_argv/VarScan -jar $ReadCfg{$Method_SNV} \\
	$ReadCfg{VarScanOpts} somatic \\
	$Process_Path{$normal_sample}/runs/step_1-1_map/$normal_sample.$runchr_group.pileup \\
	$Process_Path{$tumor_sample}/runs/step_1-1_map/$tumor_sample.$runchr_group.pileup \\
	$Process_Path{$PairMarker{$sub_argv}}/step_2-2_snv/$runchr_group.VarScan2.snv --output-vcf || { echo "varscan $sub_argv failed"; while [ -e $Analysis_Log{$sub_argv}.tmp ]; do sleep 1; done; perl $Bin/UpdateLog.pl $Analysis_Log{$sub_argv} "somatic snv $runchr_group" 3 > $Analysis_Log{$sub_argv}.tmp && mv $Analysis_Log{$sub_argv}.tmp $Analysis_Log{$sub_argv}; exit 1; }

while [ -e $Analysis_Log{$normal_sample}.tmp ]; do sleep 1; done; perl $Bin/UpdateLog.pl $Analysis_Log{$normal_sample} "somatic snv $runchr_group" 0 > $Analysis_Log{$normal_sample}.tmp && mv $Analysis_Log{$normal_sample}.tmp $Analysis_Log{$normal_sample}

while [ -e $Analysis_Ckp{$normal_sample}.tmp ]; do sleep 1; done; perl $Bin/UpdateCkp.pl $Analysis_Ckp{$normal_sample} snv.$runchr_group $Method_SNV 4 > $Analysis_Ckp{$normal_sample}.tmp && mv $Analysis_Ckp{$normal_sample}.tmp $Analysis_Ckp{$normal_sample}

perl $Bin/JobProcessor.pl $OUTDIR/$normal_sample $Analysis_Ckp{$normal_sample} snv.$runchr_group $Method_SNV
	
SCRIPTS
				print FILE_HANDLE_SOMATIC $initial_somatic;
				close FILE_HANDLE_SOMATIC;
				undef $initial_somatic;
				
				
				#------------------------------------------------------
				# snv annotation
				#------------------------------------------------------
				open FILE_HANDLE_ANN, ">$Process_Path{$sub_argv}/shells/step_2-2-4_ann_snv_$Method_SNV.$runchr_group.sh" or die $!;
				my $initial_ann = <<"SCRIPTS";
#!/bin/sh
## annotation $Method_ANN for snv $sub_argv
#\$ -cwd -l virtual_free=20G,h_vmem=30G

$ReadCfg{Java} -Xmx10G -Xms8G -Djava.io.tmpdir=/mnt/isilon/cbttc/tester/javaTmp/$sub_argv/SelectVariants -jar $ReadCfg{$Method_SNV} \\
	-T SelectVariants \\
	-R $ReadCfg{$Method_REF} \\
	--discordance $ReadCfg{'hapmap_3.3.b37'} \\
	-V $Process_Path{$PairMarker{$sub_argv}}/step_2-1_snv/$runchr_group.muTect.pass.vcf \\
	-L $ReadCfg{"$Method_REF.root"}/newDefine/$runchr_group.intervals \\
	-o $Process_Path{$PairMarker{$sub_argv}}/step_2-1_snv/$runchr_group.muTect.rmsnp1.snv.vcf || { while [ -e $Analysis_Log{$normal_sample}.tmp ]; do sleep 1; done; perl $Bin/UpdateLog.pl $Analysis_Log{$normal_sample} "ann snv $runchr_group" 3 > $Analysis_Log{$normal_sample}.tmp && mv $Analysis_Log{$normal_sample}.tmp $Analysis_Log{$normal_sample}; exit 1; }
	
$ReadCfg{Java} -Xmx10G -Xms8G -Djava.io.tmpdir=/mnt/isilon/cbttc/tester/javaTmp/$sub_argv/SelectVariants -jar $ReadCfg{$Method_SNV} \\
	-T SelectVariants \\
	-R $ReadCfg{$Method_REF} \\
	--discordance $ReadCfg{'1000G_omni2.5.b37'} \\
	-V $Process_Path{$PairMarker{$sub_argv}}/step_2-1_snv/$runchr_group.muTect.rmsnp1.snv.vcf\\
	-L $ReadCfg{"$Method_REF.root"}/newDefine/$runchr_group.intervals \\
	-o $Process_Path{$PairMarker{$sub_argv}}/step_2-1_snv/$runchr_group.muTect.rmsnp2.snv.vcf || { while [ -e $Analysis_Log{$normal_sample}.tmp ]; do sleep 1; done; perl $Bin/UpdateLog.pl $Analysis_Log{$normal_sample} "ann snv $runchr_group" 3 > $Analysis_Log{$normal_sample}.tmp && mv $Analysis_Log{$normal_sample}.tmp $Analysis_Log{$normal_sample}; exit 1; }
	
$ReadCfg{Java} -Xmx10G -Xms8G -Djava.io.tmpdir=/mnt/isilon/cbttc/tester/javaTmp/$sub_argv/SelectVariants -jar $ReadCfg{$Method_SNV} \\
	-T SelectVariants \\
	-R $ReadCfg{$Method_REF} \\
	--discordance $ReadCfg{'1000G_p1.snps.hq'} \\
	-V $Process_Path{$PairMarker{$sub_argv}}/step_2-1_snv/$runchr_group.muTect.rmsnp2.snv.vcf \\
	-L $ReadCfg{"$Method_REF.root"}/newDefine/$runchr_group.intervals \\
	-o $Process_Path{$PairMarker{$sub_argv}}/step_2-1_snv/$runchr_group.muTect.rmsnp.snv.vcf || { while [ -e $Analysis_Log{$normal_sample}.tmp ]; do sleep 1; done; perl $Bin/UpdateLog.pl $Analysis_Log{$normal_sample} "ann snv $runchr_group" 3 > $Analysis_Log{$normal_sample}.tmp && mv $Analysis_Log{$normal_sample}.tmp $Analysis_Log{$normal_sample}; exit 1; }

$ReadCfg{Java} -Xms4G -Xmx8G -Djava.io.tmpdir=/mnt/isilon/cbttc/tester/javaTmp/$sub_argv/snpEff -jar $ReadCfg{snpEff} \\
	-c $ReadCfg{snpEffCfg} \\
	-i vcf -o vcf $ReadCfg{snpEffdb} \\
	-v -lof -nextprot -motif \\
	$Process_Path{$PairMarker{$sub_argv}}/step_2-1_snv/$runchr_group.muTect.rmsnp.snv.vcf > $Process_Path{$PairMarker{$sub_argv}}/step_2-1_snv/$runchr_group.muTect.rmsnp.snv.snpeff.vcf || { while [ -e $Analysis_Log{$normal_sample}.tmp ]; do sleep 1; done; perl $Bin/UpdateLog.pl $Analysis_Log{$normal_sample} "ann snv $runchr_group" 3 > $Analysis_Log{$normal_sample}.tmp && mv $Analysis_Log{$normal_sample}.tmp $Analysis_Log{$normal_sample}; exit 1; }

$ReadCfg{Java} -Xms4G -Xmx8G -Djava.io.tmpdir=/mnt/isilon/cbttc/tester/javaTmp/snpSift -jar $ReadCfg{snpSift} \\
	dbnsfp -v $ReadCfg{dbNSFP} \\
$Process_Path{$PairMarker{$sub_argv}}/step_2-1_snv/$runchr_group.muTect.rmsnp.snv.snpeff.vcf >$Process_Path{$PairMarker{$sub_argv}}/step_2-1_snv/$runchr_group.muTect.dbnsfp.snv.vcf || { while [ -e $Analysis_Log{$normal_sample}.tmp ]; do sleep 1; done; perl $Bin/UpdateLog.pl $Analysis_Log{$normal_sample} "ann snv $runchr_group" 3 > $Analysis_Log{$normal_sample}.tmp && mv $Analysis_Log{$normal_sample}.tmp $Analysis_Log{$normal_sample}; exit 1; }

while [ -e $Analysis_Log{$normal_sample}.tmp ]; do sleep 1; done; perl $Bin/UpdateLog.pl $Analysis_Log{$normal_sample} "ann snv $runchr_group" 0 > $Analysis_Log{$normal_sample}.tmp && mv $Analysis_Log{$normal_sample}.tmp $Analysis_Log{$normal_sample}

while [ -e $Analysis_Ckp{$normal_sample}.tmp ]; do sleep 1; done; perl $Bin/UpdateCkp.pl $Analysis_Ckp{$normal_sample} snv.$runchr_group $Method_SNV 5 > $Analysis_Ckp{$normal_sample}.tmp && mv $Analysis_Ckp{$normal_sample}.tmp $Analysis_Ckp{$normal_sample}

# perl $Bin/JobProcessor.pl $OUTDIR/$sub_argv $Analysis_Ckp{$sub_argv} snv.$runchr_group $Method_SNV
SCRIPTS
				print FILE_HANDLE_ANN $initial_ann;
				close FILE_HANDLE_ANN;
				undef $initial_ann;
			}
		}
	}
	undef $sub_argv;
}

sub _idvScripts {
	my $sub_argv = shift;
	foreach my $runchr_group ( @RunChr_Group ) {
		#--------------------------------------------------------------
		# initial idv scripts
		#--------------------------------------------------------------
		switch ( $Method_IDV ) {
			case ('GATK') {
				#------------------------------------------------------
				# gatk HaplotypeCaller
				#------------------------------------------------------
				open FILE_HANDLE_CALL, ">$Process_Path{$sub_argv}/shells/step_2-2-1_call_idv_$Method_IDV.$runchr_group.sh" or die $!;
				my $initial_call = <<"SCRIPTS";
#!/bin/sh
## call idv for $sub_argv
#\$ -cwd -l virtual_free=1G,h_vmem=1G

echo 'waiting'

SCRIPTS
				print FILE_HANDLE_CALL $initial_call;
				close FILE_HANDLE_CALL;
				undef $initial_call;
	
	
				#------------------------------------------------------
				# idv vqsr training
				#------------------------------------------------------
				open FILE_HANDLE_VQSR, ">$Process_Path{$sub_argv}/shells/step_2-2-2_vqsr_idv_$Method_IDV.$runchr_group.sh" or die $!;
				my $initial_vqsr_on = <<"SCRIPTS";
#!/bin/sh
## BaseRecalibrator for idv $sub_argv
#\$ -cwd -l virtual_free=30G,h_vmem=50G

$ReadCfg{Java} -Xms8G -Xmx10G -Djava.io.tmpdir=/mnt/isilon/cbttc/tester/javaTmp/$sub_argv/SelectVariants -jar $ReadCfg{GATK} \\
	-T SelectVariants \\
	-R $ReadCfg{$Method_REF} \\
	-V $Process_Path{$sub_argv}/runs/step_2-1_snv/$sub_argv.$runchr_group.HaplotypeCaller.vcf \\
	-L $ReadCfg{"$Method_REF.root"}/newDefine/$runchr_group.intervals \\
	-selectType INDEL \\
	-o $Process_Path{$sub_argv}/runs/step_2-2_idv/$sub_argv.$runchr_group.idv.vcf || { while [ -e $Analysis_Log{$sub_argv}.tmp ]; do sleep 1; done; perl $Bin/UpdateLog.pl $Analysis_Log{$sub_argv} "vqsr idv $runchr_group" 3 > $Analysis_Log{$sub_argv}.tmp && mv $Analysis_Log{$sub_argv}.tmp $Analysis_Log{$sub_argv}; exit 1; }

$ReadCfg{Java} -Xmx10G -Xms8G -Djava.io.tmpdir=/mnt/isilon/cbttc/tester/javaTmp/$sub_argv/SelectVariants -jar $ReadCfg{GATK} \\
	-T VariantRecalibrator \\
	-R $ReadCfg{$Method_REF} \\
	-input $Process_Path{$sub_argv}/runs/step_2-2_idv/$sub_argv.$runchr_group.idv.vcf \\
	-recalFile $Process_Path{$sub_argv}/runs/step_2-2_idv/$sub_argv.$runchr_group.idv.recal \\
	-tranchesFile $Process_Path{$sub_argv}/runs/step_2-2_idv/$sub_argv.$runchr_group.idv.tranches \\
	-L $ReadCfg{"$Method_REF.root"}/newDefine/$runchr_group.intervals \\
	-nt 16 \\
	-mode INDEL \\
	$ReadCfg{'VariantRecalibrator.idv'} \\
	-resource:mills,known=false,training=true,truth=true,prior=12.0 $ReadCfg{'mills.indel'} \\
	-resource:1000G,known=true,training=true,truth=false,prior=10.0 $ReadCfg{'1000g.indel'} \\
	-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $ReadCfg{'snp138.GRCh37'} || { while [ -e $Analysis_Log{$sub_argv}.tmp ]; do sleep 1; done; perl $Bin/UpdateLog.pl $Analysis_Log{$sub_argv} "vqsr idv $runchr_group" 3 > $Analysis_Log{$sub_argv}.tmp && mv $Analysis_Log{$sub_argv}.tmp $Analysis_Log{$sub_argv}; exit 1; }

$ReadCfg{Java} -Xmx10G -Xms8G -Djava.io.tmpdir=/mnt/isilon/cbttc/tester/javaTmp/$sub_argv/ApplyRecalibration -jar $ReadCfg{GATK} \\
	-T ApplyRecalibration \\
	-R $ReadCfg{$Method_REF} \\
	-input $Process_Path{$sub_argv}/runs/step_2-2_idv/$sub_argv.$runchr_group.idv.vcf \\
	-recalFile $Process_Path{$sub_argv}/runs/step_2-2_idv/$sub_argv.$runchr_group.idv.recal \\
	-tranchesFile $Process_Path{$sub_argv}/runs/step_2-2_idv/$sub_argv.$runchr_group.idv.tranches \\
	-L $ReadCfg{"$Method_REF.root"}/newDefine/$runchr_group.intervals \\
	-o $Process_Path{$sub_argv}/runs/step_2-2_idv/$sub_argv.$runchr_group.recalibrated.idv.vcf \\
	--ts_filter_level 99.0 \\
	-mode INDEL --excludeFiltered || { while [ -e $Analysis_Log{$sub_argv}.tmp ]; do sleep 1; done; perl $Bin/UpdateLog.pl $Analysis_Log{$sub_argv} "vqsr idv $runchr_group" 3 > $Analysis_Log{$sub_argv}.tmp && mv $Analysis_Log{$sub_argv}.tmp $Analysis_Log{$sub_argv}; exit 1; }
	
while [ -e $Analysis_Log{$sub_argv}.tmp ]; do sleep 1; done; perl $Bin/UpdateLog.pl $Analysis_Log{$sub_argv} "vqsr idv $runchr_group" 0 > $Analysis_Log{$sub_argv}.tmp && mv $Analysis_Log{$sub_argv}.tmp $Analysis_Log{$sub_argv}

while [ -e $Analysis_Ckp{$sub_argv}.tmp ]; do sleep 1; done; perl $Bin/UpdateCkp.pl $Analysis_Ckp{$sub_argv} idv.$runchr_group $Method_IDV 4 > $Analysis_Ckp{$sub_argv}.tmp && mv $Analysis_Ckp{$sub_argv}.tmp $Analysis_Ckp{$sub_argv}

perl $Bin/JobProcessor.pl $OUTDIR/$sub_argv $Analysis_Ckp{$sub_argv} idv.$runchr_group $Method_IDV
SCRIPTS
				my $initial_vqsr_off = <<"SCRIPTS";
#!/bin/sh
## BaseRecalibrator for idv $sub_argv
#\$ -cwd -l virtual_free=30G,h_vmem=50G

$ReadCfg{Java} -Xms8G -Xmx10G -Djava.io.tmpdir=/mnt/isilon/cbttc/tester/javaTmp/$sub_argv/SelectVariants -jar $ReadCfg{GATK} \\
	-T SelectVariants \\
	-R $ReadCfg{$Method_REF} \\
	-V $Process_Path{$sub_argv}/runs/step_2-1_snv/$sub_argv.$runchr_group.HaplotypeCaller.vcf \\
	-L $ReadCfg{"$Method_REF.root"}/newDefine/$runchr_group.intervals \\
	-selectType INDEL \\
	-o $Process_Path{$sub_argv}/runs/step_2-2_idv/$sub_argv.$runchr_group.idv.vcf || { while [ -e $Analysis_Log{$sub_argv}.tmp ]; do sleep 1; done; perl $Bin/UpdateLog.pl $Analysis_Log{$sub_argv} "vqsr idv $runchr_group" 3 > $Analysis_Log{$sub_argv}.tmp && mv $Analysis_Log{$sub_argv}.tmp $Analysis_Log{$sub_argv}; exit 1; }

ln $Process_Path{$sub_argv}/runs/step_2-2_idv/$sub_argv.$runchr_group.idv.vcf $Process_Path{$sub_argv}/runs/step_2-2_idv/$sub_argv.$runchr_group.recalibrated.idv.vcf
	
while [ -e $Analysis_Log{$sub_argv}.tmp ]; do sleep 1; done; perl $Bin/UpdateLog.pl $Analysis_Log{$sub_argv} "vqsr idv $runchr_group" 0 > $Analysis_Log{$sub_argv}.tmp && mv $Analysis_Log{$sub_argv}.tmp $Analysis_Log{$sub_argv}

while [ -e $Analysis_Ckp{$sub_argv}.tmp ]; do sleep 1; done; perl $Bin/UpdateCkp.pl $Analysis_Ckp{$sub_argv} idv.$runchr_group $Method_IDV 4 > $Analysis_Ckp{$sub_argv}.tmp && mv $Analysis_Ckp{$sub_argv}.tmp $Analysis_Ckp{$sub_argv}

perl $Bin/JobProcessor.pl $OUTDIR/$sub_argv $Analysis_Ckp{$sub_argv} idv.$runchr_group $Method_IDV
SCRIPTS
				if ( !defined $ReadCfg{'GATK.VQSR'} ) {$ReadCfg{'GATK.VQSR'} = 'off';}
				
				if ( $ReadCfg{'GATK.VQSR'} eq 'on' ){
					print FILE_HANDLE_VQSR $initial_vqsr_on;
				} elsif ( $ReadCfg{'GATK.VQSR'} eq 'off' ) {
					print FILE_HANDLE_VQSR $initial_vqsr_off;
				} else {
					die "bad vqsr setting in config file, only \'on\' & \'off\' allow to use, fix to run\n";
				}
				close FILE_HANDLE_VQSR;
				undef $initial_vqsr_off;
				undef $initial_vqsr_on;
				
				
				#------------------------------------------------------
				# idv annotation
				#------------------------------------------------------
				open FILE_HANDLE_ANN, ">$Process_Path{$sub_argv}/shells/step_2-2-4_ann_idv_$Method_IDV.$runchr_group.sh" or die $!;
				my $initial_ann;
				switch ( $Species ) {
					case 'human' {
						$initial_ann = <<"SCRIPTS";
#!/bin/sh
## annotation $Method_ANN for idv $sub_argv
#\$ -cwd -l virtual_free=20G,h_vmem=30G

$ReadCfg{Java} -Xmx10G -Xms8G -Djava.io.tmpdir=/mnt/isilon/cbttc/tester/javaTmp/$sub_argv/SelectVariants -jar $ReadCfg{$Method_IDV} \\
	-T SelectVariants \\
	-R $ReadCfg{$Method_REF} \\
	--discordance $ReadCfg{'mills.indel'} \\
	-V $Process_Path{$sub_argv}/runs/step_2-2_idv/$sub_argv.$runchr_group.recalibrated.idv.vcf \\
	-L $ReadCfg{"$Method_REF.root"}/newDefine/$runchr_group.intervals \\
	-o $Process_Path{$sub_argv}/runs/step_2-2_idv/$sub_argv.$runchr_group.rmidp1.idv.vcf || { while [ -e $Analysis_Log{$sub_argv}.tmp ]; do sleep 1; done; perl $Bin/UpdateLog.pl $Analysis_Log{$sub_argv} "ann idv $runchr_group" 3 > $Analysis_Log{$sub_argv}.tmp && mv $Analysis_Log{$sub_argv}.tmp $Analysis_Log{$sub_argv}; exit 1; }
	
$ReadCfg{Java} -Xmx10G -Xms8G -Djava.io.tmpdir=/mnt/isilon/cbttc/tester/javaTmp/$sub_argv/SelectVariants -jar $ReadCfg{$Method_IDV} \\
	-T SelectVariants \\
	-R $ReadCfg{$Method_REF} \\
	--discordance $ReadCfg{'1000g.indel'} \\
	-V $Process_Path{$sub_argv}/runs/step_2-2_idv/$sub_argv.$runchr_group.rmidp1.idv.vcf \\
	-L $ReadCfg{"$Method_REF.root"}/newDefine/$runchr_group.intervals \\
	-o $Process_Path{$sub_argv}/runs/step_2-2_idv/$sub_argv.$runchr_group.rmidp.idv.vcf || { while [ -e $Analysis_Log{$sub_argv}.tmp ]; do sleep 1; done; perl $Bin/UpdateLog.pl $Analysis_Log{$sub_argv} "ann idv $runchr_group" 3 > $Analysis_Log{$sub_argv}.tmp && mv $Analysis_Log{$sub_argv}.tmp $Analysis_Log{$sub_argv}; exit 1; }

$ReadCfg{Java} -Xms4G -Xmx8G -Djava.io.tmpdir=/mnt/isilon/cbttc/tester/javaTmp/$sub_argv/snpEff -jar $ReadCfg{snpEff} \\
	-c $ReadCfg{snpEffCfg} \\
	$ReadCfg{'snpEff.human'} \\
	$Process_Path{$sub_argv}/runs/step_2-2_idv/$sub_argv.$runchr_group.rmidp.idv.vcf > $Process_Path{$sub_argv}/runs/step_2-2_idv/$sub_argv.idv.snpeff.vcf || { while [ -e $Analysis_Log{$sub_argv}.tmp ]; do sleep 1; done; perl $Bin/UpdateLog.pl $Analysis_Log{$sub_argv} "ann idv $runchr_group" 3 > $Analysis_Log{$sub_argv}.tmp && mv $Analysis_Log{$sub_argv}.tmp $Analysis_Log{$sub_argv}; exit 1; }

while [ -e $Analysis_Log{$sub_argv}.tmp ]; do sleep 1; done; perl $Bin/UpdateLog.pl $Analysis_Log{$sub_argv} "ann idv $runchr_group" 0 > $Analysis_Log{$sub_argv}.tmp && mv $Analysis_Log{$sub_argv}.tmp $Analysis_Log{$sub_argv}

while [ -e $Analysis_Ckp{$sub_argv}.tmp ]; do sleep 1; done; perl $Bin/UpdateCkp.pl $Analysis_Ckp{$sub_argv} idv.$runchr_group $Method_IDV 5 > $Analysis_Ckp{$sub_argv}.tmp && mv $Analysis_Ckp{$sub_argv}.tmp $Analysis_Ckp{$sub_argv}

perl $Bin/JobProcessor.pl $OUTDIR/$sub_argv $Analysis_Ckp{$sub_argv} idv.$runchr_group $Method_IDV
SCRIPTS
					}
					case 'mouse' {
						$initial_ann = <<"SCRIPTS";
#!/bin/sh
## annotation $Method_ANN for idv $sub_argv
#\$ -cwd -l virtual_free=20G,h_vmem=30G

$ReadCfg{Java} -Xmx10G -Xms8G -Djava.io.tmpdir=/mnt/isilon/cbttc/tester/javaTmp/$sub_argv/SelectVariants -jar $ReadCfg{$Method_IDV} \\
	-T SelectVariants \\
	-R $ReadCfg{$Method_REF} \\
	--discordance $ReadCfg{'sanger.indel'} \\
	-V $Process_Path{$sub_argv}/runs/step_2-2_idv/$sub_argv.$runchr_group.recalibrated.idv.vcf \\
	-L $ReadCfg{"$Method_REF.root"}/newDefine/$runchr_group.intervals \\
	-o $Process_Path{$sub_argv}/runs/step_2-2_idv/$sub_argv.$runchr_group.rmidp.idv.vcf || { while [ -e $Analysis_Log{$sub_argv}.tmp ]; do sleep 1; done; perl $Bin/UpdateLog.pl $Analysis_Log{$sub_argv} "ann idv $runchr_group" 3 > $Analysis_Log{$sub_argv}.tmp && mv $Analysis_Log{$sub_argv}.tmp $Analysis_Log{$sub_argv}; exit 1; }

$ReadCfg{Java} -Xms4G -Xmx8G -Djava.io.tmpdir=/mnt/isilon/cbttc/tester/javaTmp/$sub_argv/snpEff -jar $ReadCfg{snpEff} \\
	-c $ReadCfg{snpEffCfg} \\
	$ReadCfg{'snpEff.mouse'} \\
	$Process_Path{$sub_argv}/runs/step_2-2_idv/$sub_argv.$runchr_group.rmidp.idv.vcf > $Process_Path{$sub_argv}/runs/step_2-2_idv/$sub_argv.idv.snpeff.vcf || { while [ -e $Analysis_Log{$sub_argv}.tmp ]; do sleep 1; done; perl $Bin/UpdateLog.pl $Analysis_Log{$sub_argv} "ann idv $runchr_group" 3 > $Analysis_Log{$sub_argv}.tmp && mv $Analysis_Log{$sub_argv}.tmp $Analysis_Log{$sub_argv}; exit 1; }

while [ -e $Analysis_Log{$sub_argv}.tmp ]; do sleep 1; done; perl $Bin/UpdateLog.pl $Analysis_Log{$sub_argv} "ann idv $runchr_group" 0 > $Analysis_Log{$sub_argv}.tmp && mv $Analysis_Log{$sub_argv}.tmp $Analysis_Log{$sub_argv}

while [ -e $Analysis_Ckp{$sub_argv}.tmp ]; do sleep 1; done; perl $Bin/UpdateCkp.pl $Analysis_Ckp{$sub_argv} idv.$runchr_group $Method_IDV 5 > $Analysis_Ckp{$sub_argv}.tmp && mv $Analysis_Ckp{$sub_argv}.tmp $Analysis_Ckp{$sub_argv}

perl $Bin/JobProcessor.pl $OUTDIR/$sub_argv $Analysis_Ckp{$sub_argv} idv.$runchr_group $Method_IDV
SCRIPTS
					}
				}
				print FILE_HANDLE_ANN $initial_ann;
				close FILE_HANDLE_ANN;
				undef $initial_ann;
			}
			
			case ('VarScan') {
				die "not cancer mode, bad snv method setting!" unless $CancerMode == 1;
				my $normal_sample = $NormalSample{$PairMarker{$sub_argv}};
				my $tumor_sample = $TumorSample{$PairMarker{$sub_argv}};
				#------------------------------------------------------
				# mpileup && check cancer pair status 
				#------------------------------------------------------
				open FILE_HANDLE_CALL, ">$Process_Path{$sub_argv}/shells/step_2-2-1_call_idv_$Method_IDV.$runchr_group.sh" or die $!;
				my $initial_call = <<"SCRIPTS";
#!/bin/sh
## call somatic InDel $PairMarker{$sub_argv}
#\$ -cwd -l virtual_free=1G,h_vmem=1G

$ReadCfg{Samtools} mpileup -f $ReadCfg{$Method_REF} $Process_Path{$sub_argv}/runs/step_1-1_map/$sub_argv.$runchr_group.bam > $Process_Path{$sub_argv}/runs/step_1-1_map/$sub_argv.$runchr_group.pileup || { echo "mpileup $sub_argv failed"; while [ -e $Analysis_Log{$sub_argv}.tmp ]; do sleep 1; done; perl $Bin/UpdateLog.pl $Analysis_Log{$sub_argv} "call idv $runchr_group" 3 > $Analysis_Log{$sub_argv}.tmp && mv $Analysis_Log{$sub_argv}.tmp $Analysis_Log{$sub_argv}; exit 1; }

while [ -e $Analysis_Ckp{$sub_argv}.tmp ]; do sleep 1; done; perl $Bin/UpdateCkp.pl $Analysis_Ckp{$sub_argv} idv.$runchr_group $Method_IDV 3 > $Analysis_Ckp{$sub_argv}.tmp && mv $Analysis_Ckp{$sub_argv}.tmp $Analysis_Ckp{$sub_argv}

# check pair status
perl $Bin/CheckPair.pl $Analysis_Ckp{$normal_sample} $Analysis_Ckp{$tumor_sample} idv.$runchr_group 3 || { echo "Pair mate $PairLink{$sub_argv} is not ready" && exit 1; }

while [ -e $Analysis_Log{$normal_sample}}.tmp ]; do sleep 1; done; perl $Bin/UpdateLog.pl $Analysis_Log{$normal_sample} "call idv $runchr_group" 0 > $Analysis_Log{$normal_sample}.tmp && mv $Analysis_Log{$normal_sample}.tmp $Analysis_Log{$normal_sample}
while [ -e $Analysis_Log{$tumor_sample}.tmp ]; do sleep 1; done; perl $Bin/UpdateLog.pl $Analysis_Log{$tumor_sample} "call idv $runchr_group" 0 > $Analysis_Log{$tumor_sample}.tmp && mv $Analysis_Log{$tumor_sample}.tmp $Analysis_Log{$tumor_sample}

# next step
perl $Bin/JobProcessor.pl $OUTDIR/$normal_sample $Analysis_Ckp{$normal_sample} idv.$runchr_group $Method_IDV

SCRIPTS
				print FILE_HANDLE_CALL $initial_call;
				close FILE_HANDLE_CALL;
				undef $initial_call;
				
				
				#------------------------------------------------------
				# VarScan for both samples
				#------------------------------------------------------
				open FILE_HANDLE_SOMATIC, ">$Process_Path{$sub_argv}/shells/step_2-2-3_somatic_idv_$Method_IDV.$runchr_group.sh" or die $!;
				my $initial_somatic = <<"SCRIPTS";
#!/bin/sh
## call somatic InDel $PairMarker{$sub_argv}
#\$ -cwd -l virtual_free=20G,h_vmem=35G -pe smp 8

$ReadCfg{Java} -Xms5G -Xmx15G -Djava.io.tmpdir=/mnt/isilon/cbttc/tester/javaTmp/$sub_argv/VarScan -jar $ReadCfg{$Method_IDV} \\
	somatic $Process_Path{$normal_sample}/runs/step_1-1_map/$normal_sample.$runchr_group.pileup \\
	$Process_Path{$tumor_sample}/runs/step_1-1_map/$tumor_sample.$runchr_group.pileup \\
	$Process_Path{$PairMarker{$sub_argv}}/step_2-2_idv/$runchr_group.VarScan2 \\
	$ReadCfg{VarScanOpts} --output-vcf || { echo "varscan $sub_argv failed"; while [ -e $Analysis_Log{$sub_argv}.tmp ]; do sleep 1; done; perl $Bin/UpdateLog.pl $Analysis_Log{$sub_argv} "somatic idv $runchr_group" 3 > $Analysis_Log{$sub_argv}.tmp && mv $Analysis_Log{$sub_argv}.tmp $Analysis_Log{$sub_argv}; exit 1; }

$Rm_Mark $Process_Path{$normal_sample}/runs/step_1-1_map/$normal_sample.$runchr_group.pileup $Process_Path{$tumor_sample}/runs/step_1-1_map/$tumor_sample.$runchr_group.pileup 
	
mv $Process_Path{$PairMarker{$sub_argv}}/step_2-2_idv/$runchr_group.VarScan2.snp.vcf $Process_Path{$PairMarker{$sub_argv}}/step_2-1_snv/$runchr_group.VarScan2.snv.vcf
perl $Bin/FormatVarScan.pl $Process_Path{$PairMarker{$sub_argv}}/step_2-2_idv/$runchr_group.VarScan2.indel.vcf > $Process_Path{$PairMarker{$sub_argv}}/step_2-2_idv/$runchr_group.VarScan2.idv.vcf

while [ -e $Analysis_Log{$normal_sample}.tmp ]; do sleep 1; done; perl $Bin/UpdateLog.pl $Analysis_Log{$normal_sample} "somatic idv $runchr_group" 0 > $Analysis_Log{$normal_sample}.tmp && mv $Analysis_Log{$normal_sample}.tmp $Analysis_Log{$normal_sample}

while [ -e $Analysis_Ckp{$normal_sample}.tmp ]; do sleep 1; done; perl $Bin/UpdateCkp.pl $Analysis_Ckp{$normal_sample} idv.$runchr_group $Method_IDV 4 > $Analysis_Ckp{$normal_sample}.tmp && mv $Analysis_Ckp{$normal_sample}.tmp $Analysis_Ckp{$normal_sample}

perl $Bin/JobProcessor.pl $OUTDIR/$normal_sample $Analysis_Ckp{$normal_sample} idv.$runchr_group $Method_IDV
	
SCRIPTS
				print FILE_HANDLE_SOMATIC $initial_somatic;
				close FILE_HANDLE_SOMATIC;
				undef $initial_somatic;
				
				
				#------------------------------------------------------
				# idv annotation
				#------------------------------------------------------
				open FILE_HANDLE_ANN, ">$Process_Path{$sub_argv}/shells/step_2-2-4_ann_idv_$Method_IDV.$runchr_group.sh" or die $!;
				my $initial_ann;
				switch ( $Species ) {
					case 'human' {
						$initial_ann = <<"SCRIPTS";
#!/bin/sh
## annotation $Method_ANN for idv $PairMarker{$sub_argv}
#\$ -cwd -l virtual_free=20G,h_vmem=30G -pe smp 8

$ReadCfg{Java} -Xmx10G -Xms8G -Djava.io.tmpdir=/mnt/isilon/cbttc/tester/javaTmp/$sub_argv/SelectVariants -jar $ReadCfg{GATK} \\
	-T SelectVariants \\
	-R $ReadCfg{$Method_REF} \\
	--discordance $ReadCfg{'mills.indel'} \\
	-V $Process_Path{$PairMarker{$sub_argv}}/step_2-2_idv/$runchr_group.VarScan2.idv.vcf \\
	-L $ReadCfg{"$Method_REF.root"}/newDefine/$runchr_group.intervals \\
	-o $Process_Path{$PairMarker{$sub_argv}}/step_2-2_idv/$runchr_group.rmidp1.idv.vcf || { while [ -e $Analysis_Log{$normal_sample}.tmp ]; do sleep 1; done; perl $Bin/UpdateLog.pl $Analysis_Log{$normal_sample} "ann idv $runchr_group" 3 > $Analysis_Log{$normal_sample}.tmp && mv $Analysis_Log{$normal_sample}.tmp $Analysis_Log{$normal_sample}; exit 1; }
	
$ReadCfg{Java} -Xmx10G -Xms8G -Djava.io.tmpdir=/mnt/isilon/cbttc/tester//$sub_argv/SelectVariants -jar $ReadCfg{GATK} \\
	-T SelectVariants \\
	-R $ReadCfg{$Method_REF} \\
	--discordance $ReadCfg{'1000g.indel'} \\
	-V $Process_Path{$PairMarker{$sub_argv}}/step_2-2_idv/$runchr_group.rmidp1.idv.vcf \\
	-L $ReadCfg{"$Method_REF.root"}/newDefine/$runchr_group.intervals \\
	-o $Process_Path{$PairMarker{$sub_argv}}/step_2-2_idv/$runchr_group.rmidp.idv.vcf || { while [ -e $Analysis_Log{$normal_sample}.tmp ]; do sleep 1; done; perl $Bin/UpdateLog.pl $Analysis_Log{$normal_sample} "ann idv $runchr_group" 3 > $Analysis_Log{$normal_sample}.tmp && mv $Analysis_Log{$normal_sample}.tmp $Analysis_Log{$normal_sample}; exit 1; }

$ReadCfg{Java} -Xms4G -Xmx8G -Djava.io.tmpdir=/mnt/isilon/cbttc/tester/javaTmp/$sub_argv/snpEff -jar $ReadCfg{snpEff} \\
	-c $ReadCfg{snpEffCfg} \\
	$ReadCfg{'snpEff.human'} \\
	$Process_Path{$PairMarker{$sub_argv}}/step_2-2_idv/$runchr_group.rmidp.idv.vcf > $Process_Path{$PairMarker{$sub_argv}}/step_2-2_idv/$runchr_group.idv.snpeff.vcf || { while [ -e $Analysis_Log{$normal_sample}.tmp ]; do sleep 1; done; perl $Bin/UpdateLog.pl $Analysis_Log{$normal_sample} "ann idv $runchr_group" 3 > $Analysis_Log{$normal_sample}.tmp && mv $Analysis_Log{$normal_sample}.tmp $Analysis_Log{$normal_sample}; exit 1; }

while [ -e $Analysis_Log{$normal_sample}.tmp ]; do sleep 1; done; perl $Bin/UpdateLog.pl $Analysis_Log{$normal_sample} "ann idv $runchr_group" 0 > $Analysis_Log{$normal_sample}.tmp && mv $Analysis_Log{$normal_sample}.tmp $Analysis_Log{$normal_sample}

while [ -e $Analysis_Ckp{$normal_sample}.tmp ]; do sleep 1; done; perl $Bin/UpdateCkp.pl $Analysis_Ckp{$normal_sample} idv.$runchr_group $Method_IDV 5 > $Analysis_Ckp{$normal_sample}.tmp && mv $Analysis_Ckp{$normal_sample}.tmp $Analysis_Ckp{$normal_sample}

# perl $Bin/JobProcessor.pl $OUTDIR/$normal_sample $Analysis_Ckp{$normal_sample} idv.$runchr_group $Method_IDV
SCRIPTS
			}
			
					case 'mouse' {
						$initial_ann = <<"SCRIPTS";
#!/bin/sh
## annotation $Method_ANN for idv $PairMarker{$sub_argv}
#\$ -cwd -l virtual_free=20G,h_vmem=30G -pe smp 8

$ReadCfg{Java} -Xms4G -Xmx8G -Djava.io.tmpdir=/mnt/isilon/cbttc/tester/javaTmp/$sub_argv/snpEff -jar $ReadCfg{snpEff} \\
	-c $ReadCfg{snpEffCfg} \\
	$ReadCfg{'snpEff.mouse'} \\
	$Process_Path{$PairMarker{$sub_argv}}/step_2-1_snv/$runchr_group.VarScan2.snv.vcf > $Process_Path{$PairMarker{$sub_argv}}/step_2-1_snv/$runchr_group.VarScan2.snv.snpeff.vcf || { while [ -e $Analysis_Log{$normal_sample}.tmp ]; do sleep 1; done; perl $Bin/UpdateLog.pl $Analysis_Log{$normal_sample} "ann idv $runchr_group" 3 > $Analysis_Log{$normal_sample}.tmp && mv $Analysis_Log{$normal_sample}.tmp $Analysis_Log{$normal_sample}; exit 1; }

$ReadCfg{Java} -Xms4G -Xmx8G -Djava.io.tmpdir=/mnt/isilon/cbttc/tester/javaTmp/$sub_argv/snpEff -jar $ReadCfg{snpEff} \\
	-c $ReadCfg{snpEffCfg} \\
	$ReadCfg{'snpEff.mouse'} \\
	$Process_Path{$PairMarker{$sub_argv}}/step_2-2_idv/$runchr_group.VarScan2.idv.vcf > $Process_Path{$PairMarker{$sub_argv}}/step_2-2_idv/$runchr_group.VarScan2.idv.snpeff.vcf || { while [ -e $Analysis_Log{$normal_sample}.tmp ]; do sleep 1; done; perl $Bin/UpdateLog.pl $Analysis_Log{$normal_sample} "ann idv $runchr_group" 3 > $Analysis_Log{$normal_sample}.tmp && mv $Analysis_Log{$normal_sample}.tmp $Analysis_Log{$normal_sample}; exit 1; }

while [ -e $Analysis_Log{$normal_sample}.tmp ]; do sleep 1; done; perl $Bin/UpdateLog.pl $Analysis_Log{$normal_sample} "ann idv $runchr_group" 0 > $Analysis_Log{$normal_sample}.tmp && mv $Analysis_Log{$normal_sample}.tmp $Analysis_Log{$normal_sample}

while [ -e $Analysis_Ckp{$normal_sample}.tmp ]; do sleep 1; done; perl $Bin/UpdateCkp.pl $Analysis_Ckp{$normal_sample} idv.$runchr_group $Method_IDV 5 > $Analysis_Ckp{$normal_sample}.tmp && mv $Analysis_Ckp{$normal_sample}.tmp $Analysis_Ckp{$normal_sample}

# perl $Bin/JobProcessor.pl $OUTDIR/$normal_sample $Analysis_Ckp{$normal_sample} idv.$runchr_group $Method_IDV
SCRIPTS
					}
				}
				print FILE_HANDLE_ANN $initial_ann;
				close FILE_HANDLE_ANN;
				undef $initial_ann;
			}
		
			case ('Strelka') {
				die "not cancer mode, bad snv method setting!" unless $CancerMode == 1;
				my $normal_sample = $NormalSample{$PairMarker{$sub_argv}};
				my $tumor_sample = $TumorSample{$PairMarker{$sub_argv}};
				#------------------------------------------------------
				# pre-Strelka
				#------------------------------------------------------
				open FILE_HANDLE_CALL, ">$Process_Path{$sub_argv}/shells/step_2-2-1_call_idv_$Method_IDV.$runchr_group.sh" or die $!;
				my $initial_call = <<"SCRIPTS";
#!/bin/sh
## call somatic InDel $PairMarker{$sub_argv}
#\$ -cwd -l virtual_free=0.1G,h_vmem=0.1G
echo ready
SCRIPTS
				print FILE_HANDLE_CALL $initial_call;
				close FILE_HANDLE_CALL;
				undef $initial_call;

				#------------------------------------------------------
				# Strelka for both samples
				#------------------------------------------------------
				open FILE_HANDLE_SOMATIC, ">$Process_Path{$sub_argv}/shells/step_2-2-3_somatic_idv_$Method_IDV.$runchr_group.sh" or die $!;
				my $initial_somatic = <<"SCRIPTS";
#!/bin/sh
## call somatic InDel $PairMarker{$sub_argv}
#\$ -cwd -l virtual_free=10G,h_vmem=20G -pe smp 3

source /home/zhuy/.bash_profile
perl $ReadCfg{Strelka} --tumor=$Process_Path{$tumor_sample}/runs/step_1-1_map/$tumor_sample.$runchr_group.bam --normal=$Process_Path{$normal_sample}/runs/step_1-1_map/$normal_sample.$runchr_group.bam --ref=$ReadCfg{$Method_REF} --config=$ReadCfg{StrelkaCfg} --output-dir=$Process_Path{$PairMarker{$sub_argv}}/step_2-2_idv/$runchr_group.Strelka
make -C $Process_Path{$PairMarker{$sub_argv}}/step_2-2_idv/$runchr_group.Strelka || { echo "Strelka $sub_argv failed"; while [ -e $Analysis_Log{$sub_argv}.tmp ]; do sleep 1; done; perl $Bin/UpdateLog.pl $Analysis_Log{$sub_argv} "somatic idv $runchr_group" 3 > $Analysis_Log{$sub_argv}.tmp && mv $Analysis_Log{$sub_argv}.tmp $Analysis_Log{$sub_argv}; exit 1; }
	
while [ -e $Analysis_Log{$normal_sample}.tmp ]; do sleep 1; done; perl $Bin/UpdateLog.pl $Analysis_Log{$normal_sample} "somatic idv $runchr_group" 0 > $Analysis_Log{$normal_sample}.tmp && mv $Analysis_Log{$normal_sample}.tmp $Analysis_Log{$normal_sample}

while [ -e $Analysis_Ckp{$normal_sample}.tmp ]; do sleep 1; done; perl $Bin/UpdateCkp.pl $Analysis_Ckp{$normal_sample} idv.$runchr_group $Method_IDV 4 > $Analysis_Ckp{$normal_sample}.tmp && mv $Analysis_Ckp{$normal_sample}.tmp $Analysis_Ckp{$normal_sample}

perl $Bin/JobProcessor.pl $OUTDIR/$normal_sample $Analysis_Ckp{$normal_sample} idv.$runchr_group $Method_IDV
	
SCRIPTS
				print FILE_HANDLE_SOMATIC $initial_somatic;
				close FILE_HANDLE_SOMATIC;
				undef $initial_somatic;

				#------------------------------------------------------
				# idv annotation
				#------------------------------------------------------
				open FILE_HANDLE_ANN, ">$Process_Path{$sub_argv}/shells/step_2-2-4_ann_idv_$Method_IDV.$runchr_group.sh" or die $!;
				my $initial_ann;
				switch ( $Species ) {
					case 'human' {
						$initial_ann = <<"SCRIPTS";
#!/bin/sh
## annotation $Method_ANN for idv $PairMarker{$sub_argv}
#\$ -cwd -l virtual_free=10G,h_vmem=15G -pe smp 3

$ReadCfg{Java} -Xms4G -Xmx8G -Djava.io.tmpdir=/mnt/isilon/cbttc/tester/javaTmp/$sub_argv/snpEff -jar $ReadCfg{snpEff} \\
	-c $ReadCfg{snpEffCfg} \\
	$ReadCfg{'snpEff.human'} \\
	$Process_Path{$PairMarker{$sub_argv}}/step_2-2_idv/$runchr_group.Strelka/results/passed.somatic.snvs.vcf > $Process_Path{$PairMarker{$sub_argv}}/step_2-2_idv/$runchr_group.snv.snpeff.vcf

$ReadCfg{Java} -Xms4G -Xmx8G -Djava.io.tmpdir=/mnt/isilon/cbttc/tester/javaTmp/$sub_argv/snpEff -jar $ReadCfg{snpEff} \\
	-c $ReadCfg{snpEffCfg} \\
	$ReadCfg{'snpEff.human'} \\
	$Process_Path{$PairMarker{$sub_argv}}/step_2-2_idv/$runchr_group.Strelka/results/passed.somatic.indels.vcf > $Process_Path{$PairMarker{$sub_argv}}/step_2-2_idv/$runchr_group.idv.snpeff.vcf || { while [ -e $Analysis_Log{$normal_sample}.tmp ]; do sleep 1; done; perl $Bin/UpdateLog.pl $Analysis_Log{$normal_sample} "ann idv $runchr_group" 3 > $Analysis_Log{$normal_sample}.tmp && mv $Analysis_Log{$normal_sample}.tmp $Analysis_Log{$normal_sample}; exit 1; }

while [ -e $Analysis_Log{$normal_sample}.tmp ]; do sleep 1; done; perl $Bin/UpdateLog.pl $Analysis_Log{$normal_sample} "ann idv $runchr_group" 0 > $Analysis_Log{$normal_sample}.tmp && mv $Analysis_Log{$normal_sample}.tmp $Analysis_Log{$normal_sample}

while [ -e $Analysis_Ckp{$normal_sample}.tmp ]; do sleep 1; done; perl $Bin/UpdateCkp.pl $Analysis_Ckp{$normal_sample} idv.$runchr_group $Method_IDV 5 > $Analysis_Ckp{$normal_sample}.tmp && mv $Analysis_Ckp{$normal_sample}.tmp $Analysis_Ckp{$normal_sample}

# perl $Bin/JobProcessor.pl $OUTDIR/$normal_sample $Analysis_Ckp{$normal_sample} idv.$runchr_group $Method_IDV
SCRIPTS
					}
				}
				print FILE_HANDLE_ANN $initial_ann;
				close FILE_HANDLE_ANN;
				undef $initial_ann;
			}
		}
	}	
	undef $sub_argv;
}

sub _cnvScripts {
	my $sub_argv = shift;

	foreach my $runchr_group ( @RunChr_Group ) {
		#--------------------------------------------------------------
		# initial sv scripts
		#--------------------------------------------------------------
		switch ( $Method_CNV ) {
			case ('CNVnator') {
				#------------------------------------------------------
				# CNVnator tree
				#------------------------------------------------------
				open FILE_HANDLE_CALL, ">$Process_Path{$sub_argv}/shells/step_2-3-1_call_cnv_$Method_CNV.$runchr_group\_tree.sh" or die $!;
				my $initial_cnv_call = <<"SCRIPTS";
#!/bin/sh
## cnv for $sub_argv
#\$ -cwd -l virtual_free=10G,h_vmem=35G -pe smp 8
export ROOTSYS=/mnt/isilon/cbttc/tools/root
export PATH=\$PATH:\$ROOTSYS/bin
export SHLIB_PATH=\$SHLIB_PATH:\$ROOTSYS/lib
export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:\$ROOTSYS/lib:/cm/shared/apps/gcc/4.7.0/lib64
export LIBPATH=\$LIBPATH:\$ROOTSYS/lib
{ $ReadCfg{CNVnator} -root $Process_Path{$sub_argv}/runs/step_2-3_cnv/$sub_argv.$runchr_group.root -chrom $Group_Content{$runchr_group} -tree $Process_Path{$sub_argv}/runs/step_1-1_map/$sub_argv.$runchr_group.bam -unique; } || { while [ -e $Analysis_Log{$sub_argv}.tmp ]; do sleep 1; done; perl $Bin/UpdateLog.pl $Analysis_Log{$sub_argv} "call cnv tree $runchr_group" 3 > $Analysis_Log{$sub_argv}.tmp && mv $Analysis_Log{$sub_argv}.tmp $Analysis_Log{$sub_argv}; exit 1; }
	
while [ -e $Analysis_Log{$sub_argv}.tmp ]; do sleep 1; done; perl $Bin/UpdateLog.pl $Analysis_Log{$sub_argv} "call cnv tree $runchr_group" 0 > $Analysis_Log{$sub_argv}.tmp && mv $Analysis_Log{$sub_argv}.tmp $Analysis_Log{$sub_argv}

while [ -e $Analysis_Ckp{$sub_argv}.tmp ]; do sleep 1; done; perl $Bin/UpdateCkp.pl $Analysis_Ckp{$sub_argv} cnv.$runchr_group $Method_CNV 2 > $Analysis_Ckp{$sub_argv}.tmp && mv $Analysis_Ckp{$sub_argv}.tmp $Analysis_Ckp{$sub_argv}

perl $Bin/JobProcessor.pl $OUTDIR/$sub_argv $Analysis_Ckp{$sub_argv} cnv.$runchr_group $Method_CNV
SCRIPTS
				print FILE_HANDLE_CALL $initial_cnv_call;
				close FILE_HANDLE_CALL;
				undef $initial_cnv_call;
			
			
				#------------------------------------------------------
				# CNVnator histogram
				#------------------------------------------------------
				open FILE_HANDLE_CALL, ">$Process_Path{$sub_argv}/shells/step_2-3-2_call_cnv_$Method_CNV.$runchr_group\_histogram.sh" or die $!;
				$initial_cnv_call = <<"SCRIPTS";
#!/bin/sh
## cnv for $sub_argv
#\$ -cwd -l virtual_free=2G,h_vmem=10G -pe smp 1
export ROOTSYS=/mnt/isilon/cbttc/tools/root
export PATH=\$PATH:\$ROOTSYS/bin
export SHLIB_PATH=\$SHLIB_PATH:\$ROOTSYS/lib
export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:\$ROOTSYS/lib:/cm/shared/apps/gcc/4.7.0/lib64
export LIBPATH=\$LIBPATH:\$ROOTSYS/lib
{ $ReadCfg{CNVnator} -root $Process_Path{$sub_argv}/runs/step_2-3_cnv/$sub_argv.$runchr_group.root -his $ReadLen{$sub_argv} -d $ReadCfg{"$Method_REF.root"}/original; } || { while [ -e $Analysis_Log{$sub_argv}.tmp ]; do sleep 1; done; perl $Bin/UpdateLog.pl $Analysis_Log{$sub_argv} "call cnv histogram" 3 > $Analysis_Log{$sub_argv}.tmp && mv $Analysis_Log{$sub_argv}.tmp $Analysis_Log{$sub_argv}; exit 1; }
	
while [ -e $Analysis_Log{$sub_argv}.tmp ]; do sleep 1; done; perl $Bin/UpdateLog.pl $Analysis_Log{$sub_argv} "call cnv histogram $runchr_group" 0 > $Analysis_Log{$sub_argv}.tmp && mv $Analysis_Log{$sub_argv}.tmp $Analysis_Log{$sub_argv}

while [ -e $Analysis_Ckp{$sub_argv}.tmp ]; do sleep 1; done; perl $Bin/UpdateCkp.pl $Analysis_Ckp{$sub_argv} cnv.$runchr_group $Method_CNV 3 > $Analysis_Ckp{$sub_argv}.tmp && mv $Analysis_Ckp{$sub_argv}.tmp $Analysis_Ckp{$sub_argv}

perl $Bin/JobProcessor.pl $OUTDIR/$sub_argv $Analysis_Ckp{$sub_argv} cnv.$runchr_group $Method_CNV
SCRIPTS
				print FILE_HANDLE_CALL $initial_cnv_call;
				close FILE_HANDLE_CALL;
				undef $initial_cnv_call;
				
				
				#------------------------------------------------------
				# CNVnator stat
				#------------------------------------------------------
				open FILE_HANDLE_CALL, ">$Process_Path{$sub_argv}/shells/step_2-3-3_call_cnv_$Method_CNV.$runchr_group\_stat.sh" or die $!;
				$initial_cnv_call = <<"SCRIPTS";
#!/bin/sh
## cnv for $sub_argv
#\$ -cwd -l virtual_free=2G,h_vmem=20G -pe smp 1
export ROOTSYS=/mnt/isilon/cbttc/tools/root
export PATH=\$PATH:\$ROOTSYS/bin
export SHLIB_PATH=\$SHLIB_PATH:\$ROOTSYS/lib
export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:\$ROOTSYS/lib:/cm/shared/apps/gcc/4.7.0/lib64
export LIBPATH=\$LIBPATH:\$ROOTSYS/lib
{ $ReadCfg{CNVnator} -root $Process_Path{$sub_argv}/runs/step_2-3_cnv/$sub_argv.$runchr_group.root -stat $ReadLen{$sub_argv}; } || { while [ -e $Analysis_Log{$sub_argv}.tmp ]; do sleep 1; done; perl $Bin/UpdateLog.pl $Analysis_Log{$sub_argv} "call cnv stat $runchr_group" 3 > $Analysis_Log{$sub_argv}.tmp && mv $Analysis_Log{$sub_argv}.tmp $Analysis_Log{$sub_argv}; exit 1; }
	
while [ -e $Analysis_Log{$sub_argv}.tmp ]; do sleep 1; done; perl $Bin/UpdateLog.pl $Analysis_Log{$sub_argv} "call cnv stat $runchr_group" 0 > $Analysis_Log{$sub_argv}.tmp && mv $Analysis_Log{$sub_argv}.tmp $Analysis_Log{$sub_argv}

while [ -e $Analysis_Ckp{$sub_argv}.tmp ]; do sleep 1; done; perl $Bin/UpdateCkp.pl $Analysis_Ckp{$sub_argv} cnv.$runchr_group $Method_CNV 4 > $Analysis_Ckp{$sub_argv}.tmp && mv $Analysis_Ckp{$sub_argv}.tmp $Analysis_Ckp{$sub_argv}

perl $Bin/JobProcessor.pl $OUTDIR/$sub_argv $Analysis_Ckp{$sub_argv} cnv.$runchr_group $Method_CNV
SCRIPTS
				print FILE_HANDLE_CALL $initial_cnv_call;
				close FILE_HANDLE_CALL;
				undef $initial_cnv_call;
				
				
				#------------------------------------------------------
				# CNVnator partition
				#------------------------------------------------------
				open FILE_HANDLE_CALL, ">$Process_Path{$sub_argv}/shells/step_2-3-4_call_cnv_$Method_CNV.$runchr_group\_partition.sh" or die $!;
				$initial_cnv_call = <<"SCRIPTS";
#!/bin/sh
## cnv for $sub_argv
#\$ -cwd -l virtual_free=2G,h_vmem=20G
export ROOTSYS=/mnt/isilon/cbttc/tools/root
export PATH=\$PATH:\$ROOTSYS/bin
export SHLIB_PATH=\$SHLIB_PATH:\$ROOTSYS/lib
export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:\$ROOTSYS/lib:/cm/shared/apps/gcc/4.7.0/lib64
export LIBPATH=\$LIBPATH:\$ROOTSYS/lib
{ $ReadCfg{CNVnator} -root $Process_Path{$sub_argv}/runs/step_2-3_cnv/$sub_argv.$runchr_group.root -partition $ReadLen{$sub_argv}; } || { while [ -e $Analysis_Log{$sub_argv}.tmp ]; do sleep 1; done; perl $Bin/UpdateLog.pl $Analysis_Log{$sub_argv} "call cnv partition $runchr_group" 3 > $Analysis_Log{$sub_argv}.tmp && mv $Analysis_Log{$sub_argv}.tmp $Analysis_Log{$sub_argv}; exit 1; }
	
while [ -e $Analysis_Log{$sub_argv}.tmp ]; do sleep 1; done; perl $Bin/UpdateLog.pl $Analysis_Log{$sub_argv} "call cnv partition $runchr_group" 0 > $Analysis_Log{$sub_argv}.tmp && mv $Analysis_Log{$sub_argv}.tmp $Analysis_Log{$sub_argv}

while [ -e $Analysis_Ckp{$sub_argv}.tmp ]; do sleep 1; done; perl $Bin/UpdateCkp.pl $Analysis_Ckp{$sub_argv} cnv.$runchr_group $Method_CNV 5 > $Analysis_Ckp{$sub_argv}.tmp && mv $Analysis_Ckp{$sub_argv}.tmp $Analysis_Ckp{$sub_argv}

perl $Bin/JobProcessor.pl $OUTDIR/$sub_argv $Analysis_Ckp{$sub_argv} cnv.$runchr_group $Method_CNV
SCRIPTS
				print FILE_HANDLE_CALL $initial_cnv_call;
				close FILE_HANDLE_CALL;
				undef $initial_cnv_call;
				
				
				#------------------------------------------------------
				# CNVnator call
				#------------------------------------------------------
				open FILE_HANDLE_CALL, ">$Process_Path{$sub_argv}/shells/step_2-3-5_call_cnv_$Method_CNV.$runchr_group\_call.sh" or die $!;
				$initial_cnv_call = <<"SCRIPTS";
#!/bin/sh
## cnv for $sub_argv
#\$ -cwd -l virtual_free=2,h_vmem=20G -pe smp 3
export ROOTSYS=/mnt/isilon/cbttc/tools/root
export PATH=\$PATH:\$ROOTSYS/bin
export SHLIB_PATH=\$SHLIB_PATH:\$ROOTSYS/lib
export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:\$ROOTSYS/lib:/cm/shared/apps/gcc/4.7.0/lib64
export LIBPATH=\$LIBPATH:\$ROOTSYS/lib
{ $ReadCfg{CNVnator} -root $Process_Path{$sub_argv}/runs/step_2-3_cnv/$sub_argv.$runchr_group.root -call $ReadLen{$sub_argv} > $Process_Path{$sub_argv}/runs/step_2-3_cnv/$sub_argv.$runchr_group.raw.cnv; } || { while [ -e $Analysis_Log{$sub_argv}.tmp ]; do sleep 1; done; perl $Bin/UpdateLog.pl $Analysis_Log{$sub_argv} "call cnv call $runchr_group" 3 > $Analysis_Log{$sub_argv}.tmp && mv $Analysis_Log{$sub_argv}.tmp $Analysis_Log{$sub_argv}; exit 1; }
	
while [ -e $Analysis_Log{$sub_argv}.tmp ]; do sleep 1; done; perl $Bin/UpdateLog.pl $Analysis_Log{$sub_argv} "call cnv call $runchr_group" 0 > $Analysis_Log{$sub_argv}.tmp && mv $Analysis_Log{$sub_argv}.tmp $Analysis_Log{$sub_argv}

while [ -e $Analysis_Ckp{$sub_argv}.tmp ]; do sleep 1; done; perl $Bin/UpdateCkp.pl $Analysis_Ckp{$sub_argv} cnv.$runchr_group $Method_CNV 6 > $Analysis_Ckp{$sub_argv}.tmp && mv $Analysis_Ckp{$sub_argv}.tmp $Analysis_Ckp{$sub_argv}

perl $Bin/JobProcessor.pl $OUTDIR/$sub_argv $Analysis_Ckp{$sub_argv} cnv.$runchr_group $Method_CNV
SCRIPTS
				print FILE_HANDLE_CALL $initial_cnv_call;
				close FILE_HANDLE_CALL;
				undef $initial_cnv_call;
			}
		}
	}
}

sub _svScripts {
	my $sub_argv = shift;

	#--------------------------------------------------------------
	# initial sv scripts
	#--------------------------------------------------------------
	switch ( $Method_SV ) {
		case ('CREST') {
			#------------------------------------------------------
			# structural variation calling
			#------------------------------------------------------
			open FILE_HANDLE_CALL, ">$Process_Path{$sub_argv}/shells/step_2-4-1_call_sv_$Method_SV.sh" or die $!;
			my $initial_sv_call = <<"SCRIPTS";
#!/bin/sh
## stru for idv $sub_argv
#\$ -cwd -l virtual_free=30G,h_vmem=50G -pe smp 8

{ extractSClip.pl -i $Process_Path{$sub_argv}/runs/step_1-1_map/$sub_argv.bam --ref_genome $ReadCfg{$Method_REF} && CREST.pl -f $Process_Path{$sub_argv}/runs/step_1-1_map/$sub_argv.bam.cover -d $Process_Path{$sub_argv}/runs/step_1-1_map/$sub_argv.bam --ref_genome $ReadCfg{$Method_REF} -t $ReadCfg{$Method_REF}.2bit; } || { while [ -e $Analysis_Log{$sub_argv}.tmp ]; do sleep 1; done; perl $Bin/UpdateLog.pl $Analysis_Log{$sub_argv} "call sv" 3 > $Analysis_Log{$sub_argv}.tmp && mv $Analysis_Log{$sub_argv}.tmp $Analysis_Log{$sub_argv}; exit 1; }
	
while [ -e $Analysis_Log{$sub_argv}.tmp ]; do sleep 1; done; perl $Bin/UpdateLog.pl $Analysis_Log{$sub_argv} "call sv" 0 > $Analysis_Log{$sub_argv}.tmp && mv $Analysis_Log{$sub_argv}.tmp $Analysis_Log{$sub_argv}

while [ -e $Analysis_Ckp{$sub_argv}.tmp ]; do sleep 1; done; perl $Bin/UpdateCkp.pl $Analysis_Ckp{$sub_argv} sv $Method_SV 3 > $Analysis_Ckp{$sub_argv}.tmp && mv $Analysis_Ckp{$sub_argv}.tmp $Analysis_Ckp{$sub_argv}

perl $Bin/JobProcessor.pl $OUTDIR/$sub_argv $Analysis_Ckp{$sub_argv} sv $Method_SV
SCRIPTS
			print FILE_HANDLE_CALL $initial_sv_call;
			close FILE_HANDLE_CALL;
			undef $initial_sv_call;
		}
		case ('Delly'){
			#------------------------------------------------------
			# delly call DEL
			#------------------------------------------------------
			foreach my $type ('DEL', 'DUP', 'INV', 'TRA') {
				open FILE_HANDLE_CALL, ">$Process_Path{$sub_argv}/shells/step_2-4-1_call_sv_$Method_SV\_$type.sh" or die $!;
				my $initial_sv_call = <<"SCRIPTS";
#!/bin/sh
## stru for idv $sub_argv
#\$ -cwd -l virtual_free=30G,h_vmem=50G -pe smp 8
export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:/mnt/isilon/cbttc/tools/boost-1.55/lib:/mnt/isilon/cbttc/tools/bamtools-master/lib:/cm/shared/apps/gcc/4.7.0/lib64
{ $ReadCfg{Delly} -t $type -o $Process_Path{$sub_argv}/runs/step_2-4_sv/$sub_argv.$type.vcf -g $ReadCfg{$Method_REF} $ReadCfg{"delly.excl.$Method_REF"} $Process_Path{$sub_argv}/runs/step_1-1_map/$sub_argv.bam; echo 'sv $type ok'; } || { echo 'sv $type failed'; while [ -e $Analysis_Log{$sub_argv}.tmp ]; do sleep 1; done; perl $Bin/UpdateLog.pl $Analysis_Log{$sub_argv} "call sv $type" 3 > $Analysis_Log{$sub_argv}.tmp && mv $Analysis_Log{$sub_argv}.tmp $Analysis_Log{$sub_argv}; exit 1; }
	
while [ -e $Analysis_Log{$sub_argv}.tmp ]; do sleep 1; done; perl $Bin/UpdateLog.pl $Analysis_Log{$sub_argv} "call sv $type" 0 > $Analysis_Log{$sub_argv}.tmp && mv $Analysis_Log{$sub_argv}.tmp $Analysis_Log{$sub_argv}

while [ -e $Analysis_Ckp{$sub_argv}.tmp ]; do sleep 1; done; perl $Bin/UpdateCkp.pl $Analysis_Ckp{$sub_argv} sv $Method_SV 2 > $Analysis_Ckp{$sub_argv}.tmp && mv $Analysis_Ckp{$sub_argv}.tmp $Analysis_Ckp{$sub_argv}

SCRIPTS
				print FILE_HANDLE_CALL $initial_sv_call;
				close FILE_HANDLE_CALL;
				undef $initial_sv_call;
			}
		}
	}
}

=Description
	Sequencing Data Analysis Pipeline for CBTTC
	Yuankun Zhu, ZhuY@email.chop.edu

Options:
	-m	input manifest files list file [forced]
	-c	analysis type Config file [forced]
	-o	output dir
	-u  update scripts only [0/1 for debug]
	-help|?               output help information to screen

Example
  perl SeqAnPipeline.pl -m <fq file manifeset> -c <analysis type config> -o <out put dir>

=cut
