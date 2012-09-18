#!/usr/bin/perl
use Getopt::Long;
use strict;
use POSIX;
use File::Basename;
use File::Path;
use File::Copy;
use File::HomeDir;
use File::Temp qw/ tempdir /;
use FindBin;

sub getTime{
	my @months = qw(Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec);
	my @weekDays = qw(Sun Mon Tue Wed Thu Fri Sat Sun);
	my ($second, $minute, $hour, $dayOfMonth, $month, $yearOffset, $dayOfWeek, $dayOfYear, $daylightSavings) = localtime();
	my $year = 1900 + $yearOffset;
	my $theTime = "$hour:$minute:$second, $weekDays[$dayOfWeek] $months[$month] $dayOfMonth, $year";
	return $theTime; 

}

sub execCommand {
	my $retval;
	my ($command) = @_;
	print LOGFILE getTime() ."\t\t". $command . "\n";
	print "Exec at " . getTime() . " $command...\n";
	#$retval  = `$command`;
	system($command);
	if ( $? == -1 || $? == 1 ) {
		print "command failed: $!";
	}
	chomp $retval;
	print $retval . "\n";
	return $retval;
}


sub fasta_length {
	my ($filename) = @_;
	my $lines  = `cat $filename | grep -v ">" | wc -l`;
	my $chars  = `cat $filename | grep -v ">" | wc -c`;
	my $length = $chars - $lines;
	return $length;
}

sub checkExists {
	my ($filename, $type) = @_;
	if ($type eq "dir") {
		(-d $filename)       or die( "Folder $filename does not exist!");
	} elsif ($type eq "file") {
		(-e $filename)       or die( "File $filename does not exist!");
	} else {
		die ("checkExists: unknown type ($type)\n");
	}
	return $filename;
}


sub usage{
	die("Usage: \nMandatory options:\n\t--mode: can be all:stage1:stage2:make_gcbins:make_exp:make_win:segment:finish:plot_doc_multi\n\t--work_dir: output directory\n\t--ref_folder: location of companion package\nOptional options:\n\t--chr: limit run to one chromosome\n\t--par : parallelize over chroms, only valid if no chr specified (valid for stage2)\n\t--map_list: text files listing all the bam files\n\t--short_chr : chr names in bam files do not contain \'chr\' prefix, e.g. \'10\' instead of \'chr10\' (valid in mode <= stage2)\n\t--masks: prefix of masking files (def: 76mer)\n\t--win_size : size of windows to be used (def: 1mil)\n\t--nogc: disable gc correction (not supported for external use)\n\t--gcbins : file to use for gcbins instead of generating a new one (relevant to stage <= make_exp)\n\t--num_gc_bins : number of bins to use for gc correction (relevant to stage <= make_exp) (def: 40)\n\t--pval : desired significance level for calls (def: 0.05)\n\t--normal_dir : distribution to use for calculation pvals for calls\n\t--minlogratio : when normal_dir is not specified, the minimum logratio to consider normal\n\t--maxlogratio : the maximum ratio to consider normal\n\t--threads (not currently supported)(def: 8)\nSpecial options for ploting:\n\t--calls2plot : callfile to plot\n\t--cgh_file: for plotting cgh data (not currently supported)\n\t--snp_file: for plotting snp data\n\t--samples: for plotting multiple samples. String containing \"Sample1Dir Sample1Name ColToPlot1 Sample2Dir...\"\n");
}

my $mode;
my $work_dir;
my %options=();
my $threads = 8;
my $size_gc_win  = 1000;
my $num_gc_bins  = 40;
my $map_list;
my $chr;
my $chr_names_filename;
my $def_ref_folder = File::HomeDir->my_home . "/data/hg19comp/";
my $ref_folder;
my $Rfolder = File::HomeDir->my_home . "/bin/";
my $exec_folder = $FindBin::Bin; 
my $masks_prefix = "76mer";
my $alt_work_dir;
my $cgh_file;
my $nogc = 0;
my $short_chr = 0;
my $samples2plot;
my $temp_dir = tempdir();
my $snp_file;
my $pval = 0.05;
my $gcbins;
my $par = 0;
my $calls2plot;
my $normal_dir;
my $normal_sd = 0;
my $window_size = 1000000; 
my $min_ratio;
my $max_ratio;





my $message = "";
foreach (@ARGV) { $message .= "$_ " }

GetOptions ('mode=s' => \$mode, 'map_list=s' => \$map_list, 'work_dir=s' => \$work_dir, 'threads=i' => \$threads, 'chr=s' => \$chr, 'ref_folder=s' => \$ref_folder, 'masks=s' => \$masks_prefix, 'nogc' => \$nogc, 'alt_work_dir=s' => \$alt_work_dir, 'cgh_file=s' => \$cgh_file, 'snp_file=s' => \$snp_file, 'short_chr' =>\$short_chr, 'samples=s' => \$samples2plot, 'pval=f' => \$pval, 'minlogratio=f' => \$minlogratio, 'maxlogratio=f' => \$maxlogratio, 'normal_dir=s' => \$normal_dir, 'par' => \$par, 'gcbins=s' => \$gcbins, 'num_gc_bins=i' => \$num_gc_bins, 'calls2plot=s' => \$calls2plot, 'win_size=i' => \$window_size) or usage();

open ( LOGFILE, ">>log.txt") or die ("$0 : failed to open log file for output: $!\n");
print LOGFILE getTime() . "\t $0 called with $message\n";

if (defined($work_dir)) {
	mkpath($work_dir);
}

if (!defined($ref_folder) && (-d $def_ref_folder)) {
	$ref_folder = $def_ref_folder;
}
if (!defined($mode) || !defined($ref_folder)) {
	usage();
}

(-d $ref_folder)       or die( "Folder $ref_folder does not exist!");
my $chr_len_file = checkExists("$ref_folder/chrom_lengths", "file");

my $mask_file            = "$ref_folder/gem_mappability/$masks_prefix/$chr.mask";
my $scov_file            = "$work_dir/$chr/$chr.scov";
my $exp_file             = "$work_dir/$chr/$chr.exp"; 
my $win_file             = "$work_dir/$chr/$chr.win"; 
my $all_chromosomes_file = checkExists("$ref_folder/allchr.txt", "file");
my $fasta_file           = "$ref_folder/fasta_files_folder/$chr.fa";


if ($mode eq "misc") { 
	goto MISC;
} elsif ($mode eq "plot_doc_multi") {
	goto PLOT_DOC_MULTI;
}



if (!defined($work_dir) && ($mode ne "plot_doc_multi")) {
	usage();
}

if (!defined($chr) && $mode ne "stage1") {
	if ($par) {
		open (CMDFILE, ">rjobs.txt") or die ("$0 : failed to open rjobs.txt for output: $!\n");
	}
	open(CHR_NAMES, $all_chromosomes_file) or die("File $all_chromosomes_file does not exist");
	while (<CHR_NAMES>) {
		chomp;
		$chr = $_;
		if ($par) {
			print CMDFILE "$0 $message --chr $chr &\n";
		} else {
			execCommand("$0  $message --chr $chr  ");
		}

	}

	if ($par) {
		print CMDFILE "wait\n";
		close CMDFILE;
		execCommand("cat rjobs.txt | sh");
	}

	if ($mode eq "stage2" || $mode eq "make_exp" || $mode eq "make_win" || $mode eq "segment" || $mode eq "finish") {
		execCommand("cat $work_dir/chr*/chr*.win | sort -k1,1 -k2n,2 | uniq | grep -v chrX | grep -v chrY > $work_dir/all.windows");
		execCommand("cat $work_dir/chr*/chr*.segments | $exec_folder/drop_first_col | sort -k1,1 -k2n,2 | uniq | grep -v chrX | grep -v chrY > $work_dir/all.segments");

		if (!defined($samples2plot)) {
			$samples2plot= "$work_dir Sample LogRatio";
		}
		goto PLOT_DOC_MULTI;
	}

	exit;
}

my $chr_len = `cat $chr_len_file | awk '{ if (\$1 == "$chr") print \$2; }'`;


if ($mode eq "stage1") {
	goto STAGE1;
} elsif ($mode eq "stage2") {
	goto STAGE2;
} elsif ($mode eq "make_gcbins") {
	goto MAKE_GCBINS;
} elsif ($mode eq "make_exp") {
	goto MAKE_EXP;
} elsif ($mode eq "make_win") {
	goto MAKE_WIN;
} elsif ($mode eq "segment") {
	goto SEGMENT;
} elsif ($mode eq "plot_doc_multi") {
	goto PLOT_DOC_MULTI;
} elsif ($mode eq "finish") {
	goto FINISH;
} else {
	print STDERR "Unknown mode: $mode.\n";
	die();
}


####################################
###   STAGE1                     ### 
###   Files used:                ###
###   Files created:             ###
####################################
STAGE1:
mkpath("$work_dir/mapping_files");
chdir("$work_dir/mapping_files");
if (!defined($map_list)) {
	die("Error: map_list not provided.\n");
}
checkExists($map_list, "file");
my $param = "bla";
if ($short_chr) {
	$param = "short_chr";
}
execCommand("$exec_folder/maps2bin $map_list $all_chromosomes_file bam $param");

if ($mode eq "stage1") {
	exit(0);
}


####################################
###   STAGE2                     ### 
###   Files used:                ###
###   Files created:             ###
####################################
STAGE2:
if (!defined($chr)) {
	die("Error: chromosome must be specified for stage2.\n");
}

mkpath("$work_dir/$chr");

my $rmap_file = "$work_dir/mapping_files/$chr.rmap";
if (!(-e $rmap_file)) {
	$rmap_file = "$alt_work_dir/mapping_files/$chr.rmap";
}
checkExists($rmap_file, "file");

my $chr_lbl = $chr;
if ($short_chr) { 
	$chr_lbl = substr $chr, 3; 
}

execCommand("$exec_folder/make_simple_coverage $rmap_file $scov_file $chr_lbl $chr_len");



####################################
###   MAKE_GCBINS                ### 
###   Files used:                ###
###   Files created:             ###
####################################
MAKE_GCBINS:
mkpath("$work_dir/$chr");
checkExists($mask_file, "file");
if ($nogc) {
	$fasta_file = "$ref_folder/fasta_files_folder/random/$chr.fa";
}

if (!(-e $scov_file)) {
	$scov_file = "$alt_work_dir/$chr/$chr.scov";
}
checkExists($scov_file, "file");
checkExists($fasta_file, "file");


execCommand("$exec_folder/make_gcbins $fasta_file $mask_file $scov_file $work_dir/$chr/$chr.gcbins $size_gc_win $num_gc_bins");



####################################
###   MAKE_EXP                   ### 
###   Files used:                ###
###   Files created:             ###
####################################
MAKE_EXP:
mkpath("$work_dir/$chr");
checkExists($mask_file, "file");
if ($nogc) {
	$fasta_file = "$ref_folder/fasta_files_folder/random/$chr.fa";
}
checkExists($fasta_file, "file");

if (!(-e $scov_file)) {
	$scov_file = "$alt_work_dir/$chr/$chr.scov";
}
checkExists($scov_file, "file");
 
my $gcbins_file = "$work_dir/$chr/$chr.gcbins";
if (defined($gcbins)) {
	checkExists($gcbins, "file");
	$gcbins_file = $gcbins
}
checkExists($gcbins_file, "file");

execCommand("$exec_folder/make_exp $fasta_file $mask_file $scov_file $gcbins_file $size_gc_win $num_gc_bins $exp_file");


####################################
###   MAKE_WIN                   ### 
###   Files used:                ###
###   Files created:             ###
####################################
MAKE_WIN:
checkExists($exp_file, "file");
checkExists($mask_file, "file");
checkExists($scov_file, "file");

my $max_total_window_size = 1.5 * $window_size;
execCommand("$exec_folder/make_masked_windows $mask_file $window_size $max_total_window_size $chr | $exec_folder/doc_walker $scov_file $exp_file full | $exec_folder/gcContent $ref_folder/fasta_files_folder/$chr.fa | awk '{ print \$1, \$2, \$3, log(\$5) / log(2), \$5, \$7, \$9, \$11, \$13, \$14 }' | $exec_folder/add_first_line \"Chr Start End LogRatio Ratio Observed Expected Masked ObservedFull GC\" | tr ' ' '\\t' > $win_file");


####################################
###   SEGMENT                    ###
###   Files used:                ###
###   Files created:             ###
####################################

SEGMENT:

execCommand("rm -f $win_file.*.pval");


if (defined($normal_dir)) {
	#my $normal_win1  = "$normal_dir/$chr/$chr.win.1.pval";
	my $normal_win1  = "$normal_dir/all.windows";
	$normal_sd   = `cat $normal_win1 | $Rfolder/Rscript $exec_folder/get_mean_sd.R Ratio | $exec_folder/item 2`; 
	chomp $normal_sd;
	print STDERR "Found sd: $normal_sd.\n";
}


#execCommand("cat $win_file    | $Rfolder/Rscript $exec_folder/add_pvals.R $normal_sd   > $win_file.1.pval");

execCommand("cat $win_file |  awk '{ if (\$2 != last + 1) curchr++;  if (\$1 != \"Chr\") \$1 = \"sub\" curchr; last = \$3; print \$0; }' | tr ' ' '\\t' > $win_file.forsegment"); 
execCommand("$Rfolder/Rscript $exec_folder/segment.R $normal_sd $pval  $win_file.forsegment $work_dir/$chr/$chr.segments $chr LogRatio $minlogratio $maxlogratio" );


#execCommand("cat $work_dir/$chr/$chr.calls.* | grep -v Call | sort -k2n,2 -k3n,3 | merge_intervals.py 0 | tr ' ' '\\t' > $work_dir/$chr/$chr.calls");

FINISH:
my $dummy;

exit();


####################################
###   PLOT_DOC_MULTI             ### 
###   Files used:                ###
###   Files created:             ###
####################################
PLOT_DOC_MULTI:

my @files = split / /, $samples2plot;

open (CMDFILE, ">rjobs.txt") or die ("$0 : failed to open rjobs.txt for output: $!\n");


open(CHR_NAMES, $all_chromosomes_file) or die("File $all_chromosomes_file does not exist");
while (<CHR_NAMES>) {
	chomp;
	$chr = $_;
	my @changed_files;
	for (my $i = 0; $i < scalar(@files); $i++) {
		if (($i % 3) == 0) {
			$changed_files[$i] = "$files[$i]/$chr/$chr.win";
		} else {
			$changed_files[$i] = $files[$i];
		}
	}
	my $whatever = join " ", @changed_files;

	if (defined($snp_file) && ($snp_file ne "nosnps")) {
		checkExists($snp_file, "file");
	} else {
		$snp_file = "nosnps";
	}

	#DEPRECATED: aCGH incorporation into plot
	#my $use_cgh = "";
	#if (defined($cgh_file)) {
	#	checkExists($cgh_file, "file");
	#	execCommand("cat $win_file | $exec_folder/lower_resolution $cgh_file > $win_file.cgh");
	#	$win2_file = "$win_file.cgh";
	#	$use_cgh = "use_cgh";
	#}



	my $chr_len = `cat $chr_len_file | awk '{ if (\$1 == "$chr") print \$2; }'`;
	my $pixel_width = int(1000 * $chr_len / 247249719); 

	my $call_file50 = "nocalls";
	if (defined($work_dir) && defined($calls2plot)) {
		$call_file50 = $calls2plot; 
	}

	print CMDFILE "$Rfolder/Rscript $exec_folder/plot_doc_multi.R $chr $pixel_width $temp_dir/compound_plots.$chr $snp_file $work_dir/$chr/$chr.segments nocalls nocalls $call_file50 $whatever &\n";
}
print CMDFILE "wait\n";
close CMDFILE;
execCommand("rm -f $temp_dir/compound_plot*");
execCommand("rm -f $work_dir/coverage_plot.png");

execCommand("cat rjobs.txt | sh");

execCommand("convert +append $temp_dir/compound_plots.chr[123].png $temp_dir/compound_plots.all1.png; 
	convert +append $temp_dir/compound_plots.chr[4567].png $temp_dir/compound_plots.all2.png; 
	convert +append $temp_dir/compound_plots.chr[89].png $temp_dir/compound_plots.chr1[0123].png $temp_dir/compound_plots.all3.png; 
	convert +append $temp_dir/compound_plots.chr1[456789].png $temp_dir/compound_plots.chr2[012].png $temp_dir/compound_plots.all4.png; 
	convert -append $temp_dir/compound_plots.all?.png $work_dir/coverage_plot.png");

execCommand("convert -append $temp_dir/compound_plots.chr[1234].png $work_dir/compound_plots1.png;
	convert -append $temp_dir/compound_plots.chr[789].png $temp_dir/compound_plots.chr10.png $work_dir/compound_plots2.png;
	convert +append $temp_dir/compound_plots.chr1[12].png $temp_dir/compound_plots.all31.png;
	convert +append $temp_dir/compound_plots.chr1[34].png $temp_dir/compound_plots.all32.png;
	convert +append $temp_dir/compound_plots.chr1[567].png $temp_dir/compound_plots.all33.png;
	convert +append $temp_dir/compound_plots.chr1[89].png $temp_dir/compound_plots.chr2[012].png  $temp_dir/compound_plots.all34.png;
	convert -append $temp_dir/compound_plots.all3?.png $work_dir/compound_plots3.png;");


exit;



MISC:

exit;


=comment
sbsJoin chr1.win /home/pmedvedev/blood/data/TCGA-25-1319-10Arun/chr1/chr1.win | item 1 2 3 4 13 9 | drop_first_line | add_first_line "Chr Start End NoCorrection GCCorrection GCContent" | tr ' ' '\t'> ~/blood/gc_correction_evidence.txt

awk '{ sum += $1; sumsq += $1*$1 } END { printf "Mean: %f, Std: %f\n", sum/NR, sqrt(sumsq/NR - (sum/NR)^2) }' filename

#processing array cgh data
less MSK_252152912680_S01_CGH-v4_10_27Aug08.txt | less -x80 | awk -vFS='\t' '{ print $16, $18, $22 }' > probes
less MSK_252152912680_S01_CGH-v4_10_27Aug08__GCN_V3.mat | sort -k1,1 > ratios
less probes  | sort -k1,1 > probes.s
join ratios probes.s | tr ':' ' ' | item 3 4 2 | awk -vFS="-" '{ print $1, $0 }' | item 1 2 2 5 | sort -k1,1 -k2n,2 > pos2ratio.txt


#processing snvs in KI
cat  ~/reads/KIBu/KIBu_SNV.vcf  | grep -v "\#" | awk '{ print "chr" $1, $2, $4, $5, $10 }'  | tr ':,' ' ' | awk '{ if ($5 == "0/1") print $0 }' |  tr ' ' '\t' > KIBu_SNV.vcf.flt1
cat KI1_SNV.vcf |  awk -vFS="\t" -vOFS="\t" '{ print $1, $2, $3, $4, $21, $22, $25, $26, $27, $29, $30 }' | awk '{ if (NR == 2) print $0; else if ($7 ==0 && $8 == 0 && $9 != "AA" && $9 != "CC" && $9 != "GG" && $9 != "TT" && ($10 + $11 >= 6) ) print "chr" $0 }' | awk '{ if (NR == 1) print $0, "bfreq"; if ($5 + $6 < 6) next; bfreq=($6 / ($5 + $6)); print  $0, bfreq  }' |  tr ' ' '\t' >snps
cat KI1_SNV.vcf |  grep -v NOVEL | awk -vFS="\t" -vOFS="\t" '{ print $1, $2, $3, $4, $21, $22, $25, $26, $27, $29, $30 }' | awk '{ if (NR == 2) print $0; else if ($7 ==0 && $8 == 0 && $9 != "AA" && $9 != "CC" && $9 != "GG" && $9 != "TT" && ($10 + $11 >= 6) ) print "chr" $0 }' | awk '{ if (NR == 1) { print $0, "bfreq"; next }  if (($5 + $6 < 10) || ($10 < 5) || ($11 < 5)) next; bfreq=($6 / ($5 + $6)); print  $0, bfreq  }' |  tr ' ' '\t' > snps3

bu <- read.table("Bu.win.txt", header=T); bu1 <- subset(bu, (bu$End - bu$Start - bu$Masked) == 999999); baseexpect <- sum(bu1$Observed) / (length(bu1$Observed));
 plot(bu1$Ratio, ylim=c(0,2)); points(bu1$Observed / baseexpect, col = "red"); points(bu1$Masked / 1000000, col = "blue"); points(bu1$GC, col="green")


 #collect gc
 cat ~/data/hg19comp/chrom_lengths | item 1 | awk '{ print "cat "$1"/"$1".gcbins | item 3 | drop_first_line | add_first_line " $1 " > gcsummary." $1}' | grep -v chrY | sh; paste gcsummary.* > gcsummary

 #collect calls
 cat chr*/chr*.segments | less -x25 | sort -k2,2 -k3n,3  | uniq | awk '{ if ($10 != 0) print $2, $3, $4,  $10 }' | tr ' ' '\t' > calls

=cut

