#!/usr/bin/perl

# Copyright 2012, Paul Medvedev
#
# This file is part of Cnfind
#
# Cnfind is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Cnfind is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with cnfind (see file gpl.txt). If not, see <http://www.gnu.org/licenses/>.

use Getopt::Long;
use strict;
use POSIX;
use File::Basename;
use File::Path;
use File::Copy;
use File::HomeDir;
use File::Temp qw/ tempdir /;
use FindBin;

my @bgPids;
my $mode;
my $work_dir;
my %options=();
my $threads = 8;
my $size_gc_win  = 1000;
my $num_gc_bins  = 40;
my $bam_files;
my $chr;
my $chr_names_filename;
my $def_ref_folder = File::HomeDir->my_home . "/data/hg19comp/";
my $ref_folder;
my @all_chroms;
my $chr_len_file;
my $all_chromosomes_file;
my $Rfolder = File::HomeDir->my_home . "/bin/";
my $exec_folder = $FindBin::Bin; 
my $masks_prefix = "76mer";
my $alt_work_dir;
my $cgh_file;
my $nogc = 0;
my $singlestage= 0;
my $rmap_short_chr = 0;
my $chr_format = "auto";
my $samples2plot;
my $fig_name;
my $temp_dir = tempdir();
my $snp_file;
my $pval = 0.05;
my $gcbins;
my $par = 0;
my $calls2plot;
my $male = 0;
my $normal_dir;
my $normal_sd = 0;
my $window_size = 1000000; 
my $minlogratio = 0;
my $maxlogratio = 1000000;
my @working_chroms;
my $normalization = "gw";
my $chr_len_file;
my $chr_len;          
my $mask_file;         
my $scov_file;       
my $exp_file;        
my $win_file;       
my $fasta_file;    
my $rmap_file;
my $plot_file_base;


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
	my ($command, $bg) = @_;
	print LOGFILE getTime() ."\t\t". $command . "\n";
	print "Exec at " . getTime() . " $command...\n";
	#$retval  = `$command`;
	if ($bg) {
		my $pid = fork();
		die "unable to fork: $!" unless defined($pid);
		if (!$pid) {  # child
			exec($command);
			die "unable to exec: $!";
		}
		push(@bgPids, $pid);
	} else {
		system($command);
		if ( $? == -1 || $? == 1 ) {
			print "command failed: $!";
		}
	}
}

sub reconvene {
	for my $pid (@bgPids) {
		waitpid $pid, 0;
	}
	@bgPids = {};
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

sub setFileNames {
	my ($chr, @tryAlt) = @_;
	mkpath("$work_dir/$chr");
	$chr_len              = `cat $chr_len_file | awk '{ if (\$1 == "$chr") print \$2; }'`;
	$mask_file            = "$ref_folder/gem_mappability/$masks_prefix/$chr.mask";
	$fasta_file           = "$ref_folder/fasta_files_folder/$chr.fa";
	$rmap_file            = "$work_dir/mapping_files/$chr.rmap";
	$scov_file            = "$work_dir/$chr/$chr.scov";
	$exp_file             = "$work_dir/$chr/$chr.exp"; 
	$win_file             = "$work_dir/$chr/$chr.win"; 
	$plot_file_base       = "$work_dir/$chr/$chr.covplot";
	foreach my $extension (@tryAlt) {
		if ($extension eq "scov") {
			if (!(-e $scov_file) && defined($alt_work_dir)) {
				$scov_file = "$alt_work_dir/$chr/$chr.scov";
			}
		}
		if ($extension eq "exp") {
			if (!(-e $exp_file) && defined($alt_work_dir)) {
				$exp_file = "$alt_work_dir/$chr/$chr.exp";
			}
		}
		if ($extension eq "rmap") {
			if (!(-e $rmap_file) && defined($alt_work_dir)) {
				$rmap_file = "$alt_work_dir/mapping_files/$chr.rmap";
			}
		}
		if ($extension eq "win") {
			if (!(-e $win_file) && defined($alt_work_dir)) {
				$win_file = "$alt_work_dir/$chr/$chr.win";
			}
		}
	}
}


sub usage{
	die("Cnfind usage: 
		Mandatory options:
		\t--mode: can be all:stage1:stage2:make_gcbins:get_rmse:renormalize:make_win:segment:plot_doc
		\t--work_dir: output directory
		\t--ref_folder: location of companion package
		General options:
		\t--bam_files: either the name of a bam file or the name of a text files listing all the bam files (one per line)
		\t--singlestage: only run the stage specified by mode, and do not continue further
		\t--chr: limit run to one chromosome
		\t--par: parallelize over chroms whenever possible. 
		\t--masks: prefix of masking files (def: 76mer)
		\t--chr_format: \"auto\", \"short\", or \"long\". Specifices whether ref names in BAM are given as \"chr10\" (long) or \"10\" (short). (valid for stage1 only) (def: auto)
		Normalization and segmentation options:
		\t--normalization: \"auto\" or \"own_chr\" or \"chr*\" or \"gw\"(def: gw)
		\t--num_gc_bins: number of bins to use for gc correction (def: 40)
		\t--win_size: size of windows to be used (def: 1mil)
		Calling options:
		\t--pval: desired significance level for calls (def: 0.05)
		\t--normal_dir: distribution to use for calculation pvals for calls
		\t--minlogratio: when normal_dir is not specified, the minimum logratio to consider normal
		\t--maxlogratio: the maximum ratio to consider normal
		\t--male: expect half the coverage on chrX and coverage on chrY (not yet supported)
		Special options for ploting:
		\t--calls2plot : callfile to plot
		\t--cgh_file: for plotting cgh data (not currently supported)
		\t--snp_file: for plotting snp data
		\t--samples: for plotting multiple samples. String containing \"Sample1Dir Sample1Name ColToPlot1 Sample2Dir...\"
		\t--fig_name: name of output figure for multiple samples
		Non-supported internal debugging options (they do not work as intended so do not use):
		\t--threads (not currently supported)(def: 8)
		\t--rmap_short_chr : chr names in rmap files are writen without \'chr\' prefix (valid in mode = stage2)
		\t--nogc: disable gc correction 
		\t--gcbins: file to use for gcbins instead of generating a new one 
		\t--alt_work_dir: take input files needed for stage from an alternate directory
		Please see README for mode details
		\n");}



my $prog_args = "";
foreach (@ARGV) { $prog_args.= "$_ " }

GetOptions (
	'mode=s' => \$mode, 
	'bam_files=s' => \$bam_files, 
	'work_dir=s' => \$work_dir, 
	'threads=i' => \$threads, 
	'chr=s' => \$chr, 
	'ref_folder=s' => \$ref_folder, 
	'masks=s' => \$masks_prefix, 
	'nogc' => \$nogc, 
	'alt_work_dir=s' => \$alt_work_dir, 
	'cgh_file=s' => \$cgh_file, 
	'snp_file=s' => \$snp_file, 
	'rmap_short_chr' =>\$rmap_short_chr, 
	'chr_format=s', \$chr_format, 
	'samples=s' => \$samples2plot, 
	'pval=f' => \$pval, 
	'minlogratio=f' => \$minlogratio, 
	'maxlogratio=f' => \$maxlogratio, 
	'normal_dir=s' => \$normal_dir, 
	'par' => \$par, 
	'gcbins=s' => \$gcbins, 
	'num_gc_bins=i' => \$num_gc_bins, 
	'calls2plot=s' => \$calls2plot, 
	'fig_name=s' => \$fig_name, 
	'male' => \$male, 
	'win_size=i' => \$window_size, 
	'singlestage' => \$singlestage, 
	'normalization=s' => \$normalization
) or usage();

open ( LOGFILE, ">>log.txt") or die ("$0 : failed to open log file for output: $!\n");
print LOGFILE getTime() . "\t $0 called with $prog_args\n";

if (!defined($ref_folder) && (-d $def_ref_folder)) {
	$ref_folder = $def_ref_folder;
}

if (!defined($mode) || !defined($ref_folder)) {
	usage();
}

if (!defined($work_dir) && ($mode ne "plot_doc")) {
	usage();
} elsif (defined($work_dir)) {
	mkpath($work_dir);
}


(-d $ref_folder)         or die( "Folder $ref_folder does not exist!");
$chr_len_file         = checkExists("$ref_folder/chrom_lengths", "file");
$all_chromosomes_file = checkExists("$ref_folder/allchr.txt", "file");
open(CHR_NAMES, $all_chromosomes_file);
while (<CHR_NAMES>) {
	chomp;
	my $locChr = $_;
	if ($male || ($locChr ne "chrY")) {
		push(@all_chroms, $locChr);
	}
}
close(CHR_NAMES);

if (defined($chr)) {
	@working_chroms = ($chr);
} else {
	@working_chroms = @all_chroms;
}



if ($mode eq "stage1") {
	goto STAGE1;
} elsif ($mode eq "stage2") {
	goto STAGE2;
} elsif ($mode eq "make_gcbins") {
	goto MAKE_GCBINS;
} elsif ($mode eq "get_rmse") {
	goto GET_RMSE;
} elsif ($mode eq "renormalize") {
	goto RENORMALIZE;
} elsif ($mode eq "make_win") {
	goto MAKE_WIN;
} elsif ($mode eq "segment") {
	goto SEGMENT;
} elsif ($mode eq "plot_doc") {
	goto PLOT_DOC;
} elsif ($mode eq "collect") {
	goto COLLECT;
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
if (!defined($bam_files)) {
	die("Error: bam_files not provided.\n");
}
checkExists($bam_files, "file");
if ($bam_files =~ /.*\.bam/) {
	print LOGFILE "Using bam file $bam_files\n";
	`echo $bam_files > $work_dir/map_list.txt`;
	$bam_files = "$work_dir/map_list.txt";
} else {
	print LOGFILE "Using map_list file $bam_files\n";
}

execCommand("$exec_folder/maps2bin $bam_files $all_chromosomes_file bam $chr_format", 0);

exit (0) if ($singlestage);


####################################
###   STAGE2                     ### 
###   Files used:                ###
###   Files created:             ###
####################################
STAGE2:
foreach $chr (@working_chroms) { 
	setFileNames($chr, "rmap");

	checkExists($rmap_file, "file");

	my $chr_lbl = $chr;
	if ($rmap_short_chr) { 
		$chr_lbl = substr $chr, 3; 
	}

	execCommand("$exec_folder/make_simple_coverage $rmap_file $scov_file $chr_lbl $chr_len", $par);
}
reconvene();

exit (0) if ($singlestage);



####################################
###   MAKE_GCBINS                ### 
###   Files used:                ###
###   Files created:             ###
####################################
MAKE_GCBINS:
foreach $chr (@working_chroms) { 
	setFileNames($chr, "scov");
	checkExists($scov_file, "file");
	checkExists($mask_file, "file");
	checkExists($fasta_file, "file");
	if ($nogc) {
		$fasta_file = "$ref_folder/fasta_files_folder/random/$chr.fa";
	}

	execCommand("$exec_folder/make_gcbins $fasta_file $mask_file $scov_file $work_dir/$chr/$chr.gcbins $size_gc_win $num_gc_bins", $par);
}
reconvene();
exit (0) if ($singlestage);

####################################
###   GET_RMSE                   ### 
###   Files used:                ###
###   Files created:             ###
####################################
GET_RMSE:
foreach $chr (@working_chroms) { 
	setFileNames($chr, "scov");
	checkExists($scov_file, "file");
	checkExists($mask_file, "file");
	checkExists($fasta_file, "file");
	if ($nogc) {
		$fasta_file = "$ref_folder/fasta_files_folder/random/$chr.fa";
	}

	my $gcbins_file = "$work_dir/$chr/$chr.gcbins";
	if (defined($gcbins)) {
		checkExists($gcbins, "file");
		$gcbins_file = $gcbins
	}
	checkExists($gcbins_file, "file");

	execCommand("$exec_folder/make_exp $fasta_file $mask_file $scov_file $gcbins_file $size_gc_win $num_gc_bins $exp_file", $par);
}
reconvene();


exit (0) if ($singlestage);

####################################
###   RENORMALIZE                ### 
###   Files used:                ###
###   Files created:             ###
####################################
RENORMALIZE:
my $best_chr;
if ($normalization eq "auto") {
	#find chr with min value in chr.exp.rmse file
	my $best_min = 100000;
	$best_chr = "chrHoward";
	foreach $chr (@all_chroms) {
		next if ($chr eq "chrX" || $chr eq "chrY");
		my $rmse = `cat $work_dir/$chr/$chr.exp.rmse`;
		if ($best_min > $rmse) {
			$best_min = $rmse;
			$best_chr = $chr;
		}
	}
} elsif ($normalization eq "gw") {
	#take the gcbins for each chr and average the labmdas to get all.gcbins
	#just to be safe, only use the autosomes
	my @tothits;
	my @totlen;
	my @rangestart;
	my $header;
	foreach $chr (@all_chroms) {
		next if ($chr eq "chrX" || $chr eq "chrY");
		open (GCBINS, "<$work_dir/$chr/$chr.gcbins") or die ("$0 : failed to open gcbins file for output: $!\n");
		while (<GCBINS>) {
			my @entry = split ' ', $_;
			#my @entry = split /\t/, $_;
			my $bin = $entry[0];
			$rangestart[$bin] = $entry[1];
			my $hits = $entry[3];
			my $len  = $entry[4];
			if ($bin eq "Bin") {
				$header = join("\t", @entry);
				next;
			}
			$tothits[$bin] += $hits;
			$totlen[$bin] += $len;
		}
		close(GCBINS);
	}
	open (GCBINS, ">$work_dir/all.gcbins") or die ("$0 : failed to open $work_dir/all.gcbins file for output: $!\n");
	print GCBINS "$header\n";
	for (my $i = 0; $i < $num_gc_bins; $i++) {
		my $lambda;
		if ($totlen[$i] == 0) {
			$lambda = -2;
		} else {
			$lambda = $tothits[$i] / $totlen[$i];
		}
		print GCBINS "$i\t$rangestart[$i]\t$lambda\t$tothits[$i]\t$totlen[$i]\n";
	}
	close (GCBINS);
}


foreach $chr (@working_chroms) { 
	setFileNames($chr);
	my $gcbins_file;
	if ($normalization eq "own_chr") {
		$gcbins_file = checkExists("$work_dir/$chr/$chr.gcbins", "file");
	} elsif ($normalization eq "auto") {
		$gcbins_file = checkExists("$work_dir/$best_chr/$best_chr.gcbins", "file");
	} elsif ($normalization eq "gw") {
		$gcbins_file = checkExists("$work_dir/all.gcbins", "file");
	} elsif ($normalization =~ /^chr[0-9XY]/) {
		$gcbins_file = checkExists("$work_dir/$normalization/$normalization.gcbins", "file");
	}
	`echo $normalization > $work_dir/$chr/$chr.normalizationr`;
	execCommand("$exec_folder/make_exp $fasta_file $mask_file $scov_file $gcbins_file $size_gc_win $num_gc_bins $exp_file", $par);
}

reconvene();

exit (0) if ($singlestage);

####################################
###   MAKE_WIN                   ### 
###   Files used:                ###
###   Files created:             ###
####################################
MAKE_WIN:
foreach $chr (@working_chroms) { 
	setFileNames($chr, "exp", "scov");
	checkExists($exp_file, "file");
	checkExists($mask_file, "file");
	checkExists($scov_file, "file");

	my $max_total_window_size = 1.5 * $window_size;
	execCommand("$exec_folder/make_masked_windows $mask_file $window_size $max_total_window_size $chr | $exec_folder/doc_walker $scov_file $exp_file full | $exec_folder/gcContent $ref_folder/fasta_files_folder/$chr.fa | awk '{ print \$1, \$2, \$3, log(\$5) / log(2), \$5, \$7, \$9, \$11, \$13, \$14 }' | $exec_folder/add_first_line \"Chr Start End LogRatio Ratio Observed Expected Masked ObservedFull GC\" | tr ' ' '\\t' > $win_file", $par);
}

reconvene();
exit (0) if ($singlestage);


####################################
###   SEGMENT                    ###
###   Files used:                ###
###   Files created:             ###
####################################

SEGMENT:

if (defined($normal_dir)) {
	my $normal_win1  = "$normal_dir/all.windows";
	$normal_sd   = `cat $normal_win1 | $Rfolder/Rscript $exec_folder/get_mean_sd.R Ratio | $exec_folder/item 2`; 
	chomp $normal_sd;
}

foreach $chr (@working_chroms) { 
	setFileNames($chr, "win");
	checkExists($win_file, "file");
	execCommand("cat $win_file |  awk '{ if (\$2 != last + 1) curchr++;  if (\$1 != \"Chr\") \$1 = \"sub\" curchr; last = \$3; print \$0; }' | tr ' ' '\\t' > $win_file.forsegment", 0); 
	execCommand("$Rfolder/Rscript $exec_folder/segment.R $normal_sd $pval  $win_file.forsegment $work_dir/$chr/$chr.segments $chr LogRatio $minlogratio $maxlogratio", $par );
}

reconvene();
exit (0) if ($singlestage);





####################################
###   COLLECT                    ### 
###   Files used:                ###
###   Files created:             ###
####################################
COLLECT:

execCommand("cat $work_dir/chr*/chr*.win | sort -k1,1 -k2n,2 | uniq > $work_dir/all.windows", 0);
execCommand("cat $work_dir/chr*/chr*.segments | $exec_folder/drop_first_col | sort -k1,1 -k2n,2 | uniq > $work_dir/all.segments", 0);
execCommand("cat $work_dir/all.segments | item 1 2 3 9 | grep -v Start | awk '{ if (\$4 != 0) print \$0 }' | sort -k2n,2 | ~/cnfind/merge_intervals.py 0 > $work_dir/all.calls");


exit (0) if ($singlestage);




####################################
###   PLOT_DOC                   ### 
###   Files used:                ###
###   Files created:             ###
####################################
PLOT_DOC:

if (!defined($samples2plot)) {
	defined($work_dir) || die("For mode plot_doc, either samples or work_dir needs to be specified");
	$samples2plot= "$work_dir Sample black";
}

my @files = split ' ', $samples2plot;
my $num_samples = @files / 3;

foreach $chr (@working_chroms) { 
	setFileNames($chr);
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

	my $chr_len = `cat $chr_len_file | awk '{ if (\$1 == "$chr") print \$2; }'`;
	my $pixel_width = int(1000 * $chr_len / 247249719); #TODO unhardcode chr1 len

	my $call_file50 = "nocalls";
	if (defined($work_dir) && defined($calls2plot)) {
		$call_file50 = $calls2plot; 
	}
	my $segment_file = "nosegments";
	if (defined($work_dir)){
		$segment_file = "$work_dir/$chr/$chr.segments";
	}

	execCommand("$Rfolder/Rscript $exec_folder/plot_doc.R $chr $pixel_width $plot_file_base $snp_file $segment_file nocalls nocalls $call_file50 $whatever\n", $par);
}
reconvene();

if (!defined($chr)) {
	if ($num_samples == 1 && !defined($fig_name)) {
		$fig_name = "$work_dir/all.covplot.png";
	}
	if (!defined($fig_name)) {
		usage();
	}
	`rm -f $fig_name`; 

	if ($male) { 
		execCommand("convert +append $work_dir/chr[123]/chr*.covplot.png                                                                  $temp_dir/compound_plots.all1.png; 
			convert +append $work_dir/chr[4567]/chr*.covplot.png                                                                          $temp_dir/compound_plots.all2.png; 
			convert +append $work_dir/chr[89]/chr*.covplot.png  $work_dir/chr1[0123]/chr*.covplot.png                                     $temp_dir/compound_plots.all3.png; 
			convert +append $work_dir/chr[456789]/chr*.covplot.pn $work_dir/chr2[012]/chr*.covplot.png $work_dir/chr[XY]/chr*.covplot.png $temp_dir/compound_plots.all4.png; 
			convert -append $temp_dir/compound_plots.all?.png $fig_name", 0);
	} else {
		execCommand("convert +append $work_dir/chr[123]/chr*.covplot.png $temp_dir/compound_plots.all1.png; 
			convert +append $work_dir/chr[4567]/chr*.covplot.png $temp_dir/compound_plots.all2.png; 
			convert +append $work_dir/chr[89]/chr*.covplot.png  $work_dir/chr1[0123]/chr*.covplot.png $temp_dir/compound_plots.all3.png; 
			convert +append $work_dir/chr[456789]/chr*.covplot.png $work_dir/chr2[012]/chr*.covplot.png $work_dir/chr[X]/chr*.covplot.png  $temp_dir/compound_plots.all4.png; 
			convert -append $temp_dir/compound_plots.all?.png $fig_name", 0);

	}


	#execCommand("convert -append $temp_dir/compound_plots.chr[1234].png $work_dir/compound_plots1.png;
	#	convert -append $temp_dir/compound_plots.chr[789].png $temp_dir/compound_plots.chr10.png $work_dir/compound_plots2.png;
	#	convert +append $temp_dir/compound_plots.chr1[12].png $temp_dir/compound_plots.all31.png;
	#	convert +append $temp_dir/compound_plots.chr1[34].png $temp_dir/compound_plots.all32.png;
	#	convert +append $temp_dir/compound_plots.chr1[567].png $temp_dir/compound_plots.all33.png;
	#	convert +append $temp_dir/compound_plots.chr1[89].png $temp_dir/compound_plots.chr2[012].png  $temp_dir/compound_plots.all34.png;
	#	convert -append $temp_dir/compound_plots.all3?.png $work_dir/compound_plots3.png;");

}
exit (0) if ($singlestage);




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

