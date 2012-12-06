#!/usr/bin/perl
use strict;

my $starting_mode = "make_win";
my $normalization = "normal";
my $work_base = "~/data/kidney/norm2norm";
my $alt_work_dir = "~/data/kidney/";

#my $starting_mode = "collect";
#my $normalization = "chr10";
#my $work_base = "~/data/kidney";



my $project_file = $ARGV[0];


my @colors = ("green", "cyan", "orange", "pink", "gray");

my @samples;
open(PROJECT, $project_file);
while (<PROJECT>) {
	chomp;
	push(@samples, $_);
}
close(PROJECT);


# first is buffy, then cfDNA
my @bufInfo = split / /, $samples[0];
my $bufBam = $bufInfo[0];
my $bufName     = $bufInfo[1];
my @cfInfo = split / /, $samples[1];
my $cfBam = $cfInfo[0];
my $cfName     = $cfInfo[1];
my @info = split / /, $samples[2];
my $metBam = $info[0];
my $metName     = $info[1];
@info = split / /, $samples[3];
my $priBam = $info[0];
my $priName     = $info[1];

print "cnfind.pl --alt_work_dir $alt_work_dir/$bufName --work_dir $work_base/$bufName --mode $starting_mode --bam_files $bufBam --normalization own_chr --minlogratio -0.15 --maxlogratio 0.15 --par \n";
print "cnfind.pl --alt_work_dir $alt_work_dir/$cfName  --work_dir $work_base/$cfName  --mode $starting_mode --bam_files $cfBam  --normalization $normalization --pval 1e-10 --normal_dir $work_base/$bufName/ --par \n";
for (my $i = 2; $i < @samples; $i++) {
	my @info = split / /, $samples[$i];
	my $bam = $info[0];
	my $name     = $info[1];
	print "cnfind.pl --alt_work_dir $alt_work_dir/$name --work_dir $work_base/$name --mode $starting_mode --bam_files $bam --normalization $normalization --minlogratio -0.15 --maxlogratio 0.15 --par \n";
}





my $tum_files;
for (my $i = 2; $i < @samples; $i++) {
	my @info = split " ", $samples[$i];
	my $name     = $info[1];
	$tum_files .= "$work_base/$name $name $colors[$i] ";
}


#combined plots
print "cnfind.pl --mode plot_doc --fig_name $work_base/all_nocalls.png --samples \"$work_base/$bufName $bufName $colors[0] $work_base/$cfName $cfName $colors[1] $tum_files\"\n";
print "cnfind.pl --mode plot_doc --fig_name $work_base/all_nobuffy.png --samples \"$work_base/$cfName $cfName $colors[1] $tum_files\" --work_dir $work_base/$cfName\n";

print "cnfind.pl --mode plot_doc --fig_name $work_base/ser_met_pri.png --samples \"$work_base/$cfName $cfName $colors[1] $work_base/$metName $metName $colors[2] $work_base/$priName $priName $colors[3]\" --work_dir $work_base/$cfName\n";

print "cnfind.pl --mode plot_doc --fig_name $work_base/all_nocalls_nobuffy.png --samples \"$work_base/$cfName $cfName $colors[1] $tum_files\" \n";
print "cnfind.pl --mode plot_doc --fig_name $work_base/tums_nocalls.png --samples \"$tum_files\" \n";

print "cnfind.pl --mode plot_doc --fig_name $work_base/$bufName.png --samples \"$work_base/$bufName $bufName $colors[0]\" --work_dir $work_base/$bufName\n";
for (my $i = 1; $i < @samples; $i++) {
	my @info = split / /, $samples[$i];
	my $name     = $info[1];
	print "cnfind.pl --mode plot_doc --fig_name $work_base/${name}_buffy.png --samples \"$work_base/$name $name $colors[$i] $work_base/$bufName $bufName $colors[0]\" --work_dir $work_base/$name\n";
}

print "Rscript ~/cnfind/figures.R $work_base/ $work_base/\n";

