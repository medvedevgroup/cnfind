#!/usr/bin/perl
use strict;

#Use this script when needing to run cnfind on multiple samples from the same individual. It will make figures.
#The script does not run anything but prints a series of commands which you can then pipe to sh or run in parallel as need be
#The project file (first argument) contains in each row a bam file and its "name"
#The first line is treated as the normal
#The second argument is the output directory

my $project_file = $ARGV[0];
my $work_base = $ARGV[1];

(-e $project_file) or die("Project file \"$project_file\" not found\n");

my $alt_work_dir;
my $singlestage;


my $starting_mode = "stage1";
my $expectedCalc = "gw";
my $normalization = "self";  
my $win_size = 10000;

$alt_work_dir = "--alt_work_dir $work_base";

#$singlestage = "--singlestage";



my @colors = ("green", "cyan", "orange", "pink");

my @samples;
open(PROJECT, $project_file);
while (<PROJECT>) {
	chomp;
	push(@samples, $_);
}
close(PROJECT);

my @normInfo     = split / /, $samples[0];
my $normBam      = $normInfo[0];
my $normName     = $normInfo[1];

if ($starting_mode eq "stage1") {
	#first run the normal all the way and everything else in only stage1. This is because stage1 takes a long time, but in stage2 the tumor calling will depend on the normal being completed. 
	for (my $i = 1; $i < @samples; $i++) {
		my @info = split / /, $samples[$i];
		my $bam = $info[0];
		my $name     = $info[1];
		print "cnfind.pl $alt_work_dir/$name --work_dir $work_base/$name --mode $starting_mode --bam_files $bam --normalization $normalization --expectedCalc $expectedCalc --par --win_size $win_size --pval 1e-10 --singlestage & \n";
	}
	print "cnfind.pl $alt_work_dir/$normName --work_dir $work_base/$normName --mode $starting_mode --bam_files $normBam --normalization self --expectedCalc $expectedCalc --par --win_size $win_size --minlogratio -0.15 --maxlogratio 0.15\n";
}


if ($starting_mode eq "stage1") {
	$starting_mode = "stage2";
}

#first run normal 
print "cnfind.pl $alt_work_dir/$normName --work_dir $work_base/$normName --mode $starting_mode --bam_files $normBam --normalization self --expectedCalc $expectedCalc --par --win_size $win_size --minlogratio -0.15 --maxlogratio 0.15 $singlestage\n";

#next run the rest of the samples
my $tum_files;
for (my $i = 1; $i < @samples; $i++) {
	my @info = split / /, $samples[$i];
	my $bam = $info[0];
	my $name     = $info[1];
	$tum_files .= "$work_base/$name $name $colors[$i] ";
	print "cnfind.pl $alt_work_dir/$name --work_dir $work_base/$name --mode $starting_mode --bam_files $bam --normalization $normalization --expectedCalc $expectedCalc --par --win_size $win_size --pval 1e-10 $singlestage\n";
}

#TODO plots
#combined plots
print "cnfind.pl --mode plot_doc --fig_name $work_base/all_nocalls.png --samples \"$work_base/$normName $normName $colors[0] $tum_files\"\n";
print "cnfind.pl --mode plot_doc --fig_name $work_base/all_nobuffy.png --samples \"$tum_files\"\n";
print "cnfind.pl --mode plot_doc --fig_name $work_base/tums_nocalls.png --samples \"$tum_files\" \n";

print "cnfind.pl --mode plot_doc --fig_name $work_base/$normName.png --samples \"$work_base/$normName $normName $colors[0]\" --work_dir $work_base/$normName\n";
for (my $i = 1; $i < @samples; $i++) {
	my @info = split / /, $samples[$i];
	my $name     = $info[1];
	print "cnfind.pl --mode plot_doc --fig_name $work_base/${name}_$normName.png --samples \"$work_base/$name $name $colors[$i] $work_base/$normName $normName $colors[0]\" --work_dir $work_base/$name\n";
}

