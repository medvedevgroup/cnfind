#!/usr/bin/perl
use strict;

my $starting_mode = "make_win";
my $normalization = "gw";
my $work_base = "~/data/thayer/win1k/";
my $win_size = 1000;
my $alt_base = "~/data/thayer/";


my $project_file = $ARGV[0];
(-e $project_file) or die("Project file \"$project_file\" not found\n");


my @colors = ("green", "cyan", "orange", "pink");

my @samples;
open(PROJECT, $project_file);
while (<PROJECT>) {
	chomp;
	push(@samples, $_);
}
close(PROJECT);


foreach my $sample (@samples) {
	my @info = split / /, $sample;
	my $bam = $info[0];
	my $name     = $info[1];
	print "cnfind.pl --alt_work_dir $alt_base/$name --work_dir $work_base/$name --mode $starting_mode --bam_files $bam --normalization $normalization --par --win_size $win_size --minlogratio -0.15 --maxlogratio 0.15 \n";
}

# first is buffy (or, in general, a normal)
my @bufInfo = split / /, $samples[0];
my $bufBam = $bufInfo[0];
my $bufName     = $bufInfo[1];

my $tum_files;
for (my $i = 1; $i < @samples; $i++) {
	my @info = split / /, $samples[$i];
	my $name     = $info[1];
	$tum_files .= "$work_base/$name $name $colors[$i] ";
}


#combined plots
print "cnfind.pl --mode plot_doc --fig_name $work_base/all_nocalls.png --samples \"$work_base/$bufName $bufName $colors[0] $tum_files\"\n";
print "cnfind.pl --mode plot_doc --fig_name $work_base/all_nobuffy.png --samples \"$tum_files\"\n";
print "cnfind.pl --mode plot_doc --fig_name $work_base/tums_nocalls.png --samples \"$tum_files\" \n";

print "cnfind.pl --mode plot_doc --fig_name $work_base/$bufName.png --samples \"$work_base/$bufName $bufName $colors[0]\" --work_dir $work_base/$bufName\n";
for (my $i = 1; $i < @samples; $i++) {
	my @info = split / /, $samples[$i];
	my $name     = $info[1];
	print "cnfind.pl --mode plot_doc --fig_name $work_base/${name}_$bufName.png --samples \"$work_base/$name $name $colors[$i] $work_base/$bufName $bufName $colors[0]\" --work_dir $work_base/$name\n";
}


