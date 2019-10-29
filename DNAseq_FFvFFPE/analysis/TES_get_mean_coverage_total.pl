#!/usr/bin/perl -w
# Author Yen Hoang, M.Sc.
# 2013

use strict;
use POSIX qw/strftime/;

# Count mean coverage

print strftime "%Y-%m-%d %H:%M:%S", localtime(time);
print "\n";
print STDERR "./TES_get_mean_coverage_total.pl  \n";

our @sample = ("2474","2561","2640","2685","2938","3050","3356","4079","4191","CN","CT","DN","DT");

print "sample\tmean_covFF\tmean_covFFPE\tintersection\n";
for (my $o = 0; $o <= $#sample; $o++) {

	my $distribution_file = "</project/results/me/TES/pileup_distribution/$sample[$o]_FF_FFPE_intersect.distribution";
	my ($it, $total_covFF, $total_covFFPE, $mean_covFF, $mean_covFFPE) = qw/0 0 0/;

#	print "opens $distribution_file...\n";
	open DISTRIBUTION, $distribution_file or die "no $distribution_file";
	while(<DISTRIBUTION>){
		chomp();
		my @l = split(/\t/,$_);			# array @l w/o tab
		$l[0] eq "chr" and next; 
		my $chr = $l[0];
		my $pos = $l[1];
	
		my $covFF = $l[3] + $l[4] + $l[5] + $l[6] + $l[7] + $l[8] + $l[9] + $l[10];
		my $covFFPE = $l[12] + $l[13] + $l[14] + $l[15] + $l[16] + $l[17] + $l[18] + $l[19];
	
		$total_covFF += $covFF;
		$total_covFFPE += $covFFPE;
		++$it;
	}
	close(DISTRIBUTION);
	
	$mean_covFF = $total_covFF/$it;
	$mean_covFFPE = $total_covFFPE/$it;
	print "$sample[$o]\t$mean_covFF\t$mean_covFFPE\t$it\n";
}
	
print "\nDone.\n";
print strftime "%Y-%m-%d %H:%M:%S", localtime(time);
print "\n";
