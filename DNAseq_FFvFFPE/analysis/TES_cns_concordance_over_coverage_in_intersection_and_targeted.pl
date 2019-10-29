#!/usr/bin/perl -w
# Author Yen Hoang, M.Sc.
# 2013

use strict;
use POSIX qw/strftime/;

# This script counts the coverage of each nucleotide on each intersected positions and compares between FF and FFPE
# Separate count of total intersection and targeted only
# Separate count of concordant and discordant positions

print strftime "%Y-%m-%d %H:%M:%S", localtime(time);
print "\n";
print STDERR "./cns_concordance_over_coverage_in_intersection.pl \n";

my $mincov = 250;

my @sample = ("2474","2561","2640","2685","2938","3050","3356","4079","4191","CN","CT","DN","DT");


our ($FF, $FFPE, $write_total, $write_targeted, $targeted);
our ( @inbed_covFF, @inbed_covFFPE, @total_covFF, @total_covFFPE);

# Initiation
for (my $i = 0; $i <= $mincov; $i++) {       
        $inbed_covFF[$i][0] = 0;     # amount of concordance FF in targeted regions
        $inbed_covFF[$i][1] = 0;     # amount of discordance FF in targeted regions
        $inbed_covFFPE[$i][0] = 0;   # amount of concordance FFPE in targeted regions
        $inbed_covFFPE[$i][1] = 0;   # amount of discordance FFPE in targeted regions
        
        $total_covFF[$i][0] = 0;     # amount of concordance FF
        $total_covFF[$i][1] = 0;     # amount of discordance FF
        $total_covFFPE[$i][0] = 0;   # amount of concordance FFPE
        $total_covFFPE[$i][1] = 0;   # amount of discordance FFPE
}

$targeted = "</project/results/me/TES/targeted_positions.list";

my %bed;
my $it_bed = 0;

print STDERR "opens $targeted... ";
open F2, $targeted or die "no $targeted";
while(<F2>){	#array data
	chomp();
	my @b = split(/\t/,$_);		# array @l w/o tab
	if ($b[0] eq "chr") {next;}
	my $chr = $b[0];
	my $pos = $b[1];
	my $ref = $b[2];
	$bed{$chr}{$pos} = $ref;
	++$it_bed;
}
close(F2);
print STDERR "$it_bed targeted positions.\n\n";
print STDERR "sample\tpos_FF\tpos_FFPE\tintersect_total\tcon_total\tdis_total\tintersect_targeted\tcon_targeted\tdis_targeted\n";
	
for (my $o = 0; $o <= $#sample; $o++) {
        $FF = "</project/results/me/TES/$sample[$o]_FF/$sample[$o]_FF.cns.starless.pileup";	
        $FFPE = "</project/results/me/TES/$sample[$o]_FFPE/$sample[$o]_FFPE.cns.starless.pileup";	
        $write_total = ">/project/results/me/TES/cns_concordance/$sample[$o]_FF_FFPE_cns_cov";
        $write_targeted = ">/project/results/me/TES/cns_concordance/$sample[$o]_FF_FFPE_targeted_cns_cov";
        
        my ($it_FF, $it_FFPE, $it_bedsame, $it_same, $itcon, $itdis, $it_bedcon, $it_beddis) = qw/0 0 0 0 0 0 0 0/;
        my (%FF, %inbed_FF, %help_FF, %FFPE, %inbed_FFPE, %help_FFPE); 

        #######################
        #print STDERR "opens $FF...\n";
        open F1, $FF or die "no $FF";
        while(<F1>){
        	chomp();
	        my @l = split(/\t/,$_);	
	        my $chr = $l[0];
	        my $pos = $l[1];
	        my $ref = $l[2];
	        my $var = $l[3];
	        my $cov = $l[7];
	        if (exists $bed{$chr}{$pos}) { 
	                $inbed_FF{$chr}{$pos} = $var; 
	        }
	        $FF{$chr}{$pos} = $var; 
	        $help_FF{$chr}{$pos} = $cov;
	        ++$it_FF;
        }       
        close(F1);

        #######################
        #print STDERR "opens $FFPE...\n";
        open F1, $FFPE or die "no $FFPE";
        while(<F1>){
	        chomp();
	        my @l = split(/\t/,$_);	
	        my $chr = $l[0];
	        my $pos = $l[1];
	        my $ref = $l[2];
	        my $var = $l[3];
	        my $cov = $l[7];
	        if ( exists $FF{$chr}{$pos} ) {
		        if (exists $inbed_FF{$chr}{$pos	}) { 
		                $inbed_FFPE{$chr}{$pos} = $var; 
		                ++$it_bedsame;
		        }
		        $FFPE{$chr}{$pos} = $var; 
		        $help_FFPE{$chr}{$pos} = $cov;
		        ++$it_same;
		}
		++$it_FFPE;
        }
        close(F1);
	
	foreach my $chr ( sort keys %FF) {  
		if (defined $FFPE{$chr}) {
			foreach my $pos (sort keys %{$FF{$chr}}){
				if (defined $FFPE{$chr}{$pos}){
				        if ($FF{$chr}{$pos} eq $FFPE{$chr}{$pos}) {
				                for (my $i = 0; $i <= $mincov; $i += 2 ) {
					                if ( $i > $help_FF{$chr}{$pos} ) {++$total_covFF[$i][0];}
					        }
					        for (my $i = 1; $i <= $mincov; $i += 2 ) {
					        	if ( $i > $help_FFPE{$chr}{$pos} ) {++$total_covFFPE[$i][0];}
					        }
				                ++$itcon;
					} else {
					        for (my $i = 0; $i <= $mincov; $i += 2 ) {
					                if ( $i > $help_FF{$chr}{$pos} ) {++$total_covFF[$i][1];}
					        }
					        for (my $i = 1; $i <= $mincov; $i += 2 ) {
					        	if ( $i > $help_FFPE{$chr}{$pos} ) {++$total_covFFPE[$i][1];}
					        }
					        ++$itdis;
				        }
        	                }
			}
		}
	}
	
	foreach my $chr ( sort keys %inbed_FF){  
		if (defined $inbed_FFPE{$chr}){
			foreach my $pos (sort keys %{$inbed_FF{$chr}}){
				if (defined $inbed_FFPE{$chr}{$pos}){
				        if ($inbed_FF{$chr}{$pos} eq $inbed_FFPE{$chr}{$pos}) {
				                for (my $i = 0; $i <= $mincov; $i += 2 ) {
					                if ( $i > $help_FF{$chr}{$pos} ) {++$inbed_covFF[$i][0];}
					        }
					        for (my $i = 1; $i <= $mincov; $i += 2 ) {
					        	if ( $i > $help_FFPE{$chr}{$pos} ) {++$inbed_covFFPE[$i][0];}
					        }
				                ++$it_bedcon;
					} else {
					        for (my $i = 0; $i <= $mincov; $i += 2 ) {
					                if ( $i > $help_FF{$chr}{$pos} ) {++$inbed_covFF[$i][1];}
					        }
					        for (my $i = 1; $i <= $mincov; $i += 2 ) {
					        	if ( $i > $help_FFPE{$chr}{$pos} ) {++$inbed_covFFPE[$i][1];}
					        }
					        ++$it_beddis;
				        }
        	                }
			}
		}
	}
	
	open WRITE_TOTAL, $write_total or die "Cannot write into $write_total! No direction!\n";
	print WRITE_TOTAL "cov\tconFF\tconFFPE\n";
	for (my $i = 0; $i <= $mincov; $i++) {
        	my $conFF = 0;
        	if ( ($total_covFF[$i][0] + $total_covFF[$i][1]) > 0) {
        	$conFF = $total_covFF[$i][0] / ($total_covFF[$i][0] + $total_covFF[$i][1]); }
        	my $conFFPE = 0;
        	if ( ($total_covFFPE[$i][0] + $total_covFFPE[$i][1]) > 0) { 
        	$conFFPE = $total_covFFPE[$i][0] / ($total_covFFPE[$i][0] + $total_covFFPE[$i][1]);}
		print WRITE_TOTAL "$i\t$conFF\t$conFFPE\n";
	}
	close(WRITE_TOTAL);

	open WRITE_TARGET, $write_targeted or die "Cannot write into $write_targeted! No direction!\n";
	print WRITE_TARGET "cov\tconFF\tconFFPE\n";
	for (my $i = 0; $i <= $mincov; $i++) {
        	my $conFF = 0;
        	if ( ($inbed_covFF[$i][0] + $inbed_covFF[$i][1]) > 0) {
        	$conFF = $inbed_covFF[$i][0] / ($inbed_covFF[$i][0] + $inbed_covFF[$i][1]); }
        	my $conFFPE = 0;
        	if ( ($inbed_covFFPE[$i][0] + $inbed_covFFPE[$i][1]) > 0) { 
        	$conFFPE = $inbed_covFFPE[$i][0] / ($inbed_covFFPE[$i][0] + $inbed_covFFPE[$i][1]);}
		print WRITE_TARGET "$i\t$conFF\t$conFFPE\n";
	}
	close(WRITE_TARGET);
	
	print STDERR "$sample[$o]\t$it_FF\t$it_FFPE\t$it_same\t$itcon\t$itdis\t$it_bedsame\t$it_bedcon\t$it_beddis\n";
}

print STDERR "Done.\n";
print strftime "%Y-%m-%d %H:%M:%S", localtime(time);
print "\n";
