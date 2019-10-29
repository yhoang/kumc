#!/usr/bin/perl -w
# Author Yen Hoang, M.Sc.
# 2013

use strict;
use POSIX qw/strftime/;

# This script counts the coverage of each nucleotide on each targeted intersected positions and compares between FF and FFPE.
# Separate count of concordant and discordant positions

print strftime "%Y-%m-%d %H:%M:%S", localtime(time);
print "\n";
print STDERR "./cns_concordance_over_coverage_in_intersection.pl \n";

my $mincov = 250;

our @sample = ("2474","2561","2640","2685","2938","3050","3356","4079","4191","CN","CT","DN","DT");

our ($FF, $FFPE, $write_total, $write_targeted, $targeted);
our ( @inbed_covFF, @inbed_covFFPE);

for (my $i = 0; $i <= $mincov; $i++) {       
        $inbed_covFF[$i][0] = 0;     # amount of concordance FF in targeted regions
        $inbed_covFF[$i][1] = 0;     # amount of discordance FF in targeted regions
        $inbed_covFFPE[$i][0] = 0;   # amount of concordance FFPE in targeted regions
        $inbed_covFFPE[$i][1] = 0;   # amount of discordance FFPE in targeted regions
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
print STDERR "$it_bed targeted positions.\n";

print STDERR "intersect_targeted\tcon_targeted\tdis_targeted\t% (con_targeted)\n";

for (my $o = 0; $o <= $#sample; $o++) {
        $FF = "</project/results/me/WES/$sample[$o]_FFex/$sample[$o]_FFex.cns.starless.pileup";	
        $FFPE = "</project/results/me/WES/$sample[$o]_FFPEex/$sample[$o]_FFPEex.cns.starless.pileup";	
        $write_targeted = ">/project/results/me/WES/cns_concordance/$sample[$o]_FFex_FFPEex_targeted_cns_cov";
        
        my ($it_FF, $it_FFPE, $it_bedsame, $it_bedcon, $it_beddis, $div_bedcon) = qw/0 0 0 0 0 0 0 0 0/;
        my ( %inbed_FF, %help_FF, %inbed_FFPE, %help_FFPE); 

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
		my $mapQ = $l[6];
		my $cov = $l[7];
			
		$mapQ == 0 and next;
		
	        if (exists $bed{$chr}{$pos}) { 
	                $inbed_FF{$chr}{$pos} = $var; 
	        	$help_FF{$chr}{$pos} = $cov;
	        }
	        ++$it_FF;
        }       
        close(F1);

        #print STDERR "opens $FFPE...\n";
        open F1, $FFPE or die "no $FFPE";
        while(<F1>){
	        chomp();
	        my @l = split(/\t/,$_);	
	        my $chr = $l[0];
	        my $pos = $l[1];
	        my $ref = $l[2];
	        my $var = $l[3];
		my $mapQ = $l[6];
		my $cov = $l[7];
			
		$mapQ == 0 and next;
		
	        if ( exists $inbed_FF{$chr}{$pos} ) {
		        $inbed_FFPE{$chr}{$pos} = $var; 
		        $help_FFPE{$chr}{$pos} = $cov;
		        ++$it_bedsame;
		}
		++$it_FFPE;
        }
        close(F1);

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
	
	$div_bedcon = $it_bedcon / $it_bedsame *100;
	
	print STDERR "$it_bedsame\t$it_bedcon\t$it_beddis\t$div_bedcon\n";

}

print STDERR "Done.\n";
print strftime "%Y-%m-%d %H:%M:%S", localtime(time);
print "\n";
