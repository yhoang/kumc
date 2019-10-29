#!/usr/bin/perl -w
# Author Yen Hoang, M.Sc.
# 2013

use strict;
use POSIX qw/strftime/;

# Transition analysis of all intersected positions

print strftime "%Y-%m-%d %H:%M:%S", localtime(time);
print "\n";
print STDERR "./TES_get_transition_all.pl \n";

our @sample = ("2474","2561","2640","2685","2938","3050","3356","4079","4191","CN","CT","DN","DT");
our ( $INTERSECT, $WRITE, $TITLE );
$TITLE = "sample\tFF_AT\tFFPE_AT\tFF_AG\tFFPE_AG\tFF_AC\tFFPE_AC\tFF_TA\tFFPE_TA\tFF_TG\tFFPE_TG\tFF_TC\tFFPE_TC\tFF_GA\tFFPE_GA\tFF_GT\tFFPE_GT\tFF_GC\tFFPE_GC\tFF_CA\tFFPE_CA\tFF_CT\tFFPE_CT\tFF_CG\tFFPE_CG\n";

$WRITE = ">/project/results/me/TES/pileup_distribution/transition_all.txt";

open WRITE, $WRITE or die;
print WRITE $TITLE;

for (my $o = 0; $o <= $#sample; $o++) {

		$INTERSECT = "</project/results/me/TES/pileup_distribution/$sample[$o]_FF_FFPE_intersect.distribution";

		my ( @FF, @FFPE, @sumFF, @sumFFPE, @conversionFF, @conversionFFPE );
		my @nuc = ( "A", "T", "G", "C" );
		
		
		for ( my $i = 0 ; $i < 4 ; $i++ ) {	# 0:A, 1:T, 2:G, 3:C
		
			$FF[$i] = 0;
			$FFPE[$i] = 0;
			$sumFF[$i] = 0;
			$sumFFPE[$i] = 0;
			
			for ( my $j = 0 ; $j < 4 ; $j++ ) {
				$conversionFF[$i*4+$j] = 0;
				$conversionFFPE[$i*4+$j] = 0;
			}
			
		}
		
		open INTERSECTION, $INTERSECT or die "no $INTERSECT";
		while(<INTERSECTION>){
			chomp();
			my @l = split(/\t/,$_);			# array @l w/o tab
			$l[0] eq "chr" and next;
			my $chr = $l[0];
			my $pos = $l[1];
			my $ref = $l[2];
			
			for ( my $i = 0 ; $i < 4 ; $i++ ) {
				$FF[$i] = $l[3+$i] + $l[7+$i];
				$FFPE[$i] = $l[12+$i] + $l[16+$i];
				
				if ( $ref eq $nuc[$i] ) {
					$sumFF[$i] = $sumFF[$i] + $l[3+$i] + $l[7+$i];
					$sumFFPE[$i] = $sumFFPE[$i] + $l[12+$i] + $l[16+$i];
				}
			}
			
			for ( my $i = 0 ; $i < 4 ; $i++ ) {	
				if ( $ref eq $nuc[$i] ) {
					for ( my $j = 0 ; $j < 4 ; $j++ ) {
						$i == $j and next;
						$conversionFF[$i*4+$j] = $conversionFF[$i*4+$j] + $FF[$j];
						$conversionFFPE[$i*4+$j] = $conversionFFPE[$i*4+$j] + $FFPE[$j];
					}
				
				}
			}
			#print "\n";
		}
		close(INTERSECTION);
		
		print "$sample[$o]";
		
		print WRITE "$sample[$o]";
		for ( my $i = 0 ; $i < 4 ; $i++ ) {	# 0:A, 1:T, 2:G, 3:C
			for ( my $j = 0 ; $j < 4 ; $j++ ) {
				$i == $j and next;
				my $rateFF = $conversionFF[$i*4+$j] / $sumFF[$i];
				my $rateFFPE = $conversionFFPE[$i*4+$j] / $sumFFPE[$i];
				print WRITE "\t$rateFF\t$rateFFPE";
				print "\t$rateFF\t$rateFFPE";
			}
		}
		print WRITE "\n";
		print "\n";
		
}
close (WRITE);

print "\n\nDone. \n";
print strftime "%Y-%m-%d %H:%M:%S", localtime(time);
print "\n";

