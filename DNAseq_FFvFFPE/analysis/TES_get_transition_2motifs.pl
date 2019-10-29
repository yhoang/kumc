#!/usr/bin/perl -w
# Author Yen Hoang, M.Sc.
# 2013

use strict;
use POSIX qw/strftime/;

# Transition analysis of all subsequent and antecedent C's and G's

print strftime "%Y-%m-%d %H:%M:%S", localtime(time);
print "\n";
print STDERR "./TES_get_transition_2motifs.pl C>T G>A subsequent and antecedent\n";

our @sample = ("2474","2561","2640","2685","2938","3050","3356","4079","4191","CN","CT","DN","DT");

our ( $INTERSECT, $WRITE_C_pre, $WRITE_C_post, $WRITE_G_pre, $WRITE_G_post, $TITLE_C, $TITLE_G);
$TITLE_C = "sample\tFF(CA>TA)\tFFPE(CA>TA)\tFF(CT>TT)\tFFPE(CT>TT)\tFF(CG>TG)\tFFPE(CG>TG)\tFF(CC>TC)\tFFPE(CC>TC)\n";
$TITLE_G = "sample\tFF(GA>AA)\tFFPE(GA>AA)\tFF(GT>AT)\tFFPE(GT>AT)\tFF(GG>AG)\tFFPE(GG>AG)\tFF(GC>AC)\tFFPE(GC>AC)\n";

$WRITE_C_pre = ">/project/results/me/TES/pileup_distribution/C_transition_2motifs_ant.txt";
$WRITE_C_post = ">/project/results/me/TES/pileup_distribution/C_transition_2motifs_sub.txt";
$WRITE_G_pre = ">/project/results/me/TES/pileup_distribution/G_transition_2motifs_ant.txt";
$WRITE_G_post = ">/project/results/me/TES/pileup_distribution/G_transition_2motifs_sub.txt";

open WRITE_C_pre, $WRITE_C_pre or die;
open WRITE_C_post, $WRITE_C_post or die;
open WRITE_G_pre, $WRITE_G_pre or die;
open WRITE_G_post, $WRITE_G_post or die;

print WRITE_C_pre $TITLE_C;
print WRITE_C_post $TITLE_C;
print WRITE_G_pre $TITLE_G;
print WRITE_G_post $TITLE_G;

for (my $o = 0; $o <= $#sample; $o++) { 

		$INTERSECT = "</project/results/me/TES/pileup_distribution/$sample[$o]_FF_FFPE_intersect.distribution";
		
		my ( $FF_C, $FF_T, $FFPE_C, $FFPE_T, $FF_G, $FF_A, $FFPE_G, $FFPE_A, $C_match, $G_match, $prenuc, $it ) = qw/0 0 0 0 0 0 0 0 0 0 0 0/;
		my ( @A_postmotif, @T_postmotif, @G_postmotif, @C_postmotif, @A_premotif, @T_premotif, @G_premotif, @C_premotif );
		my @nuc = ( "A", "T", "G", "C" );
		
		for ( my $i = 0 ; $i < 8 ; $i++ ) {	# 0:CA, 1:CT, 2:CG, 3:CC, 4:TA, 5:TT, 6:TG, 7:TC
			$A_premotif[$i] = 0;
			$T_premotif[$i] = 0;
			$G_premotif[$i] = 0;
			$C_premotif[$i] = 0;
			
			$A_postmotif[$i] = 0;
			$T_postmotif[$i] = 0;
			$G_postmotif[$i] = 0;
			$C_postmotif[$i] = 0;
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
				### look for post motifs CA/CT/CG/CC   , where C>T occured
				if ( $C_match == 1 && $ref eq $nuc[$i] ) {
				
					$C_postmotif[$i*2] = $C_postmotif[$i*2] + $FF_C;
					$C_postmotif[$i*2 + 1] = $C_postmotif[$i*2 + 1] + $FFPE_C;
					
					$T_postmotif[$i*2] = $T_postmotif[$i*2] + $FF_T;
					$T_postmotif[$i*2 + 1] = $T_postmotif[$i*2 + 1] + $FFPE_T;
				}
				
				### look for post motifs GA/GT/GG/GC   , where G>A occured
				if ( $G_match == 1 && $ref eq $nuc[$i] ) {
				
					$G_postmotif[$i*2] = $G_postmotif[$i*2] + $FF_G;
					$G_postmotif[$i*2 + 1] = $G_postmotif[$i*2 + 1] + $FFPE_G;
					
					$A_postmotif[$i*2] = $A_postmotif[$i*2] + $FF_A;
					$A_postmotif[$i*2 + 1] = $A_postmotif[$i*2 + 1] + $FFPE_A;
				}
				
			}
			$C_match = 0; $G_match = 0;
			
			### add depth of A,T,G,C
			if ( $ref eq "C" ) { # C>T
				
				$FF_C = $l[6] + $l[10];
				$FF_T = $l[4] + $l[8];
				
				$FFPE_C = $l[15] + $l[19];
				$FFPE_T = $l[13] + $l[17];
				
				$C_match = 1;
				
				### look for pre motifs AC/TC/GC/CC, where C>T occured
				$C_premotif[$prenuc*2] = $C_premotif[$prenuc*2] + $FF_C;
				$C_premotif[$prenuc*2 + 1] = $C_premotif[$prenuc*2 + 1] + $FFPE_C;
				
				$T_premotif[$prenuc*2] = $T_premotif[$prenuc*2] + $FF_T;
				$T_premotif[$prenuc*2 + 1] = $T_premotif[$prenuc*2 + 1] + $FFPE_T;
				
			} elsif ( $ref eq "G" ) { # G>A
			
				$FF_G = $l[5] + $l[9];
				$FF_A = $l[3] + $l[7];
				
				$FFPE_G = $l[14] + $l[18];
				$FFPE_A = $l[12] + $l[16];
				
				$G_match = 1;
				
				### look for pre motifs AG/TG/GG/GG, where G>A occured
				$G_premotif[$prenuc*2] = $G_premotif[$prenuc*2] + $FF_G;
				$G_premotif[$prenuc*2 + 1] = $G_premotif[$prenuc*2 + 1] + $FFPE_G;
				
				$A_premotif[$prenuc*2] = $A_premotif[$prenuc*2] + $FF_A;
				$A_premotif[$prenuc*2 + 1] = $A_premotif[$prenuc*2 + 1] + $FFPE_A;
			}
			
			### save nucleotide for next round pre motifs
			for ( my $i = 0 ; $i < 4 ; $i++ ) {
				if ($ref eq $nuc[$i] ) { $prenuc = $i; }
			}
			
			++$it;
		}
		close(INTERSECTION);
		
		print "...processed $sample[$o].\n";
		#print "$sample[$o]\t$it";
		
		
		### subsequent motifs
		print WRITE_C_post "$sample[$o]";
		for ( my $i = 0 ; $i < 8 ; $i++ ) {	# 0:CA, 1:CT, 2:CG, 3:CC, 4:TA, 5:TT, 6:TG, 7:TC
			my $conversion_rate = $T_postmotif[$i] / $C_postmotif[$i];
			print WRITE_C_post "\t".$conversion_rate;
		}
		print WRITE_C_post "\n";
		
		print WRITE_G_post "$sample[$o]";
		
		for ( my $i = 0 ; $i < 8 ; $i++ ) {	# 0:GA, 1:GT, 2:GG, 3:CG, 4:AA, 5:AT, 6:AG, 7:AC
			my $conversion_rate = $A_postmotif[$i] / $G_postmotif[$i];
			print WRITE_G_post"\t".$conversion_rate;
		}
		print WRITE_G_post "\n";


		### antecedent motifs 
		print WRITE_C_pre "$sample[$o]";
		for ( my $i = 0 ; $i < 8 ; $i++ ) {	# 0:AC, 1:TC, 2:GC, 3:CC, 4:AT, 5:TT, 6:GT, 7:CT
			my $conversion_rate = $T_premotif[$i] / $C_premotif[$i];
			print WRITE_C_pre "\t".$conversion_rate;
		}
		print WRITE_C_pre "\n";
		
		print WRITE_G_pre "$sample[$o]";
		
		for ( my $i = 0 ; $i < 8 ; $i++ ) {	# 0:AG, 1:TG, 2:GG, 3:CG, 4:AA, 5:TA, 6:GA, 7:CA
			my $conversion_rate = $A_premotif[$i] / $G_premotif[$i];
			print WRITE_G_pre"\t".$conversion_rate;
		}
		print WRITE_G_pre "\n";
		
}
close (WRITE_C_pre);
close (WRITE_G_pre);
close (WRITE_C_post);
close (WRITE_G_post);

print "\n\nDone. \n";
print strftime "%Y-%m-%d %H:%M:%S", localtime(time);
print "\n";

