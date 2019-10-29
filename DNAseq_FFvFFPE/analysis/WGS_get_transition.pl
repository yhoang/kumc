#!/usr/bin/perl -w
# Author Yen Hoang, M.Sc.
# 2013

use strict;
use POSIX qw/strftime/;

# Transition analysis of all intersected positions

print strftime "%Y-%m-%d %H:%M:%S", localtime(time);
print "\n";
print STDERR "./WGS_get_transition.pl \n";

our @sample = ("CN");

our ( $INTERSECT, $WRITE_TRANSITION, $TITLE);
$TITLE = "sample\tFF(C)\tFF(T)\tFF(C2T_ratio)\tFFPE(C)\tFFPE(T)\tFFPE(C2T_ratio)\tFF(G)\tFF(A)\tFF(G2A_ratio)\tFFPE(G)\tFFPE(A)\tFFPE(G2A_ratio)\n";
print "sample\tFF(C)\tFF(T)\tFF(C2T_ratio)\tFFPE(C)\tFFPE(T)\tFFPE(C2T_ratio)\tFF(G)\tFF(A)\tFF(G2A_ratio)\tFFPE(G)\tFFPE(A)\tFFPE(G2A_ratio)\n";

$WRITE_TRANSITION = ">/project/results/me/WGS/pileup_distribution/transition.txt";
open WRITE_TV, $WRITE_TRANSITION or die;
print WRITE_TV $TITLE;

for (my $o = 0; $o <= $#sample; $o++) {
		$INTERSECT = "</project/results/me/WGS/pileup_distribution/$sample[$o]_FF_FFPE_intersect.distribution2";
		
		my ( $FF_C, $FF_T, $FFPE_C, $FFPE_T, $FF_G, $FF_A, $FFPE_G, $FFPE_A, $Tv_C_FF, $Tv_C_FFPE, $Tv_G_FF, $Tv_G_FFPE ) = qw/0 0 0 0 0 0 0 0 0 0 0 0/;
		 
		open INTERSECTION, $INTERSECT or die "no $INTERSECT";
		while(<INTERSECTION>){
			chomp();
			my @l = split(/\t/,$_);			# array @l w/o tab
			$l[0] eq "chr" and next;
			my $chr = $l[0];
			my $pos = $l[1];
			my $ref = $l[2];
			
		#	my $covFF = ($l[3] + $l[4] + $l[5] + $l[6]) + ($l[7] + $l[8] + $l[9] + $l[10]);
		#	my $covFFPE = ($l[12] + $l[13] + $l[14] + $l[15]) + ($l[16] + $l[17] + $l[18] + $l[19]);
			
			if ( $ref eq "C" ) { # C>T
			
				$FF_C = $FF_C + $l[6] + $l[10];
				$FF_T = $FF_T + $l[4] + $l[8];
				
				$FFPE_C = $FFPE_C + $l[15] + $l[19];
				$FFPE_T = $FFPE_T + $l[13] + $l[17];
				
			} elsif ( $ref eq "G" ) { # G>A
			
				$FF_G = $FF_G + $l[5] + $l[9];
				$FF_A = $FF_A + $l[3] + $l[7];
				
				$FFPE_G = $FFPE_G + $l[14] + $l[18];
				$FFPE_A = $FFPE_A + $l[12] + $l[16];
			}
		}
		close(INTERSECTION);

		$Tv_C_FF = $FF_T / $FF_C ;
		$Tv_C_FFPE = $FFPE_T / $FFPE_C ;
		$Tv_G_FF = $FF_A / $FF_G ;
		$Tv_G_FFPE = $FFPE_A / $FFPE_G ;		

		print "$sample[$o]\t$FF_C\t$FF_T\t$Tv_C_FF\t$FFPE_C\t$FFPE_T\t$Tv_C_FFPE\t$FF_G\t$FF_A\t$Tv_G_FF\t$FFPE_G\t$FFPE_A\t$Tv_G_FFPE\n";
		print WRITE_TV "$sample[$o]\t$FF_C\t$FF_T\t$Tv_C_FF\t$FFPE_C\t$FFPE_T\t$Tv_C_FFPE\t$FF_G\t$FF_A\t$Tv_G_FF\t$FFPE_G\t$FFPE_A\t$Tv_G_FFPE\n";
}
close (WRITE_TV);

print "Done. \n";
print strftime "%Y-%m-%d %H:%M:%S", localtime(time);
print "\n";

