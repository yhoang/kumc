#!/usr/bin/perl -w
# Author Yen Hoang, M.Sc.
# 2013

use strict;
use POSIX qw/strftime/;

# This script counts the coverage of each nucleotide on each intersected positions and compares between FF and FFPE.
# Separate count of concordant and discordant positions

print strftime "%Y-%m-%d %H:%M:%S", localtime(time);
print "\n";
print STDERR "./WES_cns_concordance_over_coverage_in_intersection.pl \n";

my $mincov = 250;

our @sample = ("2474","2561","2640","2685","2938","3050","3356","4079","4191","CN","CT","DN","DT");
our @samplesize = (491266053,520469241,899117689,759536886,741913358,1021718843,412412592,672208956,530050859,656871468,614447869,187309470,189583989); 
our ($FF, $FFPE, $write_total, $write_targeted, $targeted);
our (@total_covFF, @total_covFFPE);

print STDERR "sample\tpos_FF\tpos_FFPE\tintersect_total\tcon_total\tdis_total\n";
for (my $o = 0; $o <= $#sample; $o++) {
	#WES
    $FF = "</project/results/me/WES/$sample[$o]_FFex/$sample[$o]_FFex.cns.starless.pileup";	
    $FFPE = "</project/results/me/WES/$sample[$o]_FFPEex/$sample[$o]_FFPEex.cns.starless.pileup";	
    $write_total = ">/project/results/me/WES/cns_concordance/$sample[$o]_FFex_FFPEex_cns_cov";

	my ($it_FF, $it_FFPE, $it_same, $itcon, $itdis, $pointerPos, $pointerHelp, $loop, $it_loop) = qw/0 0 0 0 0 0 0 0 0 0 0 0/;
	
	for (my $i = 0; $i <= $mincov; $i++) {       
        	$total_covFF[$i][0] = 0;     # amount of concordance FF
        	$total_covFF[$i][1] = 0;     # amount of discordance FF
        	$total_covFFPE[$i][0] = 0;   # amount of concordance FFPE
        	$total_covFFPE[$i][1] = 0;   # amount of discordance FFPE
	}
       	
    $loop = 1;
	
	do {
		my (%FF, %help_FF);
		++$it_loop;
		$pointerHelp = 0;
		
        	# print STDERR "opens $FF...\n";
        	open FF_FILE, $FF or die "no $FF";
		seek FF_FILE, $pointerPos, 0;
	        while(<FF_FILE>){
        		chomp();
		        my @l = split(/\t/,$_);	
		        my $chr = $l[0];
		        my $pos = $l[1];
		        my $ref = $l[2];
		        my $var = $l[3];
			    my $mapQ = $l[6];
		        my $cov = $l[7];
			
			    $mapQ == 0 and next;
			
		        $FF{$chr}{$pos} = $var; 
		        $help_FF{$chr}{$pos} = $cov;
		        
		        ++$pointerHelp;
		        ++$it_FF;
		        
		        if ( $pointerHelp >= 20000000 ) {
				    $pointerPos = tell FF_FILE;
			        #	print "Buffer on line $it_FF/$samplesize[$o] (byte $pointerPos) - Break!\n";
				    last;
			} 
			
			if ( $it_FF >= $samplesize[$o] ) {
				$loop = 0;
				$pointerPos = tell FF_FILE;
			#	print "Buffer on LAST line $it_FF/$samplesize[$o] (byte $pointerPos)!\n";
				last;
			}
        	}       
        	close(FF_FILE);
		
		$it_FFPE = 0;
	        # print STDERR "opens $FFPE...\n";
        	open FFPE_FILE, $FFPE or die "no $FFPE";
        	while(<FFPE_FILE>){
		        chomp();
		        my @l = split(/\t/,$_);	
		        my $chr = $l[0];
		        my $pos = $l[1];
		        my $ref = $l[2];
		        my $var = $l[3];
			my $mapQ = $l[6];
		        my $cov = $l[7];
			
			$mapQ == 0 and next;
			
		        
			++$it_FFPE;
			
		        if ( exists $FF{$chr}{$pos} ) {
				++$it_same;
		        	if ($FF{$chr}{$pos} eq $var) {
			                for (my $i = 0; $i <= $mincov; $i += 2 ) {
				                if ( $i > $help_FF{$chr}{$pos} ) { ++$total_covFF[$i][0]; }
				        }
				        for (my $i = 1; $i <= $mincov; $i += 2 ) {
				        	if ( $i > $cov ) { ++$total_covFFPE[$i][0]; }
				        }
			                ++$itcon;
				} else {
				        for (my $i = 0; $i <= $mincov; $i += 2 ) {
				                if ( $i > $help_FF{$chr}{$pos} ) { ++$total_covFF[$i][1]; }
				        }
				        for (my $i = 1; $i <= $mincov; $i += 2 ) {
				        	if ( $i > $cov ) { ++$total_covFFPE[$i][1]; }
				        }
				        ++$itdis;
				}
			}
		}
        	close(FFPE_FILE);
	} while ( $loop == 1 );
	
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

	print STDERR "$sample[$o]\t$it_FF\t$it_FFPE\t$it_same\t$itcon\t$itdis\n";
}

print STDERR "\nDone.\n";
print strftime "%Y-%m-%d %H:%M:%S", localtime(time);
print "\n";
