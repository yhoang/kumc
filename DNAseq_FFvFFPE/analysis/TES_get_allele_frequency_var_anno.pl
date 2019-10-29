#!/usr/bin/perl -w
# Author Yen Hoang, M.Sc.
# 2013

use strict;
use POSIX qw/strftime/;

# Compares allele frequency of SNVs calculated from pileups between FF and FFPE.

print strftime "%Y-%m-%d %H:%M:%S\n", localtime(time);
print "\n";
print STDERR "./TES_get_allele_frequency_var_anno.pl stepsize [100]\n";

my $stepsize = $ARGV[0];
#our @sample = ("2474","2561","2640","2685","2938","3050","3356","4079","4191","CN","CT","DN","DT");
our @method = ("q43cov13");

### Initiation
my @AF_FF, my @AF_FFPE;
my @AF_FF2, my @AF_FFPE2;
my $title = "AF\tcountGATK\tcountPILEUP";
my $title2 = "GATK\tPILEUP";

 
for (my $i = 0; $i <= 8000; $i++) {
        for (my $j = 0; $j <=1; $j++) {
                $AF_FF[$i][$j] = 0;
                $AF_FFPE[$i][$j] = 0;
                $AF_FF2[$i][$j] = 0;
                $AF_FFPE2[$i][$j] = 0;
        }
}

for (my $o = 0; $o <= $#sample; $o++) {
	my $max = 0; my $max2 = 0; my $it1 = 0; my $it2 = 0; my $it11 = 0; my $it22 = 0;
	print strftime "%Y-%m-%d %H:%M:%S\n", localtime(time);

	foreach my $p (@method) {
		my $FF = "</project/results/me/TES/$sample[$o]_FF/genotyper/anno/$sample[$o]_FF.bwa.ug.snp_${p}.genome_summary.csv";
		my $FFPE = "</project/results/me/TES/$sample[$o]_FFPE/genotyper/anno/$sample[$o]_FFPE.bwa.ug.snp_${p}.genome_summary.csv";
		my $distribution = "</project/results/me/TES/pileup_distribution/$sample[$o]_FF_FFPE_intersect.distribution";
		my $write3 = ">/project/results/me/TES/plots/alleleFreq/genotyper/density_data/$sample[$o]_AF2_FF_${p}\_$stepsize";
		my $write4 = ">/project/results/me/TES/plots/alleleFreq/genotyper/density_data/$sample[$o]_AF2_FFPE_${p}\_$stepsize";
	
		my %SNV_FF; 
		print "opens $FF...\n";
		open F1, $FF or die "no $FF"; 
		while(<F1>){
			chomp();
			my @l = split(/\t/,$_);			# array @l w/o tab
			if ($l[0] eq 'Func') {next;}
			my $chr = $l[21];
			my $pos = $l[22];
			my $var = $l[25];
			$AF_FF2[$it1][2] = "$chr\t$pos\t$var";
			if ($var eq "T") {$var = 1;}
			elsif ($var eq "G") {$var = 2;}
			elsif ($var eq "C") {$var = 3;}
			else {$var = 0;}
			$SNV_FF{$chr}{$pos} = $var;
			
			my $freq = substr $l[31], 1, -1;
			$freq = $freq * $stepsize;
			++$AF_FF[int($freq)][0];
			
			$AF_FF2[$it1][0] = int($freq);
			++$it1;
		}
		close(F1);
		
		
		my %SNV_FFPE; 
		print "opens $FFPE...\n";
		open F2, $FFPE or die "no $FFPE"; 
		while(<F2>){
			chomp();
			my @l = split(/\t/,$_);			# array @l w/o tab
			if ($l[0] eq 'Func') {next;}
			my $chr = $l[21];
			my $pos = $l[22];
			my $var = $l[25];
			$AF_FFPE2[$it1][2] = "$chr\t$pos\t$var";
			if ($var eq "T") {$var = 1;}
			elsif ($var eq "G") {$var = 2;}
			elsif ($var eq "C") {$var = 3;}
			else {$var = 0;}
			$SNV_FFPE{$chr}{$pos} = $var;
			
			my $freq = substr $l[31], 1, -1;
			$freq = $freq * $stepsize;
			++$AF_FFPE[int($freq)][0];
			
			$AF_FFPE2[$it2][0] = int($freq);
			$AF_FFPE2[$it2][2] = "$chr\t$pos\t$var";
			++$it2;
		}
		close(F2);
		
			       
		print "opens $distribution...\n";
		open DISTR, $distribution or die "no $distribution";
		while(<DISTR>) {
		        chomp();
			my @l = split(/\t/,$_);			# array @l w/o tab
			if ($l[0] eq "chr") {next;}
			my $chr = $l[0];
			my $pos = $l[1];
			
			my $covFF = $l[3] + $l[4] + $l[5] + $l[6] + $l[7] + $l[8] + $l[9] + $l[10];
			my $covFFPE = $l[12] + $l[13] + $l[14] + $l[15] + $l[16] + $l[17] + $l[18] + $l[19];
			
			if (exists $SNV_FF{$chr}{$pos} && $covFF != 0) {
			        my $var = $l[3+$SNV_FF{$chr}{$pos}] + $l[7+$SNV_FF{$chr}{$pos}];
			        my $freq = $var/$covFF*$stepsize;
			        ++$AF_FF[int($freq)][1];
			        
			        $AF_FF2[$it11][1] = int($freq);
			        ++$it11;
			 #       print "FF\t$var\t$covFF\t$SNV_FF{$chr}{$pos}\t$freq\t$AF_FF[int($freq)][1]\t$chr\t$pos\n";
			        
			}
			if (exists $SNV_FFPE{$chr}{$pos} && $covFFPE != 0) {
			        my $var = $l[12+$SNV_FFPE{$chr}{$pos}] + $l[16+$SNV_FFPE{$chr}{$pos}];
			        my $freq = $var/$covFFPE*$stepsize;
			        ++$AF_FFPE[int($freq)][1];
			        
			        $AF_FFPE2[$it22][1] = int($freq);
			        ++$it22;
			}
		}
		close(DISTR);
		
		$max = $it1; $max2 = $it2;
		#if ($it1 > $it2) {$max = $it1;} else {$max = $it2;}
		#if ($it11 > $it22 ) {$max2 = $it11;} else {$max = $it22;}
		
		open WRITE, $write1 or die;
		print WRITE "$title\n";
		for (my $i = 0; $i <= 100; $i+=(100/$stepsize)) {
			print WRITE "$i\t$AF_FF[$i][0]\t$AF_FF[$i][1]\t$AF_FF[$i][2]\n";
		}
		close(WRITE);
		print STDERR "Written in $write1. \n";
		
		open WRITE2, $write2 or die;
		print WRITE2 "$title\n";
		for (my $i = 0; $i <= 100; $i+=(100/$stepsize)) {
			print WRITE2 "$i\t$AF_FFPE[$i][0]\t$AF_FFPE[$i][1]\t$AF_FFPE[$i][2]\n";
		}
		close(WRITE2);
		print  STDERR "Written in $write2.\n";

		open WRITE3, $write3 or die;
		print WRITE3 "$title2\n";
		for (my $i = 0; $i < $max; $i++) {
			my $help1 = $AF_FF2[$i][0]/($stepsize/100);
			my $help2 = $AF_FF2[$i][1]/($stepsize/100);
			print WRITE3 "$help1\t$help2\n";#\t$AF_FF2[$i][2]\n";
		}
		close(WRITE3);
		print STDERR "Written in $write3. \n";
		
		open WRITE4, $write4 or die;
		print WRITE4 "$title2\n";
		for (my $i = 0; $i < $max2; $i++) {
			my $help1 = $AF_FFPE2[$i][0]/($stepsize/100);
			my $help2 = $AF_FFPE2[$i][1]/($stepsize/100);
			print WRITE4 "$help1\t$help2\n";#\t$AF_FFPE2[$i][2]\n"
		}
		close(WRITE4);
		print  STDERR "Written in $write4.\n\n";
	}
}

print "Done.\n"; 
print strftime "%Y-%m-%d %H:%M:%S", localtime(time);
print "\n";

