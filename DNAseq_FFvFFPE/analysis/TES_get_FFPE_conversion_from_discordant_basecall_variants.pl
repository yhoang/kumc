#!/usr/bin/perl -w
# Author Yen Hoang, M.Sc.
# 2013

use strict;
use POSIX qw/strftime/;

# Conversion rate of SNVs in targeted intersected discordant positions in FFPE Files
# i.e. A>T 30%, A>G 30%, A>C 40%
# Looking at base call only (not heterozygosity)

print strftime "%Y-%m-%d %H:%M:%S", localtime(time);
print "\n";
print "./TES_get_FFPE_conversion_from_discordant_basecall_variants.pl \n";

our $p = ("q43cov13");
our @sample = ("2474","2561","2640","2685","2938","3050","3356","4079","4191","CN","CT","DN","DT");

my @nuc = ('A','T','G','C');
my @distribution;
for ( my $i = 0; $i <=9; $i++) {
        for (my $j = 0; $j <=16; $j++) {
                $distribution[$i][$j] = 0;
        }
}

my $title1 = "conversion from nuc1\>nuc2 (actual count)";
my $title2 = "\n\nconversion from nuc1\>nuc2 (frequency)";
for ( my $k = 0; $k <=$#sample; $k++) {
	$title1 = $title1."\t$sample[$k]";
	$title2 = $title2."\t$sample[$k]";
 	my $dis_file = "</project/results/me/TES/SNV_comparison/genotyper/with_intersection/base_call/$sample[$k]_FF_FFPE_${p}_targeted_dis";
        
        my @conversion;
        my $it = 0;
        
        print "opens $dis_file...\n";
        open DISCORDANT_FILE, $dis_file or die "no DISCORDANT_FILE:$dis_file";
        while(<DISCORDANT_FILE>){
        	chomp();
	        my @l = split(/\t/,$_);			# array @l w/o tab
	        if ($l[0] eq "chr") {next;}
	        if ($l[0] eq "") {next;}
                my $ref = $l[2];
                my $varFF = $l[3];
                my $varFFPE = $l[14];
                if ($l[15] eq "\"het\"" || $l[15] eq "\"hom\"") {
                        $conversion[$it][0] = $ref;
                        $conversion[$it][1] = $varFFPE;
                } else {
                        $conversion[$it][0] = $varFF;
                        $conversion[$it][1] = $ref;
                }
                ++$it;
        }
        close(DISCORDANT_FILE);

        for (my $j = 0; $j <=3; $j++) { # for nuc1
                for (my $i = 0; $i <= $it; $i++) {      #for every discordant call
                        if ($conversion[$i][0] eq $nuc[$j]) {   
                                for (my $l = 0; $l <=3; $l++) { #for conversion in nuc2
                                        if ($conversion[$i][1] eq $nuc[$l]) {
                                                ++$distribution[$k][4*$j+$l];
                                        }
                                }
                        }
                }
        }
        $distribution[$k][16] = $it;
}
$title1 = $title1."\tmean\n";
$title2 = $title2."\tmean\n";
### print actual count
print "$title1";
for (my $i = 0; $i <= 3; $i++) {        #nuc1
        for (my $j = 0; $j <= 3; $j++) {#nuc2
                my $mean = 0;
                print "$nuc[$i]\>$nuc[$j]\t";
                for (my $s = 0; $s <= $#sample; $s++) {     #for sample number
                        print "$distribution[$s][4*$i+$j]\t";
                        $mean +=$distribution[$s][4*$i+$j];
                }
                $mean = $mean/($#sample+1);
                print "$mean\n";
        }
}
print "count";
for (my $s = 0; $s <= $#sample; $s++) {
        print "\t$distribution[$s][16]";
}


### print frequncy
for (my $s = 0; $s <= $#sample; $s++) {
        for (my $i = 0; $i <=15; $i++) {
                        $distribution[$s][$i] = $distribution[$s][$i]/$distribution[$s][16];
        }
}
print "$title2";
for (my $i = 0; $i <= 3; $i++) {
        for (my $j = 0; $j <= 3; $j++) {
                my $mean = 0;
                print "$nuc[$i]\>$nuc[$j]\t";
                for (my $s = 0; $s <= $#sample; $s++) {
                        print "$distribution[$s][4*$i+$j]\t";
                        $mean +=$distribution[$s][4*$i+$j];
                }
                $mean = $mean/($#sample+1);
                print "$mean\n";
        }
}
print "count";
for (my $s = 0; $s <= $#sample; $s++) {
        print "\t$distribution[$s][16]";
}

print "\n\nDone. \n";
print strftime "%Y-%m-%d %H:%M:%S", localtime(time);
print "\n";
