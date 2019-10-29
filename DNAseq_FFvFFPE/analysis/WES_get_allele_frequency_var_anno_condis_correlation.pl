#!/usr/bin/perl -w
# Author Yen Hoang, M.Sc.
# 2013

use strict;
use POSIX qw/strftime/;

# Compares allele frequency of SNVs calculated from intersected, exonic intersected, and targeted intersected pileups between FF and FFPE, respectively.
# With differentiation concordance and discordance.

print strftime "%Y-%m-%d %H:%M:%S", localtime(time);
print "\n";
print STDERR "./WES_get_allele_frequency_var_anno_condis_correlation.pl\n";

#my $threshold = "";
our @sample = ("2474","2561","2640","2685","2938","3050","3356","4079","4191","CN","CT","DN","DT");
our @method = ("q43cov13");
our @range = ("total","exonic","targeted");

### Initiation
my $title2 = "chr\tpos\tvar_FF\tAF_FF\tvar_FFPE\tAF_FFPE";

foreach my $o (@sample) {
	my $it1 = 0; my $it2 = 0;
	
	foreach my $q (@range) {
	foreach my $p (@method) {
        my $CON_FILE = "</project/results/me/WES/SNV_comparison/genotyper/with_intersection/stats/exonic/${o}_FFex_FFPEex_${p}_con.tsv";
	my $DIS_FILE = "</project/results/me/WES/SNV_comparison/genotyper/with_intersection/stats/exonic/${o}_FFex_FFPEex_${p}_dis.tsv";
        my $WRITE_AF_CON = ">/project/results/me/WES/plots/alleleFreq/genotyper/correlation_data/AF_${o}_FFex_FFPEex_${p}_${q}_con_correlation";
        my $WRITE_AF_DIS = ">/project/results/me/WES/plots/alleleFreq/genotyper/correlation_data/AF_${o}_FFex_FFPEex_${p}_${q}_dis_correlation";
	
	open WRITE_CON, $WRITE_AF_CON or die;
	print WRITE_CON "$title2\n";
        print "opens $CON_FILE...\n";
        open CON_FILE, $CON_FILE or die "no $CON_FILE"; 
        while(<CON_FILE>){
        	chomp();
	        my @l = split(/\t/,$_);			# array @l w/o tab
	        if ($l[0] eq 'chr') {next;}
	        my ($chr, $pos, $var_FF, $var_FFPE, $cov_FF, $cov_FFPE, $varcov_FF, $varcov_FFPE, $freq_FF, $freq_FFPE) = qw/0 0 0 0 0 0 0 0 0 0/;
	        $chr = $l[0];
	        $pos = $l[1];
	        
	        $var_FF = $l[3];
	        if ($var_FF eq "T") {$var_FF = 1;}
	        elsif ($var_FF eq "G") {$var_FF = 2;}
	        elsif ($var_FF eq "C") {$var_FF = 3;}
	        else {$var_FF = 0;}
	        
		$var_FFPE = $l[10];
	        if ($var_FFPE eq "T") {$var_FFPE = 1;}
	        elsif ($var_FFPE eq "G") {$var_FFPE = 2;}
	        elsif ($var_FFPE eq "C") {$var_FFPE = 3;}
	        else {$var_FFPE = 0;}
	        
	        $cov_FF = $l[12] + $l[13] + $l[14] + $l[15] + $l[16] + $l[17] + $l[18] + $l[19];
	        $cov_FFPE = $l[21] + $l[22] + $l[23] + $l[24] + $l[25] + $l[26] + $l[27] + $l[28];
	        
	        $varcov_FF = $l[12+$var_FF] + $l[16+$var_FF];
	        $varcov_FFPE = $l[21+$var_FFPE] + $l[25+$var_FFPE];
	        
	        if ($cov_FF != 0 && $cov_FFPE != 0 ) {
	        $freq_FF = $varcov_FF / $cov_FF;
	        $freq_FFPE = $varcov_FFPE / $cov_FFPE;
	        
	   #     print "$chr\t$pos\t$l[3]\t$var_FF\t$cov_FF\t$varcov_FF\t$l[10]\t$var_FFPE\t$cov_FFPE\t$varcov_FFPE\n";
	        
		print WRITE_CON "$chr\t$pos\t$l[3]\t$freq_FF\t$l[10]\t$freq_FFPE\n";
		}
	        ++$it1;
	        
        }
        close(CON_FILE);
        close(WRITE_CON);
        
	open WRITE_DIS, $WRITE_AF_DIS or die;
	print WRITE_DIS "$title2\n";
        print "opens $DIS_FILE...\n";
        open DIS_FILE, $DIS_FILE or die "no $DIS_FILE"; 
        while(<DIS_FILE>){
        	chomp();
	        my @l = split(/\t/,$_);			# array @l w/o tab
	        if ($l[0] eq 'chr') {next;}
	        my ($chr, $pos, $var_FF, $var_FFPE, $cov_FF, $cov_FFPE, $varcov_FF, $varcov_FFPE, $freq_FF, $freq_FFPE) = qw/0 0 0 0 0 0 0 0 0 0/;
	        
	        $chr = $l[0];
	        $pos = $l[1];
	        $var_FF = $l[3];
		$var_FFPE = $l[10];
		
	        if ($var_FF eq "T") {$var_FF = 1;}
	        elsif ($var_FF eq "G") {$var_FF = 2;}
	        elsif ($var_FF eq "C") {$var_FF = 3;}
	        else {$var_FF = 0;}
		
	        if ($var_FFPE eq "T") {$var_FFPE = 1;}
	        elsif ($var_FFPE eq "G") {$var_FFPE = 2;}
	        elsif ($var_FFPE eq "C") {$var_FFPE = 3;}
	        else {$var_FFPE = 0;}
	        
	        $cov_FF = $l[12] + $l[13] + $l[14] + $l[15] + $l[16] + $l[17] + $l[18] + $l[19];
	        $cov_FFPE = $l[21] + $l[22] + $l[23] + $l[24] + $l[25] + $l[26] + $l[27] + $l[28];
	        
	        if ($cov_FF != 0 && $cov_FFPE != 0 ) {
	        	if ( $l[3] ne "" && $l[10] ne "") {
				$varcov_FF = $l[12+$var_FF] + $l[16+$var_FF];
				$varcov_FFPE = $l[21+$var_FFPE] + $l[25+$var_FFPE];
			} elsif ($l[3] eq "") {
				$varcov_FF = $l[12+$var_FFPE] + $l[16+$var_FFPE];
				$varcov_FFPE = $l[21+$var_FFPE] + $l[25+$var_FFPE];
			} else {
				$varcov_FF = $l[12+$var_FF] + $l[16+$var_FF];
				$varcov_FFPE = $l[21+$var_FF] + $l[25+$var_FF];
			}
			
			$freq_FF = $varcov_FF / $cov_FF;
			$freq_FFPE = $varcov_FFPE / $cov_FFPE;
			
			print WRITE_DIS "$chr\t$pos\t$l[3]\t$freq_FF\t$l[10]\t$freq_FFPE\n";
		}
	        ++$it2;
        }
        close(DIS_FILE);
        close(WRITE_DIS);
       
 
print STDERR "Written in $WRITE_AF_CON and $WRITE_AF_DIS. \n\n";
        
}}}

print "\nDone.\n"; 
print strftime "%Y-%m-%d %H:%M:%S", localtime(time);
print "\n";

