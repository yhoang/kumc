#!/usr/bin/perl -w
#./compare_normal_tumor.pl <ex> optional
use POSIX qw/strftime/;
use strict;

print strftime "%Y-%m-%d %H:%M:%S", localtime(time);
print STDERR "\n./compare_het_hom_filtered.pl  \n--> Filter (exonic && nonsynonymous && !tolerant) and (splicing && !annotation)\n--> Find heterozygous positions in normal which changed to homozygous in tumor!\n";

#my @aligner = ("bwa","clc");
my @method = ("gatk","sam");
my @option = ("KK","JM","JM","JMM","JMM");
my @option2 = ("KK","JM2","JM1","JM2","JM1");

for (my $m = 0; $m <=$#method; $m++){
for (my $k = 0; $k <=$#option; $k++){

my $f1 = "</home/yhoang/Godwin/sift/new/N_".$option[$k]."_bwa_$method[$m].genome_summary.csv";
my $f2 = "</home/yhoang/Godwin/sift/new/T_".$option2[$k]."_bwa_$method[$m].genome_summary.csv";

print "\nComparing in $method[$m] N_$option[$k] with T_$option2[$k]..\n";

my $write1 = ">/home/yhoang/Godwin/loss_of_heterozygosity/pLOH_bwa_$method[$m]\_N_$option[$k]\_T_$option2[$k].xls";
my $write2 = ">/home/yhoang/Godwin/loss_of_heterozygosity/het_bwa_$method[$m]\_N_$option[$k]\_T_$option2[$k].xls";

my %het; my $het = 0;
my $it1 = 0;
open F1, $f1 or die "No $f1!";
while(<F1>){
	chomp();
	my @l = split(/,/,$_);			# array @l w/o tab
	if ($l[0] eq "Func") {next;}
	++$it1;
	my $filter = 0;
#	print "l = $l[0]\n";
	if ( ($l[0] eq "\"exonic\"") || ($l[0] eq "\"exonic;splicing\"") ) {++$filter;}# print "$filter=$l[0]\t";  }
	if ( $l[2] eq "\"nonsynonymous SNV\"" ) { ++$filter;}# print "$filter=$l[2]\t$l[13] ?"; }
	if ( ($l[13] eq "D") || ($l[13] eq "NA") ) { ++$filter;}# print "$filter=$l[13]\n"; }
	if ( ($l[0] eq "\"splicing\"") && ($l[2] eq "") ) {$filter = 3;}# print "$filter=$l[0]\n";} 
	if ($filter < 3) {next;}
	my $pos = $l[21]."\t".$l[22]."\t".$l[24];
	my $code;
	if ( $l[26] eq "\"het\"") {
		# change code
		my @var_sorted = sort($l[24],$l[25]) ;  #ref and var
		if ( $var_sorted[0] eq "A") {
			if ($var_sorted[1] eq "C") { $code = "M"; }
			if ($var_sorted[1] eq "G") { $code = "R"; }
			if ($var_sorted[1] eq "T") { $code = "W"; }				
		}		
		if ( $var_sorted[0] eq "C") {
			if ($var_sorted[1] eq "G") { $code = "S"; }
			if ($var_sorted[1] eq "T") { $code = "Y"; }				
		}	
		if ( $var_sorted[0] eq "G") {
			if ($var_sorted[1] eq "T") { $code = "K"; }				
		}	
		$het{$pos}=$l[1]."\t".$code."\t".$l[26]."\t".$l[0]."\t".$l[13]."\t".$l[7]."\t".$l[8]."\t".$l[27]."\t".$l[28]."\t".$l[29];
       	        ++$het;
       	}
}
close(F1);
print "Found $het heterozygous positions in N_$option[$k]\_$method[$m] which where exonic-nonsyn or in splicing-unknown.\n";

my $it2; my %con;
my %dis; my $itdis = 0; my $itcon = 0;
open F2,$f2 or die "No $f2!";
while(<F2>){
	chomp();
	my @l = split(/,/,$_);			# array @l w/o tab
	if ($l[0] eq "Func") {next;}
	my $gene = $l[1];
	my $pos = $l[21]."\t".$l[22]."\t".$l[24];
	my $var = $l[25];
	if (exists $het{$pos} && $l[26] eq "\"hom\"") {	# heterogenizity status
		$dis{$pos}{$var} = $het{$pos}."\t".$var."\t".$l[26]."\t".$l[27]."\t".$l[28]."\t".$l[29];
		++$itdis;	
	}
	elsif (exists $het{$pos} && $l[26] eq "\"het\"") {	# heterogenizity status
		$con{$pos}{$var} = $het{$pos}."\t".$var."\t".$l[26]."\t".$l[27]."\t".$l[28]."\t".$l[29];
		++$itcon;	
	}
	++$it2;
}
close(F2);


print "Found $itdis positions that lost heterozygosity and $itcon still the same out of $it1(normal)/$it2(tumor).\n";

open WRITE1, $write1 or die; # pLOH
if ($m == 0) { # for gatk
        print WRITE1 "gene\tchr\tpos\tref\tvarNormal\tvar type\tfunction\tSIFT\t1000g2012feb_ALL\tdbSNP135\tquality scoree\tread depth\tRMS mapping quality\tvarTumor\tvar type\tquality score\tread depth\tRMS mapping quality\n"; 
} else { print WRITE1 "gene\tchr\tpos\tref\tvarNormal\tvar type\tfunction\tSIFT\t1000g2012feb_ALL\tdbSNP135\tsnp quality\ttotal reads\treads with mutation\tvarTumor\tvar type\tsnp quality\ttotal reads\treads with mutation\n"; }


foreach my $pos ( sort keys %dis ){
	foreach my $var (sort keys %{$dis{$pos}} ){
		print WRITE1 $pos."\t".$dis{$pos}{$var}."\n";
	}
}
close(WRITE1);

open WRITE2, $write2 or die;    # stayed heterozygous
if ($m == 0) { # for gatk
        print WRITE2 "gene\tchr\tpos\tref\tvarNormal\tvar type\tfunction\tSIFT\t1000g2012feb_ALL\tdbSNP135\tquality score\tread depth\tRMS mapping quality\n"; 
} else { print WRITE2 "gene\tchr\tpos\tref\tvarNormal\tvar type\tfunction\tSIFT\t1000g2012feb_ALL\tdbSNP135\tsnp quality\ttotal reads\treads with mutation\n"; }
foreach my $pos ( sort keys %con ){
	print WRITE2 $pos."\t".$het{$pos}."\n";
}
close(WRITE2);
print "Written in $write1.\nWritten in $write2.\n";

}}

print strftime "\n%Y-%m-%d %H:%M:%S", localtime(time);
print "\nDone.\n";
