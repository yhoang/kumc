#!/usr/bin/perl -w
#./compare_normal_tumor.pl <ex> optional
use POSIX qw/strftime/;
use strict;

print strftime "%Y-%m-%d %H:%M:%S", localtime(time);
print STDERR "\n./get_compound_var.pl  \n--> Filter exonic&&nonsynonymous&&nottolerant and splicing&&noannotation\n--> Gene X pos X: N(var) and T(var) same, but pos Y: N(het) and T(hom) differ!\n";

#my @aligner = ("bwa","clc");
my @method = ("gatk","sam");
my @option = ("N_KK","N_JM","N_JM","N_JMM","N_JMM","N_JMM");
my @option2 = ("T_KK","T_JM2","T_JM1","T_JM2","T_JM1","N_JM");

for (my $m = 0; $m <=$#method; $m++){
for (my $k = 0; $k <=$#option; $k++){

my $f1 = "</home/yhoang/Godwin/sift/new/".$option[$k]."_bwa_$method[$m].genome_summary.csv";
my $f2 = "</home/yhoang/Godwin/sift/new/".$option2[$k]."_bwa_$method[$m].genome_summary.csv";

print "\nComparing in $method[$m] $option[$k] with $option2[$k]..\n";

my $write1 = ">/home/yhoang/Godwin/compound_var/compound_var_bwa_$method[$m]\_$option[$k]\_$option2[$k].xls";

my %locN; my $locN = 0;
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
	if ( ($l[13] eq "D") || ($l[13] eq "NA" ) ) { ++$filter;}# print "$filter=$l[13]\n"; }
	if ( ($l[0] eq "\"splicing\"") && ($l[2] eq "") ) {$filter = 3;}# print "$filter=$l[0]\n";} 
	if ($filter < 3) {next;}
	my $gene = $l[1];
	my $pos = $l[21]."\t".$l[22]."\t".$l[24];
	my $var = $l[25]."\t".$l[26];
       	$locN{$gene}{$pos} = $var;
       	++$locN;
}
close(F1);
print "Found $locN positions in $option[$k]\_$method[$m] which where exonic-nonsyn or splicing-unknown.\n";

my $it2; my %con;
my %dis; my $itdis = 0; my $itcon = 0;
open F2,$f2 or die "No $f2!";
while(<F2>){
	chomp();
	my @l = split(/,/,$_);			# array @l w/o tab
	if ($l[0] eq "Func") {next;}
	my $gene = $l[1];
	my $pos = $l[21]."\t".$l[22]."\t".$l[24];
	my $var = $l[25]."\t".$l[26];
	if ( (exists $locN{$gene}{$pos}) && ($var eq $locN{$gene}{$pos}) ){
        # if gene and position and variant calling is in both normal and tumor the same
                $con{$gene} = $pos."\t".$locN{$gene}{$pos}."\t".$var."\t".$l[0]."\t".$l[13]."\t".$l[7]."\t".$l[8]."\t".$l[27]."\t".$l[28]."\t".$l[29];
                ++$itcon;
        }
        elsif ( (exists $locN{$gene}{$pos}) && ($var ne $locN{$gene}{$pos}) ){
                $dis{$gene} = $pos."\t".$locN{$gene}{$pos}."\t".$var."\t".$l[0]."\t".$l[13]."\t".$l[7]."\t".$l[8]."\t".$l[27]."\t".$l[28]."\t".$l[29];
                ++$itdis;
        }
       	++$it2;
}

close(F2);

print "Found $itcon concordant variant calling and $itdis discordant variant calling out of $it1(normal)/$it2(tumor).\n";
print "Looking for compound variants now... and write into $write1...\n";

my $compound = 0;
open WRITE1, $write1 or die;
if ($m == 0) { # for gatk
        print WRITE1 "gene\tchr\tpos\tref\tNormal\tvar type\tTumor\tvar type\tfunction\tSIFT\t1000g2012feb_ALL\tdbSNP135\tquality score\tread depth\tRMS mapping quality\n"; 
} else { print WRITE1 "gene\tchr\tpos\tref\tNormal\tvar type\tTumor\tvar type\tfunction\tSIFT\t1000g2012feb_ALL\tdbSNP135\tsnp quality\ttotal reads\treads with mutation\n"; }
foreach my $gene (sort keys %con) {
        if (exists $dis{$gene}) {
                print WRITE1 $gene."\t".$con{$gene}."\n";
                print WRITE1 $gene."\t".$dis{$gene}."\n";
                ++$compound;
        }
}
print "$compound compound variants found.\n";
}}

print strftime "\n%Y-%m-%d %H:%M:%S", localtime(time);
print "\nDone.\n";


