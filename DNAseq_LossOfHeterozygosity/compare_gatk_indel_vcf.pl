#!/usr/bin/perl -w
use strict;

print STDERR "./compare_gatk_indel_vcf.pl \n";

#my $threshold = "_pass";
my $threshold = "_q30cov10";
my $covThr = 10;

my $method = "bwamem_default";
#my @NORMAL = ("N_KK","N_JM","N_JM");
#my @TUMOR = ("T_KK","T_JM1","T_JM2");
my @NORMAL = ("N_JM","N_JM");
my @TUMOR = ("T_JM1","T_JM2");

for ( my $k = 0; $k <=$#NORMAL ; $k++) {

my @file;
$file[0] = "</project/results/Godwin/Exome/$TUMOR[$k]/genotyper/$TUMOR[$k].$method.indel.recal$threshold.vcf";
$file[1] = "</project/results/Godwin/Exome/$NORMAL[$k]/genotyper/$NORMAL[$k].$method.indel.recal$threshold.vcf";
$file[2] = ">/project/results/Godwin/Exome/InDel_comparison/$TUMOR[$k]_$NORMAL[$k].$method$threshold\_con";
$file[3] = ">/project/results/Godwin/Exome/InDel_comparison/$TUMOR[$k]_$NORMAL[$k].$method$threshold\_dis";
$file[4] = "</project/results/Godwin/Exome/$NORMAL[$k]/$NORMAL[$k].$method.cns.pileup";
$file[5] = "</project/results/Godwin/Exome/$TUMOR[$k]/$TUMOR[$k].$method.cns.pileup";

my %TUMOR_indel; my %w_tumor;
my $it1 = 0;
print "opens $file[0]...\n";
open TUMOR, $file[0] or die "no $file[0]";
while(<TUMOR>){
	chomp();
	my @l = split(/\t/,$_);			# array @l w/o tab
	my @first = split(undef,$l[0]);
	if ($first[0] eq "#" || $first[0] eq "chr") {next;}
	my $chr = $l[0];
	my $pos = $l[1];	# hashkey: chr \t pos \t ref
        my $ref = $l[3];
	my $indel = $l[4];
	
	$TUMOR_indel{$chr}{$pos} = $indel;
	$w_tumor{$chr}{$pos} = join("\t",$ref,$indel,$l[5],$l[2],$l[7],$l[8],$l[9]);	# info
	++$it1;
}
close(TUMOR);

my %NORMAL_indel; my %w_normal;
my $it2 = 0;
print "opens $file[1]...\n";
open NORMAL, $file[1] or die "no $file[1]";
while(<NORMAL>){
	chomp();
	my @l = split(/\t/,$_);			# array @l w/o tab
	my @first = split(undef,$l[0]);
	if ($first[0] eq '#' || $first[0] eq 'chr') {next;}
	# CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  bar
	my $chr = $l[0];
	my $pos = $l[1];	# hashkey: chr \t pos \t ref
        my $ref = $l[3];
	my $indel = $l[4];
	
	$NORMAL_indel{$chr}{$pos} = $indel;
	$w_normal{$chr}{$pos} = join("\t",$ref,$indel,$l[5],$l[2],$l[7],$l[8],$l[9]);	# info
	++$it2;
}
close(NORMAL);

my %compareTUMOR;
my %compareNORMAL;

print "opens $file[5]...\n";
open PILEUP_TUMOR, $file[5] or die "no $file[5]";
while(<PILEUP_TUMOR>){
	chomp();
	my @l = split(/\t/,$_);			# array @l w/o tab
	my $chr = $l[0];
	my $pos = $l[1];
	my $covTUMOR = $l[7];
		
	if ( exists $NORMAL_indel{$chr}{$pos}  && $covTUMOR >= $covThr ) {
		my $ref = $l[2];
		my $indel = $l[3];
		$compareTUMOR{$chr}{$pos} = $ref."\t".$indel."\t".$covTUMOR;
	}
}
close(PILEUP_TUMOR);

print "opens $file[4]...\n";
open PILEUP_NORMAL, $file[4] or die "no $file[4]";
while(<PILEUP_NORMAL>){
	chomp();
	my @l = split(/\t/,$_);			# array @l w/o tab
	my $chr = $l[0];
	my $pos = $l[1];
	my $covNORMAL = $l[7];
		
	if ( exists $TUMOR_indel{$chr}{$pos}  && $covNORMAL >= $covThr ) {
		my $ref = $l[2];
		my $indel = $l[3];
		$compareNORMAL{$chr}{$pos} = $ref."\t".$indel."\t".$covNORMAL;
	}
}
close(PILEUP_NORMAL);

my %con, my %dis, my %dis1, my %dis2;
my $con = 0; my $dis1 = 0, my $dis2 = 0; my $posit = 0; my $dis = 0;

foreach my $chr ( sort keys %TUMOR_indel){  	# $chr: chr \t pos \t ref	
	if (defined $NORMAL_indel{$chr}){
		foreach my $pos (sort keys %{$TUMOR_indel{$chr}}){
			if (defined $NORMAL_indel{$chr}{$pos}) {
			        ++$posit;
			        if ($TUMOR_indel{$chr}{$pos} eq $NORMAL_indel{$chr}{$pos}){
				        $con{$chr}{$pos} = join ("\t",$w_tumor{$chr}{$pos},$w_normal{$chr}{$pos});
				        ++$con;
				}
			        else {
				        $dis{$chr}{$pos} = join ("\t",$w_tumor{$chr}{$pos},$w_normal{$chr}{$pos});
				        ++$dis;
			        }
		        }
		}
        }
}

foreach my $chr ( sort keys %TUMOR_indel){  	# $chr: chr \t pos \t ref
	foreach my $pos (sort keys %{$TUMOR_indel{$chr}}){
		if  (!exists $dis{$chr}{$pos} && !exists $con{$chr}{$pos} && exists $compareNORMAL{$chr}{$pos} ) {
			$dis1{$chr}{$pos} = join ("\t",$w_tumor{$chr}{$pos},$compareNORMAL{$chr}{$pos});
			++$dis1;
		}
	}
}
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  bar
my $sum1 = $con + $dis1;
my $sum2 = $con + $dis2;
my $div1 = 0; my $div2 = 0;
if ( $posit == 0 ) {
	print "ERROR!! NO POSITIONS IN COMMON!!!\n";
} else {
	$div1 = $con/$posit*100;
	$div2 = $dis1/$posit*100;
}
open WRITE_CON, $file[2] or die;
print "$con concordances out of $posit same positions ($it1 - $TUMOR[$k], $it2 $NORMAL[$k])\n";#This is $div1 percent!\n";
print WRITE_CON "chr\tpos\tref\t$TUMOR[$k]\tdbSNP\tfilter\tinfo\tformat\tbar\tref\t$NORMAL[$k]\tdbSNP\tfilter\tinfo\tformat\tbar\n";
foreach my $chr ( sort keys %con ){
	foreach my $pos (sort keys %{$con{$chr}} ){
		print WRITE_CON $chr."\t".$pos."\t".$con{$chr}{$pos}."\n";
	}
}
close (WRITE_CON);

open WRITE_DIS, $file[3] or die;
print "$dis discordances and $dis1 positions in $TUMOR[$k] that were not called in $NORMAL[$k] but at least position has to be there ($it1 - $TUMOR[$k], $it2 $NORMAL[$k])\nThis is $div1/$div2 percent!\n";
print WRITE_DIS "chr\tpos\tref\t$TUMOR[$k]\tQUAL\tdbSNP\tinfo\tformat\tbar\tref\t$NORMAL[$k]\tQUAL\tdbSNP\tinfo\tformat\tbar\n";

foreach my $chr ( sort keys %dis ){
	foreach my $pos (sort keys %{$dis{$chr}} ){
		print WRITE_DIS $chr."\t".$pos."\t".$dis{$chr}{$pos}."\n";
	}
} print WRITE_DIS "\n";
foreach my $chr ( sort keys %dis1 ){
	foreach my $pos (sort keys %{$dis1{$chr}} ){
		print WRITE_DIS $chr."\t".$pos."\t".$dis1{$chr}{$pos}."\n";
	}
} print WRITE_DIS "\n";
foreach my $chr ( sort keys %dis2 ){
	foreach my $pos (sort keys %{$dis2{$chr}} ){
		print WRITE_DIS $chr."\t".$pos."\t".$dis2{$chr}{$pos}."\n";
	}
}
close(WRITE_DIS);
print "Written in $file[2] and $file[3]\n\n";
}
print "Done. \n";
