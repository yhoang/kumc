#!/usr/bin/perl -w
#./compare_tumor_normal.pl <ex> optional
use POSIX qw/strftime/;
use strict;

print strftime "%Y-%m-%d %H:%M:%S", localtime(time);
print STDERR "./get_trueSNVs.pl\n--> Compare tumor and normal Variant Calls and get true Variants compared with consensus.\n";

#my @aligner = ("bwa","clc");
my @method = ("gatk","sam");
my @option2 = ("KK","JM2","JM1","JM2","JM1");
my @option = ("KK","JM","JM","JMM","JMM");

for (my $m = 0; $m <=$#method; $m++){
for (my $k = 0; $k <=$#option; $k++){

my $f1 = "</home/yhoang/Godwin/sift/new/T_".$option2[$k]."_bwa_$method[$m].genome_summary.csv";
my $f2 = "</home/yhoang/Godwin/sift/new/N_".$option[$k]."_bwa_$method[$m].genome_summary.csv";
my $cns = "</home/yhoang/Godwin/annotated/N_$option[$k]_cns13";

print "\nComparing in $method[$m] T_$option2[$k] with N_$option[$k]..\n";

my $write1 = ">/home/yhoang/Godwin/SNVs/con_bwa_$method[$m]\_T_$option2[$k]\_N_$option[$k].xls";
my $write2 = ">/home/yhoang/Godwin/SNVs/truedis_bwa_$method[$m]\_T_$option2[$k]\_N_$option[$k].xls";

my %tumor; my $it1 = 0;
open F1, $f1 or die "No $f1!";
while(<F1>){
	chomp();
	my @l = split(/,/,$_);			# array @l w/o tab
	if ($l[0] eq "Func") {next;}
	++$it1;
	my $gene = $l[1];
	my $pos = $l[21]."\t".$l[22]."\t".$l[24];
	my $var = $l[25];
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
		# quality score\tread depth\tRMS mapping quality
		$tumor{$pos}{$code} = $gene."\t".$l[0]."\t".$l[13]."\t".$l[7]."\t".$l[8]."\t".$l[27]."\t".$l[28]."\t".$l[29];  
       	} else { $tumor{$pos}{$var} = $gene."\t".$l[0]."\t".$l[13]."\t".$l[7]."\t".$l[8]."\t".$l[27]."\t".$l[28]."\t".$l[29]; }
}
close(F1);

my %normal; my $it2 = 0;
open F2,$f2 or die "No $f2!";
while(<F2>){
	chomp();
	my @l = split(/,/,$_);			# array @l w/o tab
	if ($l[0] eq "Func") {next;}
	++$it2;
	my $gene = $l[1];
	my $pos = $l[21]."\t".$l[22]."\t".$l[24];
	my $var = $l[25];
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
	#	print "$l[26]\t$l[25]\/$l[24] --> $code\n";
		# gene name\tfunction\tSIFT\tquality score\tread depth\tRMS mapping quality
		$normal{$pos}{$code} = $l[27]."\t".$l[28]."\t".$l[29];
       	} else { $normal{$pos}{$var} = $l[27]."\t".$l[28]."\t".$l[29]; }
}
close(F2);

my %con, my %dis,;
my $itcon = 0, my $itdis = 0;
foreach my $pos ( sort keys %tumor){  	# $pos: chr \t pos \t ref	
	if (defined $normal{$pos}){
		foreach my $var (sort keys %{$tumor{$pos}}){
			if (defined $normal{$pos}{$var}) {
			     #   print "$tumor{$pos}{$var}\t$var\t$normal{$pos}{$var}\n";
				$con{$pos}{$var} = $tumor{$pos}{$var}."\t".$var."\t".$normal{$pos}{$var};
				++$itcon;
			}
			elsif (!exists $dis{$pos}{$var} && !exists $con{$pos}{$var}) {
				$dis{$pos}{$var} = $tumor{$pos}{$var};
				++$itdis;
			}
		}
	}
	else{
		foreach my $var (sort keys %{$tumor{$pos}}){
			if  (!exists $dis{$pos}{$var} && !exists $con{$pos}{$var}) {
			$dis{$pos}{$var} = $tumor{$pos}{$var};
			++$itdis;
			}
		}
	}
}


my $sum1 = $itcon + $itdis;
print "Found $itdis discordances and $itcon concordances out of $sum1($option[$k])/$it2($option[$k]). Looking for true discordances...\n";


##### search for true discordance in cns.pileup
my $true = 0;
my %truedis;

my %file3;
print "opens $cns...\n";
open F3, $cns or die "No $cns!";
while(<F3>){
        chomp();
        my @l = split(/\t/,$_);                 # array @l w/o tab
        my $pos = $l[0]."\t".$l[1]."\t".$l[3];  # hashkey: chr \t pos \t ref
#        my $var = $l[3];
        if (exists $dis{$pos}){
#                if ( !defined $dis{$pos}{$var}){
                $truedis{$pos} = $dis{$pos};    # it should be in cns-pileup!!
                ++$true;
#                if ($true%1000 == 0) { print "Found $true true discordances out of $itdis...\n";}
                next;
#                }
        }
}
close(F3);

open WRITE1, $write1 or die;
if ($m == 0) { print WRITE1 "chr\tpos\tref\tvarT\tgene name\tfunction\tSIFT\t1000g2012feb_ALL\tdbSNP135\tquality score\tread depth\tRMS mapping quality\tvarB\tquality score\tread depth\tRMS mapping quality\n"; }
else { print WRITE1 "chr\tpos\tref\tvarT\tgene name\tfunction\tSIFT\t1000g2012feb_ALL\tdbSNP135\tsnp quality\ttotal reads\treads with mutation\tvarB\tsnp quality\ttotal reads\treads with mutation\n"; }
foreach my $pos ( sort keys %con ){
	foreach my $var (sort keys %{$con{$pos}} ){
		print WRITE1 join ("\t",$pos,$var,$con{$pos}{$var}."\n");
	}
}
close(WRITE1);

open WRITE2, $write2 or die;
print "$true true discordances out of $itdis discordances out of $sum1($option[$k]) possible SNVs in total.\n";
if ($m == 0) { print WRITE2 "chr\tpos\tref\tvarT\tgene name\tfunction\tSIFT\t1000g2012feb_ALL\tdbSNP135\tquality score\tread depth\tRMS mapping quality\n"; }
else { print WRITE2 "chr\tpos\tref\tvarT\tgene name\tfunction\tSIFT\t1000g2012feb_ALL\tdbSNP135\tsnp quality\ttotal reads\treads with mutation\n"; }
foreach my $pos ( sort keys %truedis ){
	foreach my $var (sort keys %{$truedis{$pos}} ){
		print WRITE2 join ("\t",$pos,$var,$truedis{$pos}{$var}."\n");
	}
}
close(WRITE2);

print "Written in $write1.\nWritten in $write2.\n";

}}
print strftime "\n%Y-%m-%d %H:%M:%S", localtime(time);
print "\nDone.\n";
