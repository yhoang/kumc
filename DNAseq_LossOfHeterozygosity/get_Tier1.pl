#!/usr/bin/perl -w
use strict;
use POSIX qw/strftime/;

print strftime "%Y-%m-%d %H:%M:%S\n", localtime(time);
print STDERR "./get_Tier1.pl\n--> Tier1: Compare all 3 variant calling tools (samtools,gatk,sniper) and get the concordant calls.\n";

my @method = ("gatk","sam","sniper");
my @option2 = ("KK","JM2","JM1","JM2","JM1");   #tumor T_
my @option = ("KK","JM","JM","JMM","JMM");      #normal N_

for (my $k = 0; $k <=$#option; $k++){
my @file;
$file[0] = "</home/yhoang/Godwin/SNVs/truedis_bwa_$method[0]\_T_$option2[$k]\_N_$option[$k].xls";
$file[1] = "</home/yhoang/Godwin/SNVs/truedis_bwa_$method[1]\_T_$option2[$k]\_N_$option[$k].xls";
$file[2] = "</home/yhoang/Godwin/SNVs/truedis_bwa_$method[2]\_T_$option2[$k]\_N_$option[$k].xls";

my $write = ">/home/yhoang/Godwin/SNVs/concordance/T1_bwa_T_$option2[$k]\_N_$option[$k].xls";

my $title;

my %gatk;
my $it1 = 0;
print "opens $file[0]...\n";
open F1, $file[0] or die "no $file[0]";
while(<F1>){
	chomp();
	my @l = split(/\t/,$_);			# array @l w/o tab
	if ($l[0] eq "chr") { $title = $l[0]."\t".$l[1]."\t".$l[2]."\tvar_gatk\t".$l[4]."\t".$l[5]."\t".$l[6]."\t".$l[7]."\t".$l[8]."\t".$l[9]."\t".$l[10]."\t".$l[11]; next;}
	my $pos = $l[0]."\t".$l[1]."\t".$l[2];	# hashkey: chr \t pos \t ref
	my $var = $l[3];
	$gatk{$pos}{$var} = $l[4]."\t".$l[5]."\t".$l[6]."\t".$l[7]."\t".$l[8]."\t".$l[9]."\t".$l[10]."\t".$l[11]; 
	++$it1;
}
close(F1);

my %samtools;
my $it2 = 0;
print "opens $file[1]...\n";
open F2, $file[1] or die "no $file[1]";
while(<F2>){
	chomp();
	my @l = split(/\t/,$_);			# array @l w/o tab
	if ($l[0] eq "chr") {$title = $title."\tvar_samtools\t".$l[9]."\t".$l[10]."\t".$l[11]; next;}
	my $pos = $l[0]."\t".$l[1]."\t".$l[2];	# hashkey: chr \t pos \t ref
	my $var = $l[3];
	$samtools{$pos}{$var} =  $var."\t".$l[9]."\t".$l[10]."\t".$l[11];  
	++$it2;
}
close(F2);


my %T2;
my $itT2 = 0;
foreach my $pos ( sort keys %gatk){  	# $pos: chr \t pos \t ref	
	if (defined $samtools{$pos}){
		foreach my $var (sort keys %{$gatk{$pos}}){
			if (defined $samtools{$pos}{$var}) {
				$T2{$pos}{$var} = $gatk{$pos}{$var}."\t".$samtools{$pos}{$var};
				++$itT2;
			}
		}
	}
}
print "Found $itT2 concordances out of $it1($method[0])/$it2($method[1]).\n";

my $it3 = 0;
my %sniper;
print "opens $file[2]...\n";
open F3, $file[2] or die "no $file[2]";
while(<F3>){
	chomp();
	my @l = split(/\t/,$_);			# array @l w/o tab
	my $pos = $l[0]."\t".$l[1]."\t".$l[2];	# hashkey: chr \t pos \t ref
	my $var = $l[3];
	$sniper{$pos}{$var} = $var."\t".$l[8]."\t".$l[12]; 
	++$it3;
}
close(F3);
$title = $title."\tvar_sniper\tmapping\tcoverage"; 

my %T1;
my $itT1 = 0;
foreach my $pos ( sort keys %T2){  	# $pos: chr \t pos \t ref	
	if (defined $sniper{$pos}){
		foreach my $var (sort keys %{$T2{$pos}}){
			if (defined $sniper{$pos}{$var}) {
				$T1{$pos}{$var} = $T2{$pos}{$var}."\t".$sniper{$pos}{$var};
				++$itT1;
			}
		}
	}
}

print "In the end: Found $itT1 T1 concordances (after $itT2) out of $it1($method[0])/$it2($method[1])/$it3($method[2]).\n";

open WRITE, $write or die;
print WRITE "$title\n";
foreach my $pos ( sort keys %T1 ){
	foreach my $var (sort keys %{$T1{$pos}} ){
		print WRITE $pos."\t".$var."\t".$T1{$pos}{$var}."\n";
	}
}
close(WRITE);
print "Written in $write.\n\n";
}
print strftime "%Y-%m-%d %H:%M:%S\n", localtime(time);
print "Done.\n";
