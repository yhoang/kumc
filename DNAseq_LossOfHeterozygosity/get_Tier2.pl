#!/usr/bin/perl -w
use strict;
use POSIX qw/strftime/;

print strftime "%Y-%m-%d %H:%M:%S\n", localtime(time);
print STDERR "./get_Tier2.pl\n--> Tier2: Get somatic variants from any two calling tools (samtools,gatk,sniper).\n";

my @method = ("gatk","sam","sniper");
my @option2 = ("KK","JM2","JM1","JM2","JM1");
my @option = ("KK","JM","JM","JMM","JMM");

for (my $k = 0; $k <=$#option; $k++){
print "Get Tier2 somatic mutations from T_$option2[$k]_N_$option[$k]...\n";
my %T2;
my $T2 = 0;
for (my $m = 0; $m <=$#method; $m++){
for (my $n = 0; $n <=$#method; $n++){
if ($m == $n || $m > $n) {next;}

my $f1 = "</home/yhoang/Godwin/somatic_mutations/truedis_bwa_$method[$m]\_T_$option2[$k]\_N_$option[$k].xls";
my $f2 = "</home/yhoang/Godwin/somatic_mutations/truedis_bwa_$method[$n]\_T_$option2[$k]\_N_$option[$k].xls";

my %file1;
my $it1 = 0;
print "opens $f1...\n";
open F1, $f1 or die "no $f1";
while(<F1>){
	chomp();
	my @l = split(/\t/,$_);			# array @l w/o tab
	if ($l[0] eq "chr") {next;}
	my $pos = $l[0]."\t".$l[1]."\t".$l[2];	# hashkey: chr \t pos \t ref
	my $var = $l[3];
	$file1{$pos}{$var} = $l[4]."\t".$l[5]."\t".$l[6]."\t".$l[7]."\t".$l[8]."\t".$l[9]."\t".$l[10]."\t".$l[11]; 
	++$it1;
}
close(F1);

my $it2 = 0;
my %file2;
print "opens $f2...\n";
open F2, $f2 or die "no $f2";
while(<F2>){
	chomp();
	my @l = split(/\t/,$_);			# array @l w/o tab
	if ($l[0] eq "chr") {next;}
	my $pos = $l[0]."\t".$l[1]."\t".$l[2];	# hashkey: chr \t pos \t ref
	my $var = $l[3];
	$file2{$pos}{$var} = $var."\t".$l[9]."\t".$l[10]."\t".$l[11]; 
	++$it2;
}
close(F2);

my %con;
my $itcon = 0, my $itdis = 0;
foreach my $pos ( sort keys %file1){  	# $pos: chr \t pos \t ref	
	if (defined $file2{$pos}){
		foreach my $var (sort keys %{$file1{$pos}}){
			if (defined $file2{$pos}{$var}) {
				$con{$pos}{$var} = $file1{$pos}{$var}."\t".$file2{$pos}{$var};
				++$itcon;
			}
		}
	}
}

print "Found $itcon concordances out of $it1($method[$m])/$it2($method[$n]).\n";
my $repeat = 0;
foreach my $pos ( sort keys %con){  	# $pos: chr \t pos \t ref	
	if (defined $T2{$pos}){ ++$repeat;}
	else {
		$T2{$pos} = $con{$pos};
		++$T2;
	}
}
my $added  = $itcon - $repeat;
print "$added positions were added in T2 list ($repeat already exists).\n\n";



}}
my $write = ">/home/yhoang/Godwin/SNVs/concordance/T2_bwa_T_$option2[$k]\_N_$option[$k].xls";
open WRITE, $write or die;
print WRITE "chr\tpos\tref\tvar\tvar type\tfunction\tSIFT\t1000g2012feb_ALL\tdbSNP135\tquality score\tread depth\tRMS mapping quality\tvar\tquality\ttotal reads\treads with mutation\n";
foreach my $pos ( sort keys %T2 ){
	foreach my $var (sort keys %{$T2{$pos}} ){
		print WRITE $pos."\t".$var."\t".$T2{$pos}{$var}."\n";
	}
}
close(WRITE);
print "$T2 T2 somatic mutations found! Written in $write.\n\n";
}
print strftime "%Y-%m-%d %H:%M:%S", localtime(time);
print "Done.\n";
