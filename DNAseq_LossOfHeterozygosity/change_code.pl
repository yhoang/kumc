#!/usr/bin/perl -w
use strict;
use POSIX qw/strftime/;
print strftime "%Y-%m-%d %H:%M:%S", localtime(time);
print "\n";
print STDERR "./change_code.pl <input_method> (bwa) <coverage> <ex> (optional) \n--> Change code to be comparable with normal cns.pileup... \n\n";

my $input = $ARGV[0];
my $coverage = $ARGV[1];
my $opt = $ARGV[2];

if ($opt eq "ex") { &anno_ex($input,$coverage);} 
elsif ($input eq "bwa" || $input eq "bfast" || $input eq "bfast3" || $input eq "clc3" || $input eq "sniper_clc" || $input eq "sniper_bwa") {&anno($input,$coverage);}
elsif ($input eq "clc" || $input eq "clc2") { &clc($input,$coverage);}


sub anno {
my $method = $_[0];
my $coverage = $_[1];
my @names;
@names = ("N_KK","T_KK","N_JM","N_JMM","T_JM1","T_JM2");

my $frequency = 10; 

my $f1= "hello", my $write;
foreach my $n (@names) {
if ($method eq "sniper_clc" || $method eq "sniper_bwa") {
	$f1 = "</home/yhoang/Godwin/sniper/$n\_$method\_var$coverage";
	$write = ">/home/yhoang/Godwin/annotated/changed_code/$n\_$method\_var$coverage";
} else {
        $f1 = "</home/yhoang/Godwin/annotated/cns13/$n\_cns$coverage";
	#$f1 = "</home/yhoang/Godwin/annotated/cns13/$n\_$method\_cns$coverage";
	$write = ">/home/yhoang/Godwin/annotated/changed_code/$n\_$method\_cns$coverage";
}

my $m = 0, my $code;
my %filtered;
print "opens $f1...\n";
open F1, $f1 or die "no $f1!\n";
while(<F1>){
	chomp();
	my @l = split(/\t/,$_);			# array @l w/o tab
	if ($l[10] < 2) {next;}
	if ($l[6] eq "-") {next;}
	my @insertion = split(undef,$l[6]);
	if ($#insertion > 1) {next;}
	my $pos = $l[2]."\t".$l[3]."\t".$l[5];	# hashkey column 0,1,5: chr \t pos \t ref
	my $var = $l[6];
	if ( $l[7] eq "het") {	# heterogenizity status
		# change code
		my @var_sorted = sort($l[5],$l[6]) ;
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
		$filtered{$pos}{$code} = $l[7]."\t".$l[0]."\t".$l[1]."\t". $l[8]."\t".$l[9]."\t".$l[10];
	}
	else {
		$filtered{$pos}{$var} = $l[7]."\t".$l[0]."\t".$l[1]."\t". $l[8]."\t".$l[9]."\t".$l[10];
	}
}
close(F1);

my $it = 0;
open WRITE, $write or die;
print WRITE "chr\tpos\tref\tvar\tstatus\tloc\tgene name\tsnp qual\tcoverage\tmutant counts\n";
foreach my 
$pos ( sort keys %filtered ){
	foreach my $var (sort keys %{$filtered{$pos}} ){
		print WRITE $pos."\t".$var."\t".$filtered{$pos}{$var}."\n";
		++$it;
	}
}
close(WRITE);
print "Written in $write.\n";

print "Done. $it positions.\n";
}
return 0;

}


sub anno_ex {
my $method = $_[0];
my $coverage = $_[1];
my $frequency = 10; #frequency

my @names;
if ($method eq "sniper_clc" || $method eq "sniper_bwa") { @names = ("CT_FF","CT_FFPE");}
else { @names = ("CN_FF","CN_FFPE","CT_FF","CT_FFPE");}

my $f1, my $write;
foreach my $n (@names) {
if ($method eq "sniper_clc" || $method eq "sniper_bwa") {
	$f1 = "</home/yhoang/Godwin/sniper/$n\_$method\_var$coverage\_ex";
	$write = ">/home/yhoang/Godwin/annotated/$n\_$method\_varf$coverage\_ex";
} else {
	$f1 = "</home/yhoang/Godwin/annotated/original/$n\_$method\_varf$coverage\_ex";
	$write = ">/home/yhoang/Godwin/annotated/$n\_$method\_varf$coverage\_ex";
}

my $m = 0, my $code;
my %filtered;

print "Opens $f1...\n";
open F1, $f1 or die "no $f1!\n";
while(<F1>){
	chomp();
	my @l = split(/\t/,$_);			# array @l w/o tab
	if ($l[11] < 2) {next;}
	my $pos = $l[3]."\t".$l[4]."\t".$l[6];	# hashkey column 0,1,5: chr \t pos \t ref
	my $var = $l[7];
	if ( $l[8] eq "het") {	# heterogenizity status
		# change code
		my @var_sorted = sort($l[6],$l[7]) ;
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
		$filtered{$pos}{$code} = $l[8]."\t".$l[1]."\t".$l[2]."\t".$l[9]."\t".$l[10]."\t".$l[11];
	}
	else {
		$filtered{$pos}{$var} = $l[8]."\t".$l[1]."\t".$l[2]."\t".$l[9]."\t".$l[10]."\t".$l[11];
	}
}
close(F1);

my $it = 0;
open WRITE, $write or die;
foreach my 
$pos ( sort keys %filtered ){
	foreach my $var (sort keys %{$filtered{$pos}} ){
		print WRITE $pos."\t".$var."\t".$filtered{$pos}{$var}."\n";
		++$it;
	}
}
close(WRITE);
print "Written $it SNVs in $write.\n\n";


}
print "Done.\n";
return 0;
}


sub clc {
my @names = ("N_KK","T_KK");
my $method = $_[0];
my $coverage = $_[1];
my $frequency = 10; 

foreach my $n (@names) {
my $f1 = "</home/yhoang/Godwin/clc/SNV/$n\_$method.txt";
my $write = ">/home/yhoang/Godwin/clc/SNV_changed_code/$n\_$method.txt";
my $write2 = ">/home/yhoang/Godwin/annotated/$n\_$method\_varf$coverage";

my $m = 0, my $code;
my %filtered;

open F1, $f1 or die "no $f1!\n./change_code.pl <input_method> (clc,clc2)\n";
while(<F1>){
	chomp();
	my @l = split(/\t/,$_);			# array @l w/o tab
	if ($l[2] ne "SNV") {next;}
	my $pos = $l[0]."\t".$l[1]."\t".$l[5];	# hashkey column 0,1,5: chr \t pos \t ref
	my @var = split (/\//,$l[6]);
	if ( $#var >= 1) {
		# change code
		my @var_sorted = sort(@var) ;
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
		my @count = split (/\//,$l[7]);
		my @frequency = split (/\//,$l[9]); 
		if ( ($count[0] + $count[1]) < $coverage) {next;}	#coverage
		if ( $count[0] < 2 || $count[1] < 2) {next;}
		if ( $frequency[0] < $frequency && $frequency[1] < $frequency) {next;}
		$filtered{$pos}{$code} = $l[7]."\t".$l[8]."\t".$l[9]."\t".$l[10]."\t".$l[11]."\t".$l[12];
	}
	elsif ($l[7] >= $coverage && $l[9] >= $frequency) {
		my $pos_var = $l[6];			# variant
		$filtered{$pos}{$pos_var} = $l[7]."\t".$l[8]."\t".$l[9]."\t".$l[10]."\t".$l[11]."\t".$l[12];	# info
	}
}
close(F1);

my $it = 0;
open WRITE, $write or die;
foreach my $pos ( sort keys %filtered ){
	foreach my $var (sort keys %{$filtered{$pos}} ){
		print WRITE $pos."\t".$var."\t".$filtered{$pos}{$var}."\n";
		++$it;
	}
}
close(WRITE);
print "Written in $write.\n";

open WRITE2, $write2 or die;
foreach my $pos ( sort keys %filtered ){
	foreach my $var (sort keys %{$filtered{$pos}} ){
		print WRITE2 $pos."\t".$var."\t".$filtered{$pos}{$var}."\n";
		++$it;
	}
}
close(WRITE2);
print "Written $it SNVs in $write2.\n";
}
print " Done.\n";
return 0;
}
