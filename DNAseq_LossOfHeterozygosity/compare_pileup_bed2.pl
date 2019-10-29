#!/usr/bin/perl -w

use strict;
use POSIX qw/strftime/;
print strftime "%Y-%m-%d %H:%M:%S", localtime(time);
print "\n";
print STDERR "./compare_pileup_bed <index> (2,4,5,6,7,12)\n";
my $index = $ARGV[0];

my $f1 = "</home/jchien/Godwin/$index/MMR_bwa/cns.rmd.pileup";		
my $f2 = "</project/data/Godwin/TruSeq-Exome-Targeted-Regions-BED-file";	#array
my $write = ">/home/jchien/Godwin/$index/targeted.pileup";
my $write2 = ">/home/jchien/Godwin/$index/untargeted.pileup";

my %pileup;
my $it_pileup = 0;	
open F1, $f1 or die "no $f1";
while(<F1>){
	chomp();
	my $match = 0;
	my @l = split(/\t/,$_);			# array @l w/o tab
	my $chr = $l[0];
	my $pos = $l[1];# hashkey column 0 und 1: chr \t pos \t ref
	$pileup{$chr}{$pos} = $l[2]."\t".$l[3]."\t".$l[7]."\t".$l[5]."\t".$l[6];
	++$it_pileup;

}
close(F1);
print "$it_pileup pileup positions.\n";

my @bed; my @sorted_bed;
my $it_bed = 0;
open F2, $f2 or die "no $f2";
while(<F2>){	#array data
	chomp();
	my @BED = split(/\t/,$_);		# array @l w/o tab
#	print "$BED[0]\t$BED[1]\t$BED[2]\n";
	my @chr = split(undef,$BED[0]);
#	print "$chr[0]\t$chr[1]\t$chr[2]\n";
	my $chr2;
	if ($#chr > 3) { $chr2 = $chr[3].$chr[4];}
	else { $chr2 = $chr[3];}
	my $start = $BED[1];
	my $end = $BED[2];
	$bed[$it_bed][0] = $chr2;
	$bed[$it_bed][1] = $start;
	$bed[$it_bed][2] = $end;
	++$it_bed;
}
close(F2);
print "$it_bed positions in BED-file.\n";
#@sorted_bed =
#     map { $_->[2] } # step 3
#     sort { $a->[0] <=> $b->[0] || $a->[1] cmp $b->[1] } # step 2
#     map { [ length $_[1], $_[0], $_ ] } # step 1
#     @users;

my %target, my %untarget;
my $con = 0; my $dis = 0; my $it = 0;
foreach my $chr ( sort keys %pileup ){
	++$it; 
	my $match = 0;
	for (my $i = 0; $i < $it_bed; $i++){
		if ($chr = $bed[$i][0]){
			foreach my $pos ( sort keys %{$pileup{$chr}} ){
				if ($pos >= $bed[$i][1] && $pos <= $bed[$i]) { 
				$target{$chr} = $pileup{$chr};	# info
				++$con;	
				$match = 1;
				}
			}
		}
	if ($match == 0) { $untarget{$chr} = $pileup{$chr}; }
	if ($it%10000 == 0) {print"run $it...\n";}
	}
}



print "writes $write\n";
open WRITE, $write or die;
print WRITE "chr\tpos\tref\tvar\tcov\tSNPqual\tmapqual\n";
foreach my $chr ( sort keys %target ){
	print WRITE $chr."\t".$target{$chr}."\n";
}
close(WRITE);

print "writes $write2\n";
open WRITE2, $write2 or die;
print WRITE2 "chr\tpos\tref\tvar\tcov\tSNPqual\tmapqual\n";
foreach my $chr ( sort keys %untarget ){
	print WRITE2 $chr."\t".$untarget{$chr}."\n";
}
close(WRITE2);

my $ratio = $con/$it;

print "Done. $con target positions out of $it --> ratio = $ratio.\n";
print strftime "%Y-%m-%d %H:%M:%S", localtime(time);
print "\n";
