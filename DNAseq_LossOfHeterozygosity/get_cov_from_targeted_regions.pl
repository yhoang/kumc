#!/usr/bin/perl -w
use strict;
use POSIX qw/strftime/;
print strftime "%Y-%m-%d %H:%M:%S", localtime(time);
print "\n";
print STDERR "./get_cov_from_targeted_regions.pl <sample> (T_KK,N_KK,T_JM2,T_JM2,N_JM,N_JMM)\nGet mean coverage from single targeted regions and from single genes.\n\n";
my $sample = $ARGV[0];

my @method = ("bwa");
#my @option2 = ("FF1","FFPE1","FF2","FFPE2","FF3","FFPE3");

#foreach my $sample (@option2) {
foreach my $k (@method) {

my $f1 = my $f1 = "</home/yhoang/Godwin/annotated/changed_code/$sample\_$k\_cns13_var";	
my $f2 = "</project/data/Godwin/targeted_regions.bed";	#array
my $write = ">/home/yhoang/Godwin/targeted_regions/$sample\_$k\_cov13_targetedfiltered.xls";
my $write2 = ">/home/yhoang/Godwin/targeted_regions/$sample\_$k\_cov13_genefiltered.xls";

my %pileup;
my $it_pileup = 0;	
print "opens $f1...\n";
open F1, $f1 or die "no $f1";
while(<F1>){
	chomp();
	my @l = split(/\t/,$_);			# array @l w/o tab
	if ($l[0] eq "chr") {next;}
	my $chr = $l[0];
	my $pos = $l[1];# hashkey column 0 und 1: chr \t pos \t ref
	$pileup{$chr}{$pos} = $l[8];    # coverage
	++$it_pileup;
}
close(F1);
print "$it_pileup consensus pileup positions.\n";

my %bed; my %genename; my %position;
my $it_bed = 0;
print "opens $f2...\n";
open F2, $f2 or die "no $f2";
while(<F2>){	#array data
	chomp();
	my @b = split(/\t/,$_);		# array @l w/o tab
	my @info = split(/:/,$b[3]); #print "$info[0]\t$info[1]\t$info[2]\n";
	my $genename = $info[2];
	my $chr = $b[0];
	my $start = $b[1];
	my $end = $b[2];
	$bed{$chr}{$start} = $end;
	$genename{$chr}{$start} = $genename;
	$position{$chr}{$start} = $b[2]."\t".$b[3]."\t".$b[4];
	++$it_bed;
}
close(F2);
print "$it_bed bed-file positions.\n";

my $tmp = 0; my @count_gene; my $it_gene = 0; my $pos_gene = 0;
my @count_targeted; my $it_start = 0; my $it_count1 = 0;
my $it = 0;
$count_gene[0][0] = 0;

open WRITE, $write or die;
print WRITE "regionIT\tcov\tposition count\taverage\tgene name\tchr\tstart\tend\n";

foreach my $chr ( sort keys %bed ){
	if (defined $pileup{$chr}) {
		foreach my $start ( sort keys %{$bed{$chr}} ){
			if ( $genename{$chr}{$start} ne $tmp && $tmp ne "0") {++$it_gene; $pos_gene = 0; $count_gene[$it_gene][0] = 0;} #print "\tit_gene = $it_gene goes up because $tmp ne $genename{$chr}{$start}!\n";
                        $count_gene[$it_gene][2] = $genename{$chr}{$start};
			++$it_count1;# = $it + $it_start;	# targeted position count
			$count_targeted[$it_count1][0] = 0;
			my $pos_targeted = 0;
			foreach my $pos (sort keys %{$pileup{$chr}}){
				if ($pos >= $start && $pos <= $bed{$chr}{$start}) {
					++$pos_targeted;
					$count_targeted[$it_count1][0] = $count_targeted[$it_count1][0] + $pileup{$chr}{$pos};
					$count_gene[$it_gene][0] = $count_gene[$it_gene][0] + $pileup{$chr}{$pos};
					++$pos_gene;	
				}	
			}
			$count_targeted[$it_count1][1] = $pos_targeted;
			my $average = 0;
			if ($pos_targeted != 0) {$average = $count_targeted[$it_count1][0] / $pos_targeted;}
			print WRITE "$it_count1\t$count_targeted[$it_count1][0]\t$count_targeted[$it_count1][1]\t$average\t$genename{$chr}{$start}\t$position{$chr}{$start}\n";
			++$it_start;
			$tmp = $genename{$chr}{$start};
			$count_gene[$it_gene][1] = $pos_gene;
		}
		my $average2 = $count_gene[$it_gene][0] / $pos_gene;
		print "$pos_gene positions in $it_gene genes in chr $chr --> cov = $count_gene[$it_gene][0], average = $average2\n";
	}
	++$it;
}
close(WRITE);
print "$it_gene genes out of $it_count1 targeted positions (should be $it_bed).\n";

open WRITE2, $write2 or die;
print WRITE2 "gene name\tcov\tposition count\taverage\n";
for (my $i = 0; $i <= $it_gene; $i++) {
	my $average2 = $count_gene[$i][0] / $count_gene[$i][1];
	print WRITE2 "$count_gene[$i][2]\t$count_gene[$i][0]\t$count_gene[$i][1]\t$average2\n";
}
close(WRITE2);
print "Written in $write and $write2.\n";

}
print "Done. \n";
print strftime "%Y-%m-%d %H:%M:%S", localtime(time);
print "\n";
