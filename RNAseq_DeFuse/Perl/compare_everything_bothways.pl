#!/usr/bin/perl -w

#chmod +x compare_everything.pl
#./compare_everything.pl samples_number >dendrogram.txt
#output for dendrogram in R

#chmod +x compare_everything_bothways.pl
#./compare_everything_bothways.pl samples_number >dendrogram.txt
#output: nxn table with relative match amounts for dendrogram in R
#plus ordering interchromosomal fusions (chr1 chr2)

use strict;
use warnings;

my $i = 0;
my $noti = 0;
my $j = 0;
my $notj = 0;
my $open1 = $ARGV[0].".txt";
my @counts;	# array nxn table with relative matches

open LISTi, "<" . $open1 or die "no open1";
while (<LISTi>) {
	chomp();
	my $arrayi = $_;
	++$i;

	# for a nxn array
	# just for use if data is not complety numbered 
#	while ($i < $arrayi){
#		for ($noti=1;$noti<=36;$noti++){	
#			$counts[$i][$noti]=0;	
#		}	
#		++$i;
#	}

	open LISTj, "<" . $open1 or die "no open1";
	while (<LISTj>) {
		my $it = 0;
		chomp();
		my $arrayj = $_;
		++$j;

		# for a nxn array
		# just for use if data is not complety numbered 
#		while ($j < $arrayj){
#			$counts[$i][$j]=0;
#			++$j;
#		}

		my %file1;
		my $f1 = $arrayi;
		open F1, "<sample_number/" . $arrayi or die "no f1";
		while (<F1>) {
			chomp();
			my @l = split( /\t/, $_ );    #array @l without tabs
			if ($l[26] eq "gene_chromosome1" || $l[26] eq "\"gene_chromosome1\""){next;}	# skip header
			if ($l[27] eq "gene_chromosome2" || $l[27] eq "\"gene_chromosome2\""){next;}
			else{
				if ($l[26] eq "MT" || $l[26] eq "\"MT\"") { $l[26] = 0; }	# change to int
				if ($l[27] eq "MT" || $l[27] eq "\"MT\"") { $l[27] = 0; }
				if ($l[26] eq "X" || $l[26] eq "\"X\"") { $l[26] = 24; }
				if ($l[27] eq "X" || $l[27] eq "\"X\"") { $l[27] = 24; }

				if ($l[26] > $l[27]){	# order all chr and gene start
					my $chr = $l[26];
					$l[26] = $l[27];
					$l[27] = $chr;
					my $gene = $l[28];
					$l[28] = $l[29];
					$l[28] = $gene;
				}

				my $key = $l[30]."_".$l[31]."_".$l[32]."_".$l[33];   # hashkey loc1 loc2 name1 name2
				$file1{$key} = join( "\t",$l[26],$l[27],$l[28],$l[29]);     # hash chr1 chr2 start1 start2
			}
		}
		close(F1);

		my %file2;
		my $f2 = $arrayj;
		open F2, "<sample_number/" . $f2 or die "no f2";
		while (<F2>) {
			chomp();
			my @l = split( /\t/, $_ );    #array @l without tabs
			if ($l[26] eq "gene_chromosome1" || $l[26] eq "\"gene_chromosome1\""){next;}	# skip header
			if ($l[27] eq "gene_chromosome2" || $l[27] eq "\"gene_chromosome2\""){next;}
			else{
				if ($l[26] eq "MT" || $l[26] eq "\"MT\"") { $l[26] = 0; }	# change to int
				if ($l[27] eq "MT" || $l[27] eq "\"MT\"") { $l[27] = 0; }
				if ($l[26] eq "X" || $l[26] eq "\"X\"") { $l[26] = 24; }
				if ($l[27] eq "X" || $l[27] eq "\"X\"") { $l[27] = 24; }

				if ($l[26] > $l[27]){	# order all chr and gene start
					my $chr = $l[26];
					$l[26] = $l[27];
					$l[27] = $chr;
					my $gene = $l[28];
					$l[28] = $l[29];
					$l[28] = $gene;
				}

				my $key = $l[30]."_".$l[31]."_".$l[32]."_".$l[33];   # hashkey loc1 loc2 name1 name2
				$file2{$key} = join( "\t",$l[26],$l[27],$l[28],$l[29]);     # chr1 chr2 start1 start2
			}
		}
		close(F2);
		
		# count matches
		foreach my $key1 ( sort keys %file1 ) {
			foreach my $key2 ( sort keys %file2 ) {
				if ( $key1 eq $key2 ) {
					++$it;
				}
			}
		}
		$counts[$i][$j] = $it;
	}	
	close(LISTj);
	$j = 0;
}

# printout
my $n = 0;
my $m = 0;
my $percent;
for ($n = 1;$n <= 36;$n++){
	for ($m = 1;$m <= 36;$m++){
		$percent = $counts[$n][$m] / $counts[$n][$n];
		print "$percent\t";
	}
	print"\n";
}

close(LISTi);


