#!/usr/bin/perl -w

#chmod +x compare_everything.pl
#./compare_everything.pl samples_number >dendrogram.txt
#output for dendrogram in R

#chmod +x compare_everything_bothways.pl
#./compare_everything_bothways.pl samples_number >dendrogram.txt
#output: nxn table with relative match amounts for dendrogram in R
#plus ordering interchromosomal fusions (chr1 chr2)

#chmod +x compare_myself_bothways.pl samples_number >dendrogram.txt
#input now from samples already written out
#output: nxn table with relative match amounts for dendrogram in R
#no need of hashes
#and no need of interchromosomal comparison, because writeout_myself_bothways already did

use strict;
use warnings;

my $i = 0;
my $noti = 0;
my $j = 0;
my $notj = 0;
my $open1 = $ARGV[0].".txt";
my @counts; 	# array nxn table with relative matches

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

		my @key1;
		my $c1 = 1;
		my $f1 = $arrayi;
		open F1, "<compare_myself/" . $arrayi or die "no f1";
		while (<F1>) {
			chomp();
			my @l = split( /\t/, $_ );    #array @l without tabs
			$key1[$c1] = $l[0]."_".$l[1]."_".$l[2]."_".$l[3];        
			++$c1;
		}
		close(F1);

		my @key2;
		my $c2 = 1;
		my $f2 = $arrayj;
		open F2, "<compare_myself/" . $f2 or die "no f2";
		while (<F2>) {
			chomp();
			my @l = split( /\t/, $_ );    #array @l without tabs
			$key2[$c2] = $l[0]."_".$l[1]."_".$l[2]."_".$l[3];        
			++$c2;
		}
		close(F2);

		for (my $d1 = 1;$d1 < $c1;$d1++){
			for (my $d2 = 1;$d2 < $c2;$d2++){
				if ($key1[$d1] eq $key2[$d2]){++$it; }
			}
		}
		$counts[$i][$j] = $it;
	}	
	close(LISTj);
	$j = 0;
}

my $n = 0;
my $m = 0;
my $percent;
for ($n = 1;$n <= 36;$n++){
	for ($m = 1;$m <= 36;$m++){
		$percent  = $counts[$n][$m] / $counts[$n][$n];
		print "$percent\t";
	}
	print"\n";
}

close(LISTi);


