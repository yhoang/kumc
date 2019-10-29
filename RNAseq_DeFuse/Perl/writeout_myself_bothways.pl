#!/usr/bin/perl -w

#chmod +x writeout_myselfpl
#./writeout_myself.pl samples
#output >compare_myself/samplex.txt  writeout chr1 chr2 genestart1 genestart2

#chmod +x ./writeout_myself_bothways.pl
#./writeout_myself_bothways.pl samples
#plus ordering interchromosomal fusions (chr1 chr2)
#output: samples/sample1...36.txt  writeout chr1 chr2 genestart1 genestart2

use strict;
use warnings;

my $open1 = $ARGV[0].".txt";
my %file1;

open LISTi, "<" . $open1 or die "no open1";
while (<LISTi>) {
	chomp();
	my $arrayi = $_;
	
	open WRITER, ">compare_myself/".$arrayi.".txt" or die $!;

	open F1, "<samples/" . $arrayi or die "no f1";
	while (<F1>) {
		chomp();
		my @l = split( /\t/, $_ );    #array @l ohne ","

		if ($l[26]eq "gene_chromosome1" || $l[26] eq "\"gene_chromosome1\""){next;}	# skip header
		if ($l[27]eq "gene_chromosome2" || $l[27] eq "\"gene_chromosome2\""){next;}

		else{ 
			if ($l[26]eq "MT" || $l[26] eq "\"MT\""){$l[26]=0;}	# change to int
			if ($l[27]eq "MT" || $l[27] eq "\"MT\""){$l[27]=0;}
			if ($l[26]eq "X" || $l[26] eq "\"X\""){$l[26]=23;}
			if ($l[27]eq "X" || $l[27] eq "\"X\""){$l[27]=23;}
	
			if ($l[26] > $l[27]){	# order all chr and gene start
				my $chr = $l[26];
				$l[26] = $l[27];
				$l[27] = $chr;
				my $gene = $l[28];
				$l[28] = $l[29];
				$l[28] = $gene;
			}
			my $key1 = $l[30]."_".$l[31]."_".$l[32]."_".$l[33];    # hashkey loc1 loc2 name1 name2
			$file1{$key1} = join( "\t",$l[26],$l[27],$l[28],$l[29]); # hash chr1 chr2 start1 start2
			
			print WRITER $file1{$key1}."\n";
		}
	}
	close(F1);
	close(WRITER);
}
close(LISTi);

