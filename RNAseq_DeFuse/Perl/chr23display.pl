#!/usr/bin/perl -w

#chmod +x get_filtered.pl
#./get_filtered.pl samplelist f
#plus ordering interchromosomal fusions (chr1 chr2)
#plus ordering intrachromosomal fusions (chr1 chr1)
#output >filteredf/sample1..36.txt  writeout chr1 chr2 genestart1 genestart2


use strict;
use warnings;

my $open1 = "</home/yhoang/workspace_Mariani/DeFuse_Perl/".$ARGV[0].".txt";
my $option = "</home/yhoang/workspace_Mariani/DeFuse_Perl/samples/";
my $writer = ">/home/yhoang/workspace_Mariani/DeFuse_Perl/chr2_chr3_display.txt";
my %file1;
open WRITER, $writer or die $!;
open LISTi, $open1 or die "no open1";
while (<LISTi>) {
	chomp();
	my $arrayi = $_;
	
	
	open F1, $option . $arrayi or die "no f1";
	while (<F1>) {
		chomp();
		my @l = split( /\t/, $_ );    #array @l without tabs

		if ($l[26]eq "gene_chromosome1" || $l[26] eq "\"gene_chromosome1\""){next;}	# skip header
		if ($l[27]eq "gene_chromosome2" || $l[27] eq "\"gene_chromosome2\""){next;}
		if ($l[26]eq "Y" || $l[26] eq "\"Y\""){next;}
		if ($l[27]eq "Y" || $l[27] eq "\"Y\""){next;}		
		else{ 
			if ($l[26]eq "MT" || $l[26] eq "\"MT\""){$l[26]=0;}	# change to int
			if ($l[27]eq "MT" || $l[27] eq "\"MT\""){$l[27]=0;}
			if ($l[26]eq "X" || $l[26] eq "\"X\""){$l[26]=23;}
			if ($l[27]eq "X" || $l[27] eq "\"X\""){$l[27]=23;}
	
			if ($l[26] > $l[27]){	# order all chr and gene start
				my $chr = $l[26];
				$l[26] = $l[27];
				$l[27] = $chr;
				my $gene = $l[34];
				$l[34] = $l[35];
				$l[35] = $gene;
			}
			elsif ($l[26] == $l[27]){
				if ($l[34] > $l[35]){
			#		print "hier in $arrayi $l[26]\t$l[27]\t$l[34]\t$l[35]\t to\t";
					my $gene2 = $l[34];
					$l[34] = $l[35];
					$l[35] = $gene2;
			#		print "$l[26]\t$l[27]\t$l[34]\t$l[35].\n";
				}
			}
			if (($l[26] == 2) && ($l[27] == 3)){
			my $key1 = $l[30]."_".$l[31]."_".$l[32]."_".$l[33];    # hashkey loc1 loc2 name1 name2
			$file1{$key1} = join( "\t",$l[26],$l[27],$l[34],$l[35]); # hash chr1 chr2 start1 start2
			
			print WRITER $arrayi.": ".$file1{$key1}."\n";
			print $arrayi.": ".$file1{$key1}."\n";
			}
		}
	}
	close(F1);
	
}
close(LISTi);
close(WRITER);
