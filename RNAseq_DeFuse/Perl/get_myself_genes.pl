#!/usr/bin/perl -w

#chmod +x compare_everything.pl
#./compare_everything.pl samples_number >dendrogram.txt
#output for dendrogram in R

#chmod +x compare_everything_bothways.pl
#./compare_everything_bothways.pl samples_number >dendrogram.txt
#output: nxn table with relative match amounts for dendrogram in R
#plus ordering interchromosomal fusions (chr1 chr2)

#chmod +x compare_myself_bothways.pl 
#./compare_myself_bothways.pl samples_number >dendrogram.txt
#input now from samples already written out
#output: nxn table with relative match amounts for dendrogram in R
#no need of hashes
#and no need of interchromosomal comparison, because writeout_myself_bothways already did

#chmod +x compare_myself_bothways2.pl 
#./compare_myself_bothways2.pl samples_number >dendrogram.txt
#input now from samples already written out
#output: nxn table with relative match amounts
#and no need of hashes
#and no need of interchromosomal comparison, because writeout_myself_bothways already did
#and ordering intrachromosomal fusions (chr1 chr1)

use strict;
use warnings;

my $count = 1;
my $open1 = $ARGV[0].".txt";
my %file1;

my %genestart;	# from 'ENST_NAME'
my %genestart2;	# from 'geneID'
my $error = 0;

open geneName, "</home/yhoang/workspace_mRNAseq/ENST_Name" or die $!;
while (<geneName>) {
	chomp();
	my @Ensemble = split( /\t/, $_ );    #array @l without tabs
	my $genename =  $Ensemble[0];	# Associated Gene Name
	$genestart{$genename} = $Ensemble[1];	# Gene Start (bp)
}
close (geneName);

open geneName2, "</home/yhoang/workspace_mRNAseq/geneID" or die $!;
while (<geneName2>) {
	chomp();
	my @genelist = split( /\t/, $_ );    #array @l without tabs
	my $genename =  $genelist[0];	# ENSG Name
	$genestart2{$genename} = $genelist[1];	# Gene Start (bp)
}
close (geneName2);


open LISTi, "<" . $open1 or die "no open1";
while (<LISTi>) {
	chomp();
	my $arrayi = $_;
	
	#for label just with numbers:
	open WRITER, ">compare_myself/".$count or die $!;
	++$count;

	#or label with sample#.txt
	#open WRITER, ">compare_myself/".$arrayi.".txt" or die $!;

	open F1, "<samples/" . $arrayi or die "no f1";
	while (<F1>) {
		chomp();
		my @l = split( /\t/, $_ );    #array @l without "\t"
		
		# skip header
		if ($l[26]eq "gene_chromosome1" || $l[26] eq "\"gene_chromosome1\""){next;}
		if ($l[27]eq "gene_chromosome2" || $l[27] eq "\"gene_chromosome2\""){next;}
		if ($l[26]eq "Y" || $l[26] eq "\"Y\""){next;}
		if ($l[27]eq "Y" || $l[27] eq "\"Y\""){next;}		
	
		else{ 
			# change chr to int
			if ($l[26]eq "MT" || $l[26] eq "\"MT\""){$l[26]=0;}	
			if ($l[27]eq "MT" || $l[27] eq "\"MT\""){$l[27]=0;}
			if ($l[26]eq "X" || $l[26] eq "\"X\""){$l[26]=23;}
			if ($l[27]eq "X" || $l[27] eq "\"X\""){$l[27]=23;}

			# order all chr and gene start
			if ($l[26] > $l[27]) {	
				my $chr = $l[26];
				$l[26] = $l[27];
				$l[27] = $chr;
				my $gene = $l[32];
				$l[32] = $l[33];
				$l[33] = $gene;
			}
			elsif ($l[26] == $l[27]) {
				if ($l[34] > $l[35]) {
					my $gene2 = $l[32];
					$l[32] = $l[33];
					$l[33] = $gene2;
				}
			}

			# change ID to gene start position
			my $match1 = 0;
			my $match2 = 0;
			my $save_genename1 = $l[32];
			my $save_genename2 = $l[33];

			foreach my $genekey (keys %genestart){
				if ($l[32] eq "\"$genekey\"" || $l[32] eq $genekey){
				#	print "$l[32] changed to ";
					$l[32] = $genestart{$genekey};
					$match1 = 1;
				#	print "$l[32].\n";
				}
			}
			if ($match1 == 0) {
			#	print "error#$error in $arrayi [32]: $l[32]\n";
				foreach my $genekey (keys %genestart2){
					if ($l[22] eq $genekey){
					#	print "$l[32] changed to ";
						$l[32] = $genestart2{$genekey};
						$match1 = 1;
					#	print "$l[32].\n";
					}
				}
			}
			if ($match1 == 0) {
				++$error;
				print "still error#$error in $arrayi [32]: $l[32]\n";
			}

			foreach my $genekey (keys %genestart){
				if ($l[33] eq "\"$genekey\"" || $l[33] eq $genekey){
				#	print "$l[33] changed to ";
					$l[33] = $genestart{$genekey};
					$match2 = 1;
				#	print "$l[33].\n";
				}
			}
			if ($match2 == 0) {
		#		print "error#$error in $arrayi [33]: $l[33]\n";
				foreach my $genekey (keys %genestart2){
					if ($l[23] eq $genekey){
					#	print "$l[32] changed to ";
						$l[33] = $genestart2{$genekey};
						$match2 = 1;
					#	print "$l[32].\n";
					}
				}
			}
			if ($match2 == 0) {
				++$error;
				print "still error#$error in $arrayi [33]: $l[33]\n";
			}		
			print WRITER join ("\t",$l[26],$l[27],$l[32],$l[33],$save_genename1,$save_genename2,$l[1]."\n");
		}
	}
	close(F1);
	close(WRITER);
}
close(LISTi);

