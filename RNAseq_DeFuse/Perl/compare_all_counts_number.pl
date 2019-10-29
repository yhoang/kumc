#!/usr/bin/perl -w

#chmod +x compare_all_counts_number.pl
#./compare_all_counts_number.pl samples_number >dendro.txt

# NOT USEFUL YET!!!
# TRY TO CHANGE LIKE compare_myself_bothways.pl

use strict;
use warnings;

my @memory;
my $i     = 0;
my $run = 0;
my $open1 = $ARGV[0].".txt";
my @counts;
open LISTi, "<" . $open1 or die "no open1";
while (<LISTi>) {
	chomp();
	my $arrayi = $_;
	$memory[$i] = $arrayi;
	++$i;
	open LISTj, "<" . $open1 or die "no open1";
	while (<LISTj>) {
		my $it = 0;
		chomp();
		my $arrayj = $_;
		my $repeated = 0;
		foreach my $m (@memory) {
			if ( $arrayj eq $m ) { $repeated = 1; }
		}
		if ($repeated == 1) { next; }
		else {
			++$run;
			my %file1;
			my $f1 = $arrayi;
			open F1, "<sample_number/" . $arrayi or die "no f1";
			while (<F1>) {
				chomp();
				my @l = split( /\t/, $_ );     #array @l without tabs
				my $key1 = $l[30]."_".$l[31]."_".$l[32]."_".$l[33];                   #hashkey von Spalte 1 bis 3
				$file1{$key1} = join( "\t",$l[26],$l[27],$l[39],$l[40]);    #$file1{$key} von Spalte 11 bis 16
			}
			close(F1);

			my %file2;
			my $f2 = $arrayj;
			open F2, "<sample_number/" . $f2 or die "no f2";
			while (<F2>) {
				chomp();
				my @l = split( /\t/, $_ );    #array @l ohne ","
				my $key2 = $l[30]."_".$l[31]."_".$l[32]."_".$l[33];                   #hashkey von Spalte 1 bis 3
				$file2{$key2} = join( "\t",$l[26],$l[27],$l[39],$l[40]);    #$file1{$key} von Spalte 11 bis 16
			}
			close(F2);

			foreach my $key1 ( sort keys %file1 ) {
				if ( $key1 eq
					"gene_location1_gene_location2_gene_name1_gene_name2" )
				{
					next;
				}
				elsif ( $key1 eq "\"gene_location1\"_\"gene_location2\"_\"gene_name1\"_\"gene_name2\"" ){		
						next;
					}
				else {
					foreach my $key2 ( sort keys %file2 ) {
						if ( $key1 eq $key2 ) { ++$it; }
					}
				}
			}
			
			$counts[$run]=$it;
			print "$arrayi\t$arrayj\t$counts[$run]\n";
		}
		
	}
	close(LISTj);
}
close(LISTi);


