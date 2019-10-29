#!/usr/bin/perl -w
use strict;
use warnings;

my $arg1 = $ARGV[0];
my $open1 = $arg1.".txt";
open LISTi, "<".$open1 or die "no open1";
while(<LISTi>){
	chomp();
	my $arrayi = $_;		
	open LISTj, "<".$open1 or die "no open1"; 
	while(<LISTj>){
		my $it=0;
		chomp();
		my $arrayj = $_;
		my %file1;
		my $f1 = $arrayi;
		open F1, "<".$arrayi or die "no f1";
		while(<F1>){
			chomp();
			my @l = split(/\t/,$_);			#array @l ohne ","
			my $key1 = $l[30]."_".$l[31]."_".$l[32]."_".$l[33];		#hashkey von Spalte 1 bis 3
			$file1{$key1} = join("\t",$l[26],$l[27],$l[39],$l[40]);	#$file1{$key} von Spalte 11 bis 16
		}
		close(F1);

		my %file2;
		my $f2 = $arrayj;
		open F2, "<".$f2 or die "no f2";
		while(<F2>){
			chomp();
			my @l = split(/\t/,$_);			#array @l ohne ","
			my $key2 = $l[30]."_".$l[31]."_".$l[32]."_".$l[33];		#hashkey von Spalte 1 bis 3
			$file2{$key2} = join("\t",$l[26],$l[27],$l[39],$l[40]);	#$file1{$key} von Spalte 11 bis 16
		}
		close(F2);
		foreach my $key1 ( sort keys %file1){  	
			if ($key1 eq "gene_location1_gene_location2_gene_name1_gene_name2") {next;}
			else {
				foreach my $key2 ( sort keys %file2){  
					if ($key1 eq $key2){
						print $file1{$key1}."\n";
						++$it;
					}
				}	
			}
		}
		print "\t\t $arrayi/$arrayj: $it match(es).\n\n";
	}
	close(LISTj);
}
close(LISTi);

