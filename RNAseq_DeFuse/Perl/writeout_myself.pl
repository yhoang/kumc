#!/usr/bin/perl -w

#chmod +x writeout_myselfpl
#./writeout_myself.pl samples
#output >compare_myself/samplex.txt

use strict;
use warnings;

my @memory;
my $i     = 0;
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
		my $key1 = $l[30]."_".$l[31]."_".$l[32]."_".$l[33];    
		if ( $key1 eq "gene_location1_gene_location2_gene_name1_gene_name2" ){
			next;
		}
		elsif ( $key1 eq "\"gene_location1\"_\"gene_location2\"_\"gene_name1\"_\"gene_name2\"" ){		
			next;
		}
		else {
			if ($l[26]eq "MT" || $l[26] eq "\"MT\""){$l[26]=0;}
			if ($l[26]>$l[27]){
			$file1{$key1} = join( "\t",$l[27],$l[26],$l[29],$l[28]); 
			}
			else{
			$file1{$key1} = join( "\t",$l[26],$l[27],$l[28],$l[29]); 
			}	
			print WRITER $file1{$key1}."\n";
		}
	}
	close(F1);

	close(WRITER);
}
close(LISTi);

