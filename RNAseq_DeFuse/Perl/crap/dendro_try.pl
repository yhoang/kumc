#!/usr/bin/perl -w

#chmod +x dendro.pl
#./dendro.pl samples

#use strict;
use warnings;

my @memory;
my $i     = 0;
my $j =0;
my $notj=0;
my $not_i =0;
my $open1 = $ARGV[0].".txt";
my %file1;
open LISTi, "<" . $open1 or die "no open1";
open WRITER, ">>dendro_try.txt" or die $!;
while (<LISTi>) {
	chomp();
	my $arrayi = $_;
	++$j;
	while ($j < $arrayi){
		print "$j $arrayi\n";
		for ($notj=0;$notj<44;$notj++){
			print WRITER "100\tNA\t";
		}	
		print WRITER "\n";
		++$j;
	}
	open F1, "<sample_number/" . $arrayi or die "no f1";
	while (<F1>) {
		chomp();
		++$i;
		my @l = split( /\t/, $_ );    #array @l ohne ","
		my $key1 = $l[30]."_".$l[31]."_".$l[32]."_".$l[33];    
#		if ($key1=~/D/){
		if ( $key1 eq "gene_location1_gene_location2_gene_name1_gene_name2" ){
			next;
		}
		elsif ( $key1 eq "\"gene_location1\"_\"gene_location2\"_\"gene_name1\"_\"gene_name2\"" ){		
			next;
		}
		else {
		#	$file1{$key1} = join( "\t",$l[26],$l[27],$l[39],$l[40]); 
		#	$file1{$key1} = join($l[26],$l[27],$l[39],$l[40]); 
			if ($l[26].$l[27] eq "MTMT") {$file1{$key1}=join("\t",0,0);}	
			elsif ($l[26].$l[27] eq "XX") {$file1{$key1}=join("\t",23,23);}		
			else {$file1{$key1} = join("\t",$l[26],$l[27]);} 
			print WRITER $file1{$key1}."\t";
		}
	}
	for ($not_i=44;$not_i>=$i;$not_i--){
		print WRITER "NA\tNA\t";
	}
	close(F1);
	$i=0;
	print WRITER "\n";
}
close(WRITER);
close(LISTi);

