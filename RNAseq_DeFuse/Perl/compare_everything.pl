#!/usr/bin/perl -w

#chmod +x compare_everything.pl
#./compare_everything.pl samples_number >dendogram.txt

use strict;
use warnings;

my @memory;
my $i     = 0;
my $noti=0;
my $j =0;
my $notj=0;
my $run = 0;
my $open1 = $ARGV[0].".txt";
my @counts;
open LISTi, "<" . $open1 or die "no open1";
while (<LISTi>) {
	chomp();
	my $arrayi = $_;
	++$i;
	while ($i < $arrayi){
	#print "$i$arrayi\t";
		for ($noti=1;$noti<=36;$noti++){
	#		print  "0\t";
			$counts[$i][$noti]=0;
		}	
	#	print "\n";
		++$i;
	}
	open LISTj, "<" . $open1 or die "no open1";
	while (<LISTj>) {
		my $it = 0;
		chomp();
		my $arrayj = $_;
		++$j;
		while ($j < $arrayj){
		#	print "$j$arrayj\t";
			$counts[$i][$j]=0;
			++$j;
		}
		my %file1;
		my $f1 = $arrayi;
		open F1, "<sample_number/" . $arrayi or die "no f1";
		while (<F1>) {
			chomp();
			my @l = split( /\t/, $_ );    #array @l ohne ","
			my $key1 = $l[30]."_".$l[31]."_".$l[32]."_".$l[33];                   #hashkey von Spalte 1 bis 3
			$file1{$key1} = join( "\t",$l[26],$l[27],$l[28],$l[29]);     #$file1{$key} von Spalte 11 bis 16
		}
		close(F1);

		my %file2;
		my $f2 = $arrayj;
		open F2, "<sample_number/" . $f2 or die "no f2";
		while (<F2>) {
			chomp();
			my @l = split( /\t/, $_ );    #array @l ohne ","
			my $key2 = $l[30]."_".$l[31]."_".$l[32]."_".$l[33];                   #hashkey von Spalte 1 bis 3
			$file2{$key2} = join( "\t",$l[26],$l[27],$l[28],$l[29]);    #$file1{$key} von Spalte 11 bis 16
		}
		close(F2);

		foreach my $key1 ( sort keys %file1 ) {
		if ( $key1 eq "gene_location1_gene_location2_gene_name1_gene_name2" ){
			next;
		}
		elsif ( $key1 eq "\"gene_location1\"_\"gene_location2\"_\"gene_name1\"_\"gene_name2\"" ){		
			next;
		}
		else {
			foreach my $key2 ( sort keys %file2 ) {
				if ( $key1 eq $key2 ) {
					++$it;
				}
			}
		}
	}
	$counts[$i][$j]=$it;
	}	
		#	print "@counts\n";
	close(LISTj);
	$j=0;
}
my $n=0;
my $m=0;
my $percent;
for ($n=1;$n<=36;$n++){
	for ($m=1;$m<=36;$m++){
		$percent =$counts[$n][$m]/$counts[$n][$n];
		print "$percent\t";
	}
	print"\n";
}

}
close(LISTi);


