#!/usr/bin/perl -w

#chmod +x compare_1on1_automated.pl
#./compare_1on1_automated.pl samples >samples_1on1.txt

use strict;
use warnings;

my @memory;
my $i     = 0;
my $open1 = $ARGV[0].".txt";

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
			print $arrayi."_".$arrayj."\n";
			
			my %file1;
			my $f1 = $arrayi;
			open F1, "<samples/" . $arrayi or die "no f1";
			while (<F1>) {
				chomp();
				my @l = split( /\t/, $_ );    #array @l without tabs
				my $key1 = $l[30]."_".$l[31]."_".$l[32]."_".$l[33]; 		# hashkey loc1 loc2 name1 name2
				$file1{$key1} = join( "\t",$l[26],$l[27],$l[28],$l[29]);	# hash chr1 chr2 start1 start2
			}
			close(F1);

			my %file2;
			my $f2 = $arrayj;
			open F2, "<samples/" . $f2 or die "no f2";
			while (<F2>) {
				chomp();
				my @l = split( /\t/, $_ );    #array @l without tabs
				my $key2 = $l[30]."_".$l[31]."_".$l[32]."_".$l[33];			# hashkey loc1 loc2 name1 name2
				$file2{$key2} = join( "\t",$l[26],$l[27],$l[28],$l[29]);	# hash chr1 chr2 start1 start2
			}
			close(F2);
	
			open WRITER, ">compare_1on1/".$arrayi."_".$arrayj or die $!;
				foreach my $key1 ( sort keys %file1 ) {
					if ( $key1 eq
						"gene_location1_gene_location2_gene_name1_gene_name2" ){
						next;
					}
					elsif ( $key1 eq "\"gene_location1\"_\"gene_location2\"_\"gene_name1\"_\"gene_name2\"" ){		
						next;
					}
					else {
						foreach my $key2 ( sort keys %file2 ) {
							if ( $key1 eq $key2 ) {
								print WRITER $file1{$key1} . "\n";
								++$it;
							}
						}
					}
				}
			close(WRITER);	
		}
	}
	close(LISTj);
}
close(LISTi);
