#!/usr/bin/perl -w
# Author Yen Hoang, M.Sc.
# 2013

use strict;
use POSIX qw/strftime/;

# Get pileup distributions (A/T/C/G) of targeted FF and correlated FFPE sample in one file.
# Forward/Reverse mapQuality added.

sub is_integer {
   defined $_[0] && $_[0] =~ /^[+-]?\d+$/;
}

print strftime "%Y-%m-%d %H:%M:%S", localtime(time);
print "\n";
print "./TES_get_distribution_of_intersection_with_forrev_mapQ.pl\n--> Get positions that were in both files! \n";

our @sample = ("2474","2561","2640","2685","2938","3050","3356","4079","4191","CN","CT","DN","DT");
#our @sample = ("CN","CT","DN","DT");

my ($FF, $FFPE, $write_distribution, $targeted, $title);

$targeted = "</project/results/me/TES/targeted_positions.list";

my %bed;
my $it_bed = 0;

print STDERR "opens $targeted... ";
open F2, $targeted or die "no $targeted";
while(<F2>){	#array data
	chomp();
	my @b = split(/\t/,$_);		# array @l w/o tab
	if ($b[0] eq "chr") {next;}
	my $chr = $b[0];
	my $pos = $b[1];
	my $ref = $b[2];
	$bed{$chr}{$pos} = $ref;
	++$it_bed;
}
close(F2);
print STDERR "$it_bed targeted positions.\n";


$title = "chr\tpos\tref\tA\tT\tG\tC\ta\tt\tg\tc\tmapQ\tA\tT\tG\tC\ta\tt\tg\tc\tmapQ\n";

for (my $o = 0; $o <= $#sample; $o++) {

        $FF = "</project/results/me/TES/$sample[$o]_FF/$sample[$o]_FF.cns.starless.pileup";	
        $FFPE = "</project/results/me/TES/$sample[$o]_FFPE/$sample[$o]_FFPE.cns.starless.pileup";
        $write_distribution = ">/project/results/me/TES/pileup_distribution/$sample[$o]_FF_FFPE_targeted_intersect.distribution"; 

	my %FF; my @read = ();
	my ($it_FF, $it_FF, $it_FFPE, $intersect) = qw/0 0 0 0/;	

	print "opens $FF...\n";
	open F1, $FF or die "no $FF";
	while(<F1>){
		chomp();
		my ($chr, $pos, $ref, $mapQ, $cov, $ten, $single) = qw/0 0 0 0 0 0 0/;
		my (@l, @distribution);
		
		@l = split(/\t/,$_);
		$chr = $l[0];
		$pos = $l[1];
		$ref = $l[2];
		$mapQ = $l[6];
		$cov = $l[7];
		
		if ( exists $bed{$chr}{$pos} ) {
        	# distribution(AaTtGgCc)
        	for (my $i = 0; $i <= 7; $i++) {
        	        $distribution[$i] = 0;
        	}
        	
        	@read = split(undef,$l[8]);
        	for (my $i = 0; $i <= $#read; $i++) {
        	        if ($read[$i] eq "+" || $read[$i] eq "-") {
        	        	if ( $read[$i-1] ne "^" ) {
			                $ten = 0;
			                $single = 0;
			                if ( is_integer($read[$i+2]) ) {
			                        $ten = $read[$i+1]*10; 
			                        $single = $read[$i+2]; 
			                        $i +=$ten; 
			                        $i +=$single; 
			                        $i +=2;
			                        next;
			                } else { 
			                        $single = $read[$i+1];
		        	                $i +=$single;
		        	                $i +=1;
		        	                next; 
		        	        }
		  		}
                	}
                	elsif ($read[$i] eq ".") {
                	        if ($ref eq "A") {++$distribution[0];}
                	        elsif ($ref eq "T") {++$distribution[1];}
                	        elsif ($ref eq "G") {++$distribution[2];}
                	        elsif ($ref eq "C") {++$distribution[3];}
                	}
                	elsif ($read[$i] eq ",") {
                	        if ($ref eq "A") {++$distribution[4];}
                	        elsif ($ref eq "T") {++$distribution[5];}
                	        elsif ($ref eq "G") {++$distribution[6];}
                	        elsif ($ref eq "C") {++$distribution[7];}
                	}
                	elsif ($read[$i] eq "A") {++$distribution[0];}
                	elsif ($read[$i] eq "T") {++$distribution[1];}
                	elsif ($read[$i] eq "G") {++$distribution[2];}
                	elsif ($read[$i] eq "C") {++$distribution[3];}
                	elsif ($read[$i] eq "a") {++$distribution[4];}
                	elsif ($read[$i] eq "t") {++$distribution[5];}
                	elsif ($read[$i] eq "g") {++$distribution[6];}
                	elsif ($read[$i] eq "c") {++$distribution[7];}
       		}	
		$FF{$chr}{$pos} = $ref;
		for (my $i = 0; $i <= 7; $i++) {
	        	$FF{$chr}{$pos} = $FF{$chr}{$pos}."\t".$distribution[$i]; 
		}
        	$FF{$chr}{$pos} = $FF{$chr}{$pos}."\t".$mapQ;
		
		++$it_FF;
	}}
	close(F1);

	open WRITE1, $write_distribution or die;
	print WRITE1 "$title";
	
	print "opens $FFPE...\n";
	open F1, $FFPE or die "no $FFPE";
	while(<F1>){
		my ($chr, $pos, $ref, $mapQ, $cov, $ten, $single) = qw/0 0 0 0 0 0 0/;
		my (@l, @distribution);
		chomp();
		@l = split(/\t/,$_);
		$chr = $l[0];
		$pos = $l[1];
		if (exists $FF{$chr}{$pos}) {
			$ref = $l[2];
			$mapQ = $l[6];
			$cov = $l[7];

        		# distribution(AaTtGgCc)
		        for (my $i = 0; $i <= 7; $i++) {
        		        $distribution[$i] = 0;
        		}

	        	@read = split(undef,$l[8]);
	        	for (my $i=0; $i<=$#read; $i++) {
	        	        if ($read[$i] eq "+" || $read[$i] eq "-") {
	        	                if ( $read[$i-1] ne "^" ) {
					        $ten = 0;
					        $single = 0;
					        if ( is_integer($read[$i+2]) ) {
				                        $ten = $read[$i+1]*10; 
				                        $single = $read[$i+2]; 
				                        $i +=$ten; 
				                        $i +=$single; 
				                        $i +=2;
				                        next;
				                } else { 
				                        $single = $read[$i+1];
				                        $i +=$single;
				                        $i +=1;
				                        next; 
				                }
					}
	        	        }
	               		elsif ($read[$i] eq ".") {
	                	        if ($ref eq "A") {++$distribution[0];}
	                	        elsif ($ref eq "T") {++$distribution[1];}
	                	        elsif ($ref eq "G") {++$distribution[2];}
	                	        elsif ($ref eq "C") {++$distribution[3];}
	                	}
	                	elsif ($read[$i] eq ",") {
	                	        if ($ref eq "A") {++$distribution[4];}
	                	        elsif ($ref eq "T") {++$distribution[5];}
        	        	        elsif ($ref eq "G") {++$distribution[6];}
        	        	        elsif ($ref eq "C") {++$distribution[7];}
	                	}
	                	elsif ($read[$i] eq "A") {++$distribution[0];}
		                elsif ($read[$i] eq "T") {++$distribution[1];}
	        	        elsif ($read[$i] eq "G") {++$distribution[2];}
	        	        elsif ($read[$i] eq "C") {++$distribution[3];}
	        	        elsif ($read[$i] eq "a") {++$distribution[4];}
	        	        elsif ($read[$i] eq "t") {++$distribution[5];}
	        	        elsif ($read[$i] eq "g") {++$distribution[6];}
	        	        elsif ($read[$i] eq "c") {++$distribution[7];}
	        	}
                        ++$intersect;
                        print WRITE1 "$chr\t$pos\t$FF{$chr}{$pos}";
	                for (my $i = 0; $i <= 7; $i++) {
	                        print WRITE1 "\t$distribution[$i]"; 
	                }
	                print WRITE1 "\t$mapQ\n";
		}
		++$it_FFPE;
	}
	close(F1);
	close(WRITE1);
	
	print "Written $intersect intersect positions ($it_FF/$it_FFPE) in $write_distribution.\n\n";
}

print "Done. \n";
print strftime "%Y-%m-%d %H:%M:%S", localtime(time);
print "\n";
