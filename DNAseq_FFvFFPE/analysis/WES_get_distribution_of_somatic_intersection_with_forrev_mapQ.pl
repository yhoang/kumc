#!/usr/bin/perl -w
# Author Yen Hoang, M.Sc.
# 2013

use strict;
use POSIX qw/strftime/;

# Get pileup distributions (A/T/C/G) of somatic FF and FFPE sample in one file, respectively.
# Forward/Reverse mapQuality added.

sub is_integer {
   defined $_[0] && $_[0] =~ /^[+-]?\d+$/;
}
print strftime "%Y-%m-%d %H:%M:%S", localtime(time);
print "\n";
print "./WES_get_distribution_of_somatic_intersection_with_forrev_mapQ.pl\n--> Get positions that were in both files! \n";

my @sample = ("C","D");
my @sample2 = ("FF","FFPE");

my @samplesize = (656871468,421718167,187309470,137539076); # CN_FF, CN_FFPE, DN_FF, DN_FFPE

my ($FF, $FFPE, $write_distribution, $pointerPos);

my $title = "chr\tpos\tref\tA\tT\tG\tC\ta\tt\tg\tc\tmapQ\tA\tT\tG\tC\ta\tt\tg\tc\tmapQ\n";

for (my $o = 0; $o <= $#sample; $o++) {
		for (my $q = 0; $q <= $#sample2; $q++) {
		print strftime "%Y-%m-%d %H:%M:%S\n", localtime(time);
		$FF = "</project/results/me/WES/$sample[$o]N_$sample2[$q]ex/$sample[$o]N_$sample2[$q]ex.cns.starless.pileup";	
		$FFPE = "</project/results/me/WES/$sample[$o]T_$sample2[$q]ex/$sample[$o]T_$sample2[$q]ex.cns.starless.pileup";
		$write_distribution = ">/project/results/me/WES/pileup_distribution/$sample[$o]N_$sample[$o]T_$sample2[$q]ex_intersect.distribution"; 
	
		open WRITE1, $write_distribution or die;
		print WRITE1 "$title";
		close (WRITE1);
	
		my ($it_FF, $it_FFPE, $intersect, $pointerPos, $pointerHelp, $loop, $it_loop) = qw/0 0 0 0 0 1 0/;	
	
		do {
			++$it_loop;
		
			my %FF; my @read = ();
			$pointerHelp = 0;
			$it_FFPE = 0;
		
			print "Run $it_loop: opens $FF, point on line $it_FF/byte $pointerPos...";
			open FF_FILE, $FF or die "no $FF";
			seek FF_FILE, $pointerPos, 0;
			while(<FF_FILE>){
				chomp();
				my ($chr, $pos, $ref, $mapQ, $cov, $ten, $single) = qw/0 0 0 0 0 0 0/;
				my (@l, @distribution);
			
				@l = split(/\t/,$_);
				$chr = $l[0];
				$pos = $l[1];
				$ref = $l[2];
				$mapQ = $l[6];
				$cov = $l[7];
			
				$mapQ == 0 and next;
				
				# distribution(AaTtGgCc)
				for (my $i = 0; $i <= 7; $i++) {
				        $distribution[$i] = 0;
				}
				
				@read = split(undef,$l[8]);
				for (my $i = 0; $i <= $#read; $i++) {
				        if ( $read[$i] eq "+" || $read[$i] eq "-" ) {
				        	if ( $read[$i-1] ne "^" ) {
						        $ten = 0;
						        $single = 0;
						        if ( is_integer($read[$i+2]) ) {
						       # if ( ($read[$i+2]*1) eq $read[$i+2]) {
						                $ten = $read[$i+1]*10; 
						                $single = $read[$i+2]; 
						                $i +=$ten; 
						                $i +=$single; 
						                $i +=2;
						        } else { 
						                $single = $read[$i+1];
						                $i +=$single;
						                $i +=1;
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
				$FF{$chr}{$pos} = $ref;#."\t".$var;
				for (my $i = 0; $i <= 7; $i++) {
					$FF{$chr}{$pos} = $FF{$chr}{$pos}."\t".$distribution[$i]; 
				}
				$FF{$chr}{$pos} = $FF{$chr}{$pos}."\t".$mapQ;
			
				++$it_FF;
				++$pointerHelp;
			
				if ( $pointerHelp >= 20000000 ) {
					$pointerPos = tell FF_FILE;
					print "Buffer on line $it_FF/$samplesize[$o+$q] (byte $pointerPos) - Break!\n";
					last;
				} 
			
				if ( $it_FF >= $samplesize[$o+$q] ) {
					$loop = 0;
					$pointerPos = tell FF_FILE;
					print "Buffer on LAST line $it_FF/$samplesize[$o+$q] (byte $pointerPos)!\n";
					last;
				}
			}
			close(FF_FILE);
	
			open WRITE1, ">$write_distribution" or die;
			open FFPE_FILE, $FFPE or die "no $FFPE";
				while(<FFPE_FILE>){
				my ($chr, $pos, $ref, $mapQ, $cov, $ten, $single) = qw/0 0 0 0 0 0 0/;
				my (@l, @distribution);
				chomp();
				@l = split(/\t/,$_);
				$chr = $l[0];
				$pos = $l[1];
				$ref = $l[2];
				$mapQ = $l[6];
				$mapQ == 0 and next;
				if (exists $FF{$chr}{$pos}) {
					$cov = $l[7];

					# distribution(AaTtGgCc)
					for (my $i = 0; $i <= 7; $i++) {
				        $distribution[$i] = 0;
					}

					@read = split(undef,$l[8]);
					for (my $i=0; $i<=$#read; $i++) {
						if ( $read[$i] eq "+" || $read[$i] eq "-" ) {
							if ( $read[$i-1] ne "^" ) {
								$ten = 0;
								$single = 0;
								if ( is_integer($read[$i+2]) ) {
								        $ten = $read[$i+1]*10; 
								        $single = $read[$i+2]; 
								        $i +=$ten; 
								        $i +=$single; 
								        $i +=2;
								} else { 
								        $single = $read[$i+1];
								        $i +=$single;
								        $i +=1;
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
			close(FFPE_FILE);
			close(WRITE1);
		
		} while ( $loop == 1 );
		
		print "Written $intersect intersect positions ($it_FF/$it_FFPE) in $write_distribution.\n\n";
	}
}
print "Done. \n";
print strftime "%Y-%m-%d %H:%M:%S", localtime(time);
print "\n";
