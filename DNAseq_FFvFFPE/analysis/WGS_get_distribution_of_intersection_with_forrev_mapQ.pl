#!/usr/bin/perl -w
# Author Yen Hoang, M.Sc.
# 2013

use strict;
use POSIX qw/strftime/;

# Get pileup distributions (A/T/C/G) of FF and correlated FFPE sample in one file.
# Forward/Reverse mapQuality added.

sub is_integer {
   defined $_[0] && $_[0] =~ /^[+-]?\d+$/;
}

print strftime "%Y-%m-%d %H:%M:%S", localtime(time);
print "\n";
print "./WGS_get_distribution_of_intersection_with_forrev_mapQ.pl\n--> Get positions that were in both files! \n";

my @sample = ("CN");
my @samplesize = (2843299750);
my ($FF, $FFPE, $write_distribution, $pointerPos);

my $title = "chr\tpos\tref\tA\tT\tG\tC\ta\tt\tg\tc\tmapQ\tA\tT\tG\tC\ta\tt\tg\tc\tmapQ\n";

for (my $o = 0; $o <= $#sample; $o++) {

	print strftime "%Y-%m-%d %H:%M:%S\n", localtime(time);
	$FF = "</project/results/me/WGS/$sample[$o]_FF/$sample[$o]_FF.cns.starless.pileup";	
        $FFPE = "</project/results/me/WGS/$sample[$o]_FFPE/$sample[$o]_FFPE.cns.starless.pileup";
        $write_distribution = ">/project/results/me/WGS/pileup_distribution/$sample[$o]_FF_FFPE_intersect.distribution2"; 
	
	open WRITE1, $write_distribution or die;
	print WRITE1 "$title";
	close (WRITE1);
	
	my ($it_FF, $it_FFPE, $intersect, $pointerPos, $pointerHelp, $loop, $it_loop, $pointerPos2, $pointerHelp2, $it_pointer2) = qw/0 0 0 0 0 1 0 0/;	
	my @pointer2, my @it_FFPE;
	
	do {
		++$it_loop;
		
		my %FF; my @read = ();
		$pointerHelp = 0;
		
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
				print "Buffer on line $it_FF/$samplesize[$o]_FF (byte $pointerPos) - Break!\n";
				last;
			} 
			
			if ( $it_FF >= $samplesize[$o] ) {
				$loop = 0;
				$pointerPos = tell FF_FILE;
				print "Buffer on LAST line $it_FF/$samplesize[$o]_FF (byte $pointerPos)!\n";
				last;
			}
		}
		close(FF_FILE);
	
		#print "opens $FFPE...\n";
		
		print "\tStill on Run $it_loop: opens $FFPE, point on byte $pointerPos2...\n";
		
		$it_FFPE[$it_loop] = 0;
		if ( $it_loop > 1 ) {
			$it_pointer2 = $it_pointer2 - 3;
			$pointerPos2 = $pointer2[$it_pointer2];
		}
		open WRITE1, ">$write_distribution" or die;
		open FFPE_FILE, $FFPE or die "no $FFPE";
		seek FFPE_FILE, $pointerPos2, 0;
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
			++$it_FFPE[$it_loop];
			
			if ( ($it_FFPE[$it_loop] % 20000000) == 0 ) { ## save position on every 20,000,000
				$pointerPos2 = tell FFPE_FILE;
				$pointer2[$it_pointer2] = $pointerPos2;
				print "\tSaved line $it_FFPE[$it_loop] ($it_pointer2. byte $pointerPos2)\n";
				++$it_pointer2;
			}
			
			if ($it_FFPE[$it_loop] == 3 * 20000000) {
				$pointerPos2 = tell FFPE_FILE;
				print "\tBuffer on line $it_FFPE[$it_loop] (byte $pointerPos2) - Break!\n";
				last;
			}
		}
		close(FFPE_FILE);
		close(WRITE1);
		
	} while ( $loop == 1 );
	
	my $sum_FFPE = ($it_pointer2-1)*20000000 + $it_FFPE[$it_loop];
	print "Written $intersect intersect positions ($it_FF/$sum_FFPE) in $write_distribution.\n\n";
}
print "Done. \n";
print strftime "%Y-%m-%d %H:%M:%S", localtime(time);
print "\n";
