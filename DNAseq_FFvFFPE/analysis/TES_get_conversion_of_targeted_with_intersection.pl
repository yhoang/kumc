#!/usr/bin/perl -w
# Author Yen Hoang, M.Sc.
# 2013

use strict;
use POSIX qw/strftime/;

# Conversion rate of SNVs in targeted intersected positions
# i.e. A>T 30%, A>G 30%, A>C 40%

print strftime "%Y-%m-%d %H:%M:%S", localtime(time);
print "\n";
print STDERR "./TES_get_conversion_of_targeted_with_intersection.pl\n";


our @sample = ("2474","2561","2640","2685","2938","3050","3356","4079","4191","CN","CT","DN","DT");
our $p = ("q43cov13");

our ($FF, $FFPE, $INTERSECT, $WRITE_CON, $WRITE_DIS, $TITLE, $targeted);
$TITLE = "chr\tpos\tref\tvar_FF\thet\tsnpQ\tDepth\tQualbyDepth\tRMS mapping quality\tFunc\tGene\tExonicFunc\tAAChange\tLJB_SIFT_Pred\tvar_FFPE\thet\tsnpQ\tDepth\tQualbyDepth\tRMS mapping quality\n";

$targeted = "</project/results/me/TES/targeted_positions.list";

my %bed;
my $it_bed = 0;

print STDERR "opens $targeted... \n";
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

my @nuc = ('A','T','G','C');
my ( @conversionFF, @conversionFFPE, @distributionFF, @distributionFFPE);

for ( my $i = 0; $i <=13; $i++) {
        for (my $j = 0; $j <=16; $j++) {
                $distributionFF[$i][$j] = 0;
                $distributionFFPE[$i][$j] = 0;
        }
}

my $title1 = "FF";
my $title2 = "\n\nFFPE";


for ( my $k = 0; $k <=$#sample; $k++) {
	$title1 = $title1."\t$sample[$k]";
	$title2 = $title2."\t$sample[$k]";

    $FF = "</project/results/me/TES/$sample[$k]_FF/genotyper/anno/$sample[$k]_FF.bwa.ug.snp_${p}.genome_summary.csv";
	$FFPE = "</project/results/me/TES/$sample[$k]_FFPE/genotyper/anno/$sample[$k]_FFPE.bwa.ug.snp_${p}.genome_summary.csv";
    $INTERSECT = "</project/results/me/TES/pileup_distribution/$sample[$k]_FF_FFPE_intersect.distribution";
	
	my (%FF, %FFPE);
	
    open FF_FILE, $FF or die "no $FF";
	while(<FF_FILE>){
    	chomp();
    	my @l = split(/\t/,$_);			# array @l w/o tab
    	# col00 Func,Gene,ExonicFunc,AAChange,
    	# col04 UNIMPORTANT Conserved,SegDup,ESP6500_ALL,
    	# col07 1000g2012feb_ALL,dbSNP135,
    	# col09 OTHERS: AVSIFT,LJB_PhyloP,LJB_PhyloP_Pred,LJB_SIFT,LJB_SIFT_Pred,LJB_PolyPhen2,LJB_PolyPhen2_Pred,
    	#               LJB_LRT,LJB_LRT_Pred,LJB_MutationTaster,LJB_MutationTaster_Pred,LJB_GERP++,
    	# col21 Chr,Start,End,Ref,Obs,Zygosity,snpQ_sum,coverage_per_sample
    	# col29 mapQ,snpQ,varFreq,varFreq_mean,varCount		
    	if ($l[0] eq "Func") {next;}
    	
    	my $chr = $l[21];
    	my $pos = $l[22];
    	
    	if (exists $bed{$chr}{$pos}) {
    		my $var = $l[25];
    		$FF{$chr}{$pos} = $var;
    	}
    }
	close(FF_FILE);
		
	open FFPE_FILE, $FFPE or die "no $FFPE";
	while(<FFPE_FILE>){
		chomp();
		my @l = split(/\t/,$_);			# array @l w/o tab
		# col00 Func,Gene,ExonicFunc,AAChange,
		# col04 UNIMPORTANT Conserved,SegDup,ESP6500_ALL,
		# col07 1000g2012feb_ALL,dbSNP135,
		# col09 OTHERS: AVSIFT,LJB_PhyloP,LJB_PhyloP_Pred,LJB_SIFT,LJB_SIFT_Pred,LJB_PolyPhen2,LJB_PolyPhen2_Pred,
		#               LJB_LRT,LJB_LRT_Pred,LJB_MutationTaster,LJB_MutationTaster_Pred,LJB_GERP++,
		# col21 Chr,Start,End,Ref,Obs,Zygosity,snpQ_sum,coverage_per_sample
		# col29 mapQ,snpQ,varFreq,varFreq_mean,varCount		
		if ($l[0] eq "Func") {next;}
		my $chr = $l[21];
		my $pos = $l[22];
		if (exists $bed{$chr}{$pos}) {
			my $var = $l[25];
			$FFPE{$chr}{$pos} = $var;
		}
	}
	close(FFPE_FILE);
	
	my ($it_FF, $it_FFPE) = qw/0 0/;
	open INTERSECT_FILE, $INTERSECT or die "no $INTERSECT";
	while(<INTERSECT_FILE>){	#array data
		chomp();
		my @b = split(/\t/,$_);		# array @l w/o tab
		if ($b[0] eq "chr") {next;}
		my $chr = $b[0];
		my $pos = $b[1];
		my $ref = $b[2];
		
		if ( exists $FF{$chr}{$pos} ) {
			$conversionFF[$it_FF][0] = $ref;
			$conversionFF[$it_FF][1] = $FF{$chr}{$pos};
			++$it_FF;
		}
		
		if ( exists $FFPE{$chr}{$pos} ) {
			$conversionFFPE[$it_FFPE][0] = $ref;
			$conversionFFPE[$it_FFPE][1] = $FFPE{$chr}{$pos};
			++$it_FFPE;
		}
	}
		
        for (my $j = 0; $j <=3; $j++) { # for nuc1
                for (my $i = 0; $i <= $it_FF; $i++) {      #for every discordant call
                        if ($conversionFF[$i][0] eq $nuc[$j]) {   
                                for (my $l = 0; $l <=3; $l++) { #for conversion in nuc2
                                        if ($conversionFF[$i][1] eq $nuc[$l]) {
                                                ++$distributionFF[$k][4*$j+$l];
                                        }
                                }
                        }
                }
        }
        $distributionFF[$k][16] = $it_FF;
        
        for (my $j = 0; $j <=3; $j++) { # for nuc1
                for (my $i = 0; $i <= $it_FFPE; $i++) {      #for every discordant call
                        if ($conversionFFPE[$i][0] eq $nuc[$j]) {   
                                for (my $l = 0; $l <=3; $l++) { #for conversion in nuc2
                                        if ($conversionFFPE[$i][1] eq $nuc[$l]) {
                                                ++$distributionFFPE[$k][4*$j+$l];
                                        }
                                }
                        }
                }
        }
        $distributionFFPE[$k][16] = $it_FFPE;
}

$title1 = $title1."\tmean\n";
$title2 = $title2."\tmean\n";

### print actual count FF
print "$title1";
for (my $i = 0; $i <= 3; $i++) {        #nuc1
        for (my $j = 0; $j <= 3; $j++) {#nuc2
                my $mean = 0;
                print "$nuc[$i]\>$nuc[$j]\t";
                for (my $s = 0; $s <= $#sample; $s++) {     #for sample number
                        print "$distributionFF[$s][4*$i+$j]\t";
                        $mean +=$distributionFF[$s][4*$i+$j];
                }
                $mean = $mean/($#sample+1);
                print "$mean\n";
        }
}
print "count";
for (my $s = 0; $s <= $#sample; $s++) {
        print "\t$distributionFF[$s][16]";
}


### print actual count FFPE
print "$title2";
for (my $i = 0; $i <= 3; $i++) {        #nuc1
        for (my $j = 0; $j <= 3; $j++) {#nuc2
                my $mean = 0;
                print "$nuc[$i]\>$nuc[$j]\t";
                for (my $s = 0; $s <= $#sample; $s++) {     #for sample number
                        print "$distributionFFPE[$s][4*$i+$j]\t";
                        $mean +=$distributionFFPE[$s][4*$i+$j];
                }
                $mean = $mean/($#sample+1);
                print "$mean\n";
        }
}
print "count";
for (my $s = 0; $s <= $#sample; $s++) {
        print "\t$distributionFFPE[$s][16]";
}


print "\nDone.\n"; 
print strftime "%Y-%m-%d %H:%M:%S", localtime(time);
print "\n";
