#!/usr/bin/perl -w
use strict;

print STDERR "./compare_gatk_snp_anno.pl \n";


my $threshold = "_pass";
#my $threshold = "_q43cov13";
#my $covThr = 13;
my $covThr = 10;

my $method = "bwamem";

#my @NORMAL = ("N_L","N_M","N_M","N_N","N_N","N_R");
#my @TUMOR = ("T_L","T_M1","T_M2","T_N1","T_N2","T_R");

my @NORMAL = ("N_Q");
my @TUMOR = ("T_Q");

print "sample\ttotal_tumor\ttotal_normal\tsame positions\tconcordant\tdiscordant\ttumor_not_normal\n";
for ( my $k = 0; $k <=$#NORMAL ; $k++) {
	my @file;
	$file[0] = "</project/results/Belinsky/$TUMOR[$k]/genotyper/anno/$TUMOR[$k].$method.snp$threshold.genome_summary.csv";
	$file[1] = "</project/results/Belinsky/$NORMAL[$k]/genotyper/anno/$NORMAL[$k].$method.snp$threshold.genome_summary.csv";
	$file[2] = ">/project/results/Belinsky/SNV_comparison/$TUMOR[$k]_$NORMAL[$k].$method.anno$threshold\_con";
	$file[3] = ">/project/results/Belinsky/SNV_comparison/$TUMOR[$k]_$NORMAL[$k].$method.anno$threshold\_dis";
#	$file[4] = "</project/results/Belinsky/$TUMOR[$k]/$TUMOR[$k].cns.starless.pileup";
	$file[5] = "</project/results/Belinsky/$NORMAL[$k]/$NORMAL[$k].cns.starless.pileup";

	my (%TUMOR, %TUMOR_snp, %w_tumor, %NORMAL, %NORMAL_snp, %w_normal, %w_normalHelp);
	my $it_tumor = 0;
	#print "opens $file[0]...\n";
	open TUMOR_FILE, $file[0] or die "no $file[0]";
	while(<TUMOR_FILE>){
		chomp();
		my ( $chr, $pos, $ref, $snp, $het, $cov, $snpfreq, $func, $gene, $exfunc, $AAchange, $thousand, $dbSNP, $sift ) = qw/0 0 0 0 0 0 0 0 0 0 0 0 0/;
		my @l = split(/\t/,$_);			# array @l w/o tab
		# col00 Func	Gene	ExonicFunc	AAChange	
		# col04 Conserved	SegDup	ESP6500_ALL	
		# col07 1000g2012feb_ALL	dbSNP135	AVSIFT	LJB_PhyloP	LJB_PhyloP_Pred	LJB_SIFT	LJB_SIFT_Pred	
		# col14 LJB_PolyPhen2	LJB_PolyPhen2_Pred	LJB_LRT	LJB_LRT_Pred	LJB_MutationTaster	LJB_MutationTaster_Pred	LJB_GERP++	
		# col21 Chr	Start	End	Ref	Obs	Zygosity	
		# col27 snpQ_sum	coverage_per_sample	mapQ	snpQ	varFreq	varFreq_mean	varCount
		if ($l[0] eq "Func") {next;}
		$chr = $l[21];
		$pos = $l[22];
		$ref = $l[24];
		$snp = $l[25];
		$het = $l[26];
		$cov = $l[28];
		$snpfreq = $l[31];
		$func = $l[0];
		$gene = $l[1];
		$exfunc = $l[2];
		$AAchange = $l[3];
		$thousand = $l[7];
		$dbSNP = $l[8];
		$sift = $l[13];
		
		$TUMOR_snp{$chr}{$pos} = $snp;
		$TUMOR{$chr}{$pos} = $snp."\t".$het;
		$w_tumor{$chr}{$pos} = join("\t",$ref,$snp,$het,$cov,$snpfreq,$thousand,$dbSNP,$sift,$func,$gene,$exfunc,$AAchange);	# info
		
	#	print "$TUMOR{$chr}{$pos}\n";
		++$it_tumor;
	}
	close(TUMOR_FILE);

	my $it_normal = 0;
#	print "opens $file[1]...\n";
	open NORMAL_FILE, $file[1] or die "no $file[1]";
	while(<NORMAL_FILE>){
		chomp();
		my ( $chr, $pos, $ref, $snp, $het, $cov, $snpfreq, $func, $gene, $exfunc, $AAchange, $thousand, $dbSNP, $sift ) = qw/0 0 0 0 0 0 0 0 0 0 0 0 0 0/;
		my @l = split(/\t/,$_);			# array @l w/o tab
		# col00 Func	Gene	ExonicFunc	AAChange	
		# col04 Conserved	SegDup	ESP6500_ALL	
		# col07 1000g2012feb_ALL	dbSNP135	AVSIFT	LJB_PhyloP	LJB_PhyloP_Pred	LJB_SIFT	LJB_SIFT_Pred	
		# col14 LJB_PolyPhen2	LJB_PolyPhen2_Pred	LJB_LRT	LJB_LRT_Pred	LJB_MutationTaster	LJB_MutationTaster_Pred	LJB_GERP++	
		# col21 Chr	Start	End	Ref	Obs	Zygosity	
		# col27 snpQ_sum	coverage_per_sample	mapQ	snpQ	varFreq	varFreq_mean	varCount
		if ($l[0] eq "Func") {next;}
		$chr = $l[21];
		$pos = $l[22];
		$ref = $l[24];
		$snp = $l[25];
		$het = $l[26];
		$cov = $l[28];
		$snpfreq = $l[31];
		$func = $l[0];
		$gene = $l[1];
		$exfunc = $l[2];
		$AAchange = $l[3];
		$thousand = $l[7];
		$dbSNP = $l[8];
		$sift = $l[13];
		
		$NORMAL_snp{$chr}{$pos} = $snp;
		$NORMAL{$chr}{$pos} = $snp."\t".$het;
		$w_normal{$chr}{$pos} = join("\t",$ref,$snp,$het,$cov,$snpfreq);	# info
		$w_normalHelp{$chr}{$pos} = join("\t",$thousand,$dbSNP,$sift,$func,$gene,$exfunc,$AAchange,$ref,$snp,$het,$cov,$snpfreq);
		
	#	print "$NORMAL{$chr}{$pos}\n";
		
		++$it_normal;
	}
	close(NORMAL_FILE);

	my %compareNORMAL;
#	print "opens $file[5]...\n";
	open PILEUP_NORMAL, $file[5] or die "no $file[5]";
	while(<PILEUP_NORMAL>){
		chomp();
		my @l = split(/\t/,$_);			# array @l w/o tab
		my $chr = $l[0];
		my $pos = $l[1];
		my $covNORMAL = $l[7];
		if ( exists $TUMOR_snp{$chr}{$pos} && $covNORMAL >= $covThr ) {
			my $ref = $l[2];
			$compareNORMAL{$chr}{$pos} = $ref."\t\t\t".$covNORMAL;
		}
	}
	close(PILEUP_NORMAL);

	my %con, my %dis, my %dis1, my %dis2;
	my $con = 0; my $dis1 = 0, my $dis2 = 0; my $posit = 0; my $dis = 0;

	foreach my $chr ( sort keys %TUMOR_snp){  	# $chr: chr \t pos \t ref	
		if (defined $NORMAL_snp{$chr}){
			foreach my $pos (sort keys %{$TUMOR_snp{$chr}}){
				if (defined $NORMAL_snp{$chr}{$pos}) {
					++$posit;
					if ($TUMOR{$chr}{$pos} eq $NORMAL{$chr}{$pos}){
						$con{$chr}{$pos} = join ("\t",$w_tumor{$chr}{$pos},$w_normal{$chr}{$pos});
						++$con;
					}
					else {
						$dis{$chr}{$pos} = join ("\t",$w_tumor{$chr}{$pos},$w_normal{$chr}{$pos});
						++$dis;
					}
				}
			}
		}
	}
	foreach my $chr ( sort keys %TUMOR_snp){  	# $chr: chr \t pos \t ref	
		foreach my $pos (sort keys %{$TUMOR_snp{$chr}}){
			if  (!exists $dis{$chr}{$pos} && !exists $con{$chr}{$pos} && exists $compareNORMAL{$chr}{$pos} ) {
			$dis1{$chr}{$pos} = join ("\t",$w_tumor{$chr}{$pos},$compareNORMAL{$chr}{$pos});
			++$dis1;
			}
		}
	}

	print "$TUMOR[$k]_$NORMAL[$k]\t$it_tumor\t$it_normal\t$posit\t$con\t$dis\t$dis1\n";
	open WRITE_CON, $file[2] or die;
	print WRITE_CON "chr\tstart\tref\t$TUMOR[$k]\tHeterozygosity\tDepth\tSNV frequency\t1000g2012feb\tdbSNP135\tAV_SIFT\tFunc\tGene\tExonicFunc\tAAChange\tref\t$NORMAL[$k]\tHeterozygosity\tDepth\tSNV frequency\n";

	foreach my $chr ( sort keys %con ){
		foreach my $pos (sort keys %{$con{$chr}} ){
			print WRITE_CON $chr."\t".$pos."\t".$con{$chr}{$pos}."\n";
		}
	}
	close (WRITE_CON);

	open WRITE_DIS, $file[3] or die;
	print WRITE_DIS "chr\tstart\tref\t$TUMOR[$k]\tHeterozygosity\tDepth\tSNV frequency\t1000g2012feb\tdbSNP135\tAV_SIFT\tFunc\tGene\tExonicFunc\tAAChange\tref\t$NORMAL[$k]\tHeterozygosity\tDepth\tSNV frequency\n";

	foreach my $chr ( sort keys %dis ){
		foreach my $pos (sort keys %{$dis{$chr}} ){
			print WRITE_DIS $chr."\t".$pos."\t".$dis{$chr}{$pos}."\n";
		}
	} print WRITE_DIS "\n";
	foreach my $chr ( sort keys %dis1 ){
		foreach my $pos (sort keys %{$dis1{$chr}} ){
			print WRITE_DIS $chr."\t".$pos."\t".$dis1{$chr}{$pos}."\n";
		}
	}
#	print "Written in $file[2] and $file[3]\n\n";
}
print "Done. \n";
