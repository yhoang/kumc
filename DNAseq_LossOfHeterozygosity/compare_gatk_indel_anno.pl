#!/usr/bin/perl -w
use strict;

print STDERR "./compare_gatk_indel_anno.pl \n";


#my $threshold = "_pass";
my $threshold = "_q30cov10";
my $covThr = 10;
#my $covThr = 5;

my $method = "bwamem";

my @NORMAL = ("N_L","N_M","N_M","N_N","N_N");
my @TUMOR = ("T_L","T_M1","T_M2","T_N1","T_N2");

for ( my $k = 0; $k <=$#NORMAL ; $k++) {

	my @file;
	$file[0] = "</project/results/Belinsky/$TUMOR[$k]/genotyper/anno/$TUMOR[$k].$method.indel$threshold.genome_summary.csv";
	$file[1] = "</project/results/Belinsky/$NORMAL[$k]/genotyper/anno/$NORMAL[$k].$method.indel$threshold.genome_summary.csv";
	$file[2] = ">/project/results/Belinsky/InDel_comparison/$TUMOR[$k]_$NORMAL[$k].$method.anno$threshold\_con";
	$file[3] = ">/project/results/Belinsky/InDel_comparison/$TUMOR[$k]_$NORMAL[$k].$method.anno$threshold\_dis";
	$file[4] = "</project/results/Belinsky/$TUMOR[$k]/$TUMOR[$k].cns.pileup";
	$file[5] = "</project/results/Belinsky/$NORMAL[$k]/$NORMAL[$k].cns.pileup";

	my %TUMOR_indel, my %w_tumor;
	my $it_tumor = 0;
	print "opens $file[0]...\n";
	open TUMOR_FILE, $file[0] or die "no $file[0]";
	while(<TUMOR_FILE>){
		chomp();
		my ( $chr, $start, $end, $ref, $indel, $het, $cov, $indelfreq, $func, $gene, $exfunc, $AAchange, $thousand, $dbSNP ) = qw/0 0 0 0 0 0 0 0 0 0 0 0 0/;
		my @l = split(/\t/,$_);			# array @l w/o tab
		# col00 Func,Gene,ExonicFunc,AAChange,
		# col04 UNIMPORTANT Conserved,SegDup,ESP6500_ALL,
		# col07 1000g2012feb_ALL,dbSNP135,
		# col09 UNIMPORTANT AVSIFT,LJB_PhyloP,LJB_PhyloP_Pred,LJB_SIFT,LJB_SIFT_Pred,LJB_PolyPhen2,LJB_PolyPhen2_Pred,LJB_LRT,LJB_LRT_Pred,LJB_MutationTaster,LJB_MutationTaster_Pred,LJB_GERP++,
		# col21 Chr,Start,End,Ref,Obs,Zygosity,snpQ_sum,coverage_per_sample
		# col29 mapQ,snpQ,varFreq,varFreq_mean,varCount
		if ($l[0] eq "Func") {next;}
		$chr = $l[21];
		$start = $l[22];
		$end = $l[22];
		$ref = $l[24];
		$indel = $l[25];
		$cov = $l[28];
		$indelfreq = $l[31];
		$func = $l[0];
		$gene = $l[1];
		$exfunc = $l[2];
		$AAchange = $l[3];
		$thousand = $l[7];
		$dbSNP = $l[8];
		
		$TUMOR_indel{$chr}{$start} = $indel;
		$w_tumor{$chr}{$start} = join("\t",$end,$ref,$indel,$cov,$indelfreq,$thousand,$dbSNP,$func,$gene,$exfunc,$AAchange);	# info
		++$it_tumor;
	}
	close(TUMOR_FILE);

	my %NORMAL_indel; my %w_normal; my %w_normalHelp;
	my $it_normal = 0;
	print "opens $file[1]...\n";
	open NORMAL_FILE, $file[1] or die "no $file[1]";
	while(<NORMAL_FILE>){
		chomp();
		my ( $chr, $start, $end, $ref, $indel, $cov, $indelfreq, $func, $gene, $exfunc, $AAchange, $thousand, $dbSNP ) = qw/0 0 0 0 0 0 0 0 0 0 0 0 0/;
		my @l = split(/\t/,$_);			# array @l w/o tab
		# col00 Func,Gene,ExonicFunc,AAChange,
		# col04 UNIMPORTANT Conserved,SegDup,ESP6500_ALL,
		# col07 1000g2012feb_ALL,dbSNP135,
		# col09 UNIMPORTANT AVSIFT,LJB_PhyloP,LJB_PhyloP_Pred,LJB_SIFT,LJB_SIFT_Pred,LJB_PolyPhen2,LJB_PolyPhen2_Pred,LJB_LRT,LJB_LRT_Pred,LJB_MutationTaster,LJB_MutationTaster_Pred,LJB_GERP++,
		# col21 Chr,Start,End,Ref,Obs,Zygosity,snpQ_sum,coverage_per_sample
		# col29 mapQ,snpQ,varFreq,varFreq_mean,varCount
		if ($l[0] eq "Func") {next;}
		$chr = $l[21];
		$start = $l[22];
		$end = $l[22];
		$ref = $l[24];
		$indel = $l[25];
		$cov = $l[28];
		$indelfreq = $l[31];
		$func = $l[0];
		$gene = $l[1];
		$exfunc = $l[2];
		$AAchange = $l[3];
		$thousand = $l[7];
		$dbSNP = $l[8];
	
		$NORMAL_indel{$chr}{$start} = $indel;
		$w_normal{$chr}{$start} = join("\t",$ref,$indel,$cov,$indelfreq);	# info
		$w_normalHelp{$chr}{$start} = join("\t",$thousand,$dbSNP,$func,$gene,$exfunc,$AAchange,$ref,$indel,$cov,$indelfreq);
		++$it_normal;
	}
	close(NORMAL_FILE);

	my %compareTUMOR;
	print "opens $file[4]...\n";
	open PILEUP_TUMOR, $file[4] or die "no $file[4]";
	while(<PILEUP_TUMOR>){
		chomp();
		my @l = split(/\t/,$_);			# array @l w/o tab
		my $chr = $l[0];
		my $start = $l[1];
		my $covTUMOR = $l[7];
		if ( exists $NORMAL_indel{$chr}{$start} && $covTUMOR >= $covThr ) {
			my $ref = $l[2];
			my $indel = $l[3];
			$compareTUMOR{$chr}{$start} = $indel."\t".$covTUMOR;
		}
	}
	close(PILEUP_TUMOR);

	my %compareNORMAL;
	print "opens $file[5]...\n";
	open PILEUP_NORMAL, $file[5] or die "no $file[5]";
	while(<PILEUP_NORMAL>){
		chomp();
		my @l = split(/\t/,$_);			# array @l w/o tab
		my $chr = $l[0];
		my $start = $l[1];
		my $covNORMAL = $l[7];
		if ( exists $TUMOR_indel{$chr}{$start} && $covNORMAL >= $covThr ) {
			my $ref = $l[2];
			my $indel = $l[3];
			$compareNORMAL{$chr}{$start} = $indel."\t".$covNORMAL;
		}
	}
	close(PILEUP_NORMAL);

	my %con, my %dis, my %dis1, my %dis2;
	my $con = 0; my $dis1 = 0, my $dis2 = 0; my $posit = 0; my $dis = 0;

	foreach my $chr ( sort keys %TUMOR_indel){  	# $chr: chr \t pos \t ref	
		if (defined $NORMAL_indel{$chr}){
			foreach my $start (sort keys %{$TUMOR_indel{$chr}}){
				if (defined $NORMAL_indel{$chr}{$start}) {
					++$posit;
					if ($TUMOR_indel{$chr}{$start} eq $NORMAL_indel{$chr}{$start}){
						$con{$chr}{$start} = join ("\t",$w_tumor{$chr}{$start},$w_normal{$chr}{$start});
						++$con;
					}
					else {
						$dis{$chr}{$start} = join ("\t",$w_tumor{$chr}{$start},$w_normal{$chr}{$start});
						++$dis;
						++$dis1;
						++$dis2;
					}
				}
			}
		}
	}
	foreach my $chr ( sort keys %TUMOR_indel){  	# $chr: chr \t pos \t ref	
		foreach my $start (sort keys %{$TUMOR_indel{$chr}}){
			if  (!exists $dis{$chr}{$start} && !exists $con{$chr}{$start} && exists $compareNORMAL{$chr}{$start} ) {
			$dis1{$chr}{$start} = join ("\t",$w_tumor{$chr}{$start},"\t".$compareNORMAL{$chr}{$start});
			++$dis1;
			}
		}
	}

	#foreach my $chr ( sort keys %NORMAL_indel){  	# $chr: chr \t pos \t ref
	#	foreach my $start (sort keys %{$NORMAL_indel{$chr}}){
	#		if  (!exists $dis{$chr}{$start} && !exists $con{$chr}{$start} && exists $compareTUMOR{$chr}{$start}) {
	#			$dis2{$chr}{$start} = join ("\t","\t\t".$compareTUMOR{$chr}{$start},"\t".$w_normalHelp{$chr}{$start});
	#			++$dis2;
	#		}
	#	}
	#}


	my $sum1 = $con + $dis1;
	my $sum2 = $con + $dis2;
	my $div1 = $con/$posit*100;
	my $div2 = $dis/$posit*100;

	open WRITE_CON, $file[2] or die;
	print "$con concordances out of $posit same positions ($it_tumor - $TUMOR[$k], $it_normal $NORMAL[$k])\nThis is $div1 percent!\n";
	print WRITE_CON "chr\tstart\tend\tref\t$TUMOR[$k]\tDepth\tInDel frequency\t1000g2012feb\tdbSNP135\tFunc\tGene\tExonicFunc\tAAChange\tref\t$NORMAL[$k]\tDepth\tInDel frequency\n";

	foreach my $chr ( sort keys %con ){
		foreach my $start (sort keys %{$con{$chr}} ){
			print WRITE_CON $chr."\t".$start."\t".$con{$chr}{$start}."\n";
		}
	}
	close (WRITE_CON);

	open WRITE_DIS, $file[3] or die;
	print "$dis ($dis1/$dis2) discordances out of $posit same positions ($it_tumor - $TUMOR[$k], $it_normal $NORMAL[$k])\nThis is $div2 percent!\n";
	print WRITE_DIS "chr\tstart\tend\tref\t$TUMOR[$k]\tDepth\tInDel frequency\t1000g2012feb\tdbSNP135\tFunc\tGene\tExonicFunc\tAAChange\tref\t$NORMAL[$k]\tDepth\tInDel frequency\n";

	foreach my $chr ( sort keys %dis ){
		foreach my $start (sort keys %{$dis{$chr}} ){
			print WRITE_DIS $chr."\t".$start."\t".$dis{$chr}{$start}."\n";
		}
	} print WRITE_DIS "\n";
	foreach my $chr ( sort keys %dis1 ){
		foreach my $start (sort keys %{$dis1{$chr}} ){
			print WRITE_DIS $chr."\t".$start."\t".$dis1{$chr}{$start}."\n";
		}
	}# print WRITE_DIS "\n";
	#foreach my $chr ( sort keys %dis2 ){
	#	foreach my $start (sort keys %{$dis2{$chr}} ){
	#		print WRITE_DIS $chr."\t".$start."\t".$dis2{$chr}{$start}."\n";
	#	}
	#}
	#close(WRITE_DIS);
	print "Written in $file[2] and $file[3]\n\n";
}
print "Done. \n";
