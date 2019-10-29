#!/usr/bin/perl -w
# Author Yen Hoang, M.Sc.
# 2013

use strict;
use POSIX qw/strftime/;

# Compares annotated InDels between FF and FFPE in targeted and interesected positions.

print strftime "%Y-%m-%d %H:%M:%S", localtime(time);
print "\n";
print STDERR "./TES_compare_gatk_indel_anno.pl \n";

our @method = ("q43cov13");
our @covThr = (13);

our @sample = ("2474","2561","2640","2685","2938","3050","3356","4079","4191","CN","CT","DN","DT");

print "sample\tthreshold\ttotal_FF\ttotal_FFPE\tintersect_FF\tintersect_FFPE\tsame positions\tcon\tdis\tFN\tFP\t% con\t% con_FF\t% con_FFPE\n";

my $targeted = "</project/results/me/TES/targeted_positions.list";
my %bed;
my $it_bed = 0;

open F2, $targeted or die "no $targeted";
while(<F2>){	#array data
	chomp();
	my @b = split(/\t/,$_);		# array @l w/o tab
	if ($b[0] eq "chr") {next;}
	my $chr = $b[0];
	my $start = $b[1];
	my $ref = $b[2];
	$bed{$chr}{$start} = $ref;
	++$it_bed;
}
close(F2);

for ( my $o = 0; $o <=$#sample ; $o++) {
	for ( my $p = 0; $p <=$#method ; $p++) {
		my @file;
		$file[0] = "</project/results/me/TES/$sample[$o]_FF/genotyper/anno/$sample[$o]_FF.bwa.hapC.indel_$method[$p].genome_summary.csv";
		$file[1] = "</project/results/me/TES/$sample[$o]_FFPE/genotyper/anno/$sample[$o]_FFPE.bwa.hapC.indel_$method[$p].genome_summary.csv";
		$file[2] = ">/project/results/me/TES/InDel_comparison/targeted/$sample[$o]_FF_FFPE.hapC.indel.anno_$method[$p]_con";
		$file[3] = ">/project/results/me/TES/InDel_comparison/targeted/$sample[$o]_FF_FFPE.hapC.indel.anno_$method[$p]_dis";
		$file[4] = "</project/results/me/TES/$sample[$o]_FF/$sample[$o]_FF.cns.pileup";
		$file[5] = "</project/results/me/TES/$sample[$o]_FFPE/$sample[$o]_FFPE.cns.pileup";
		
		my ( $it_FF, $it_interFF, $it_FFPE, $it_interFFPE, $con, $dis, $FN, $FP, $posit, $sum1, $sum2, $div, $div1, $div2) = qw/0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0/; 
		my ( %FF_indel, %w_FF, %compareFF, %FFPE_indel, %w_FFPE, %w_FFPEHelp, %compareFFPE, %con, %dis, %FN, %FP);

		open FF_FILE, $file[0] or die "no $file[0]";
		while(<FF_FILE>){
			chomp();
			my ( $chr, $start, $end, $ref, $indel, $het, $cov, $indelfreq, $func, $gene, $exfunc, $AAchange, $thousand, $dbSNP ) = qw/0 0 0 0 0 0 0 0 0 0 0 0 0/;
			my @l = split(/\t/,$_);			# array @l w/o tab
			# col00 Func,Gene,ExonicFunc,AAChange,
			# col04 Conserved,SegDup,ESP6500_ALL,
			# col07 1000g2012feb_ALL,dbSNP135,
			# col09 OTHERS: AVSIFT,LJB_PhyloP,LJB_PhyloP_Pred,LJB_SIFT,LJB_SIFT_Pred,LJB_PolyPhen2,LJB_PolyPhen2_Pred,
			#               LJB_LRT,LJB_LRT_Pred,LJB_MutationTaster,LJB_MutationTaster_Pred,LJB_GERP++,
			# col21 Chr,Start,End,Ref,Obs,Zygosity,snpQ_sum,coverage_per_sample
			# col29 mapQ,snpQ,varFreq,varFreq_mean,varCount
			if ($l[0] eq "Func") {next;}
			$chr = $l[21];
			$start = $l[22];
			$end = $l[23];
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
			
			if ( exists $bed{$chr}{$start} ) {
				$FF_indel{$chr}{$start} = $indel;
				$w_FF{$chr}{$start} = join("\t",$end,$ref,$indel,$cov,$indelfreq,$thousand,$dbSNP,$func,$gene,$exfunc,$AAchange);	# info
				++$it_FF;
			}
		}
		close(FF_FILE);

		open FFPE_FILE, $file[1] or die "no $file[1]";
		while(<FFPE_FILE>){
			chomp();
			my ( $chr, $start, $end, $ref, $indel, $cov, $indelfreq, $func, $gene, $exfunc, $AAchange, $thousand, $dbSNP ) = qw/0 0 0 0 0 0 0 0 0 0 0 0 0/;
			my @l = split(/\t/,$_);			# array @l w/o tab
			# col00 Func,Gene,ExonicFunc,AAChange,
			# col04 Conserved,SegDup,ESP6500_ALL,
			# col07 1000g2012feb_ALL,dbSNP135,
			# col09 OTHERS: AVSIFT,LJB_PhyloP,LJB_PhyloP_Pred,LJB_SIFT,LJB_SIFT_Pred,LJB_PolyPhen2,LJB_PolyPhen2_Pred,
			#               LJB_LRT,LJB_LRT_Pred,LJB_MutationTaster,LJB_MutationTaster_Pred,LJB_GERP++,
			# col21 Chr,Start,End,Ref,Obs,Zygosity,snpQ_sum,coverage_per_sample
			# col29 mapQ,snpQ,varFreq,varFreq_mean,varCount
			if ($l[0] eq "Func") {next;}
			$chr = $l[21];
			$start = $l[22];
			$end = $l[23];
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
			
			if ( exists $bed{$chr}{$start} ) {
				$FFPE_indel{$chr}{$start} = $indel;
				$w_FFPE{$chr}{$start} = join("\t",$ref,$indel,$cov,$indelfreq);	# info
				$w_FFPEHelp{$chr}{$start} = join("\t",$thousand,$dbSNP,$func,$gene,$exfunc,$AAchange,$ref,$indel,$cov,$indelfreq);
				++$it_FFPE;
			}
		}
		close(FFPE_FILE);
		
		my ($start, $start_old) = qw/0 0/;
		open PILEUP_FF, $file[4] or die "no $file[4]";
		while(<PILEUP_FF>){
			chomp();
			my @l = split(/\t/,$_);			# array @l w/o tab
			my $chr = $l[0];
			$start = $l[1];
			$start == $start_old and next;
			
			$start_old = $l[1];
			my $covFF = $l[7];
			if ( exists $FFPE_indel{$chr}{$start} && $covFF >= $covThr[$p] ) {
				my $ref = $l[2];
				my $indel = $l[3];
				$compareFF{$chr}{$start} = $indel."\t".$covFF;
			}
			if ( exists $FFPE_indel{$chr}{$start} && $covFF >= $covThr[$p] ) {++$it_interFFPE;}
		}
		close(PILEUP_FF);

		open PILEUP_FFPE, $file[5] or die "no $file[5]";
		while(<PILEUP_FFPE>){
			chomp();
			my @l = split(/\t/,$_);			# array @l w/o tab
			my $chr = $l[0];
			$start = $l[1];
			$start == $start_old and next;
			
			$start_old = $l[1];
			my $covFFPE = $l[7];
			if ( exists $FF_indel{$chr}{$start} && $covFFPE >= $covThr[$p] ) {
				my $ref = $l[2];
				my $indel = $l[3];
				$compareFFPE{$chr}{$start} = $indel."\t".$covFFPE;
			}
			if ( exists $FF_indel{$chr}{$start} && $covFFPE >= $covThr[$p] ) {++$it_interFF;}
		}
		close(PILEUP_FFPE);

		foreach my $chr ( sort keys %FF_indel){  	# $chr: chr \t pos \t ref
			if (defined $FFPE_indel{$chr}){
				foreach my $start (sort keys %{$FF_indel{$chr}}){
					if (defined $FFPE_indel{$chr}{$start}) {
						++$posit;
						if ($FF_indel{$chr}{$start} eq $FFPE_indel{$chr}{$start}){
							$con{$chr}{$start} = join ("\t",$w_FF{$chr}{$start},$w_FFPE{$chr}{$start});
							++$con;
						}
						else {
							$dis{$chr}{$start} = join ("\t",$w_FF{$chr}{$start},$w_FFPE{$chr}{$start});
							++$dis;
						}
					}
				}
			}
		}
		foreach my $chr ( sort keys %FF_indel){  	# $chr: chr \t pos \t ref	
			foreach my $start (sort keys %{$FF_indel{$chr}}){
				if  (!exists $dis{$chr}{$start} && !exists $con{$chr}{$start} && exists $compareFFPE{$chr}{$start} ) {
					$FN{$chr}{$start} = join ("\t",$w_FF{$chr}{$start},"\t".$compareFFPE{$chr}{$start});
					++$FN;
				}
			}
		}

		foreach my $chr ( sort keys %FFPE_indel){  	# $chr: chr \t pos \t ref
			foreach my $start (sort keys %{$FFPE_indel{$chr}}){
				if  (!exists $dis{$chr}{$start} && !exists $con{$chr}{$start} && exists $compareFF{$chr}{$start}) {
					$FP{$chr}{$start} = join ("\t","\t\t".$compareFF{$chr}{$start},"\t".$w_FFPEHelp{$chr}{$start});
					++$FP;
				}
			}
		}

		$sum1 = $con + $dis + $FN;
		$sum2 = $con + $dis + $FP;
		
		$div  = $con/$posit*100;
		$div1 = $con/$sum1*100;
		$div2 = $con/$sum2*100;
		
		print "$sample[$o]\t$method[$p]\t$it_FF\t$it_FFPE\t$it_interFF\t$it_interFFPE\t$posit\t$con\t$dis\t$FN\t$FP\t$div\t$div1\t$div2\n";
			
			
		open WRITE_CON, $file[2] or die;
		print WRITE_CON "chr\tstart\tend\tref\t$sample[$o]_FF\tDepth\tInDel frequency\t1000g2012feb\tdbSNP135\tFunc\tGene\tExonicFunc\tAAChange\tref\t$sample[$o]_FFPE\tDepth\tInDel frequency\n";

		foreach my $chr ( sort keys %con ){
			foreach my $start (sort keys %{$con{$chr}} ){
				print WRITE_CON $chr."\t".$start."\t".$con{$chr}{$start}."\n";
			}
		}
		close (WRITE_CON);

		open WRITE_DIS, $file[3] or die;
		print WRITE_DIS "chr\tstart\tend\tref\t$sample[$o]_FF\tDepth\tInDel frequency\t1000g2012feb\tdbSNP135\tFunc\tGene\tExonicFunc\tAAChange\tref\t$sample[$o]_FFPE\tDepth\tInDel frequency\n";

		foreach my $chr ( sort keys %dis ){
			foreach my $start (sort keys %{$dis{$chr}} ){
				print WRITE_DIS $chr."\t".$start."\t".$dis{$chr}{$start}."\n";
			}
		} print WRITE_DIS "\n";
		foreach my $chr ( sort keys %FN ){
			foreach my $start (sort keys %{$FN{$chr}} ){
				print WRITE_DIS $chr."\t".$start."\t".$FN{$chr}{$start}."\n";
			}
		} print WRITE_DIS "\n";
		foreach my $chr ( sort keys %FP ){
			foreach my $start (sort keys %{$FP{$chr}} ){
				print WRITE_DIS $chr."\t".$start."\t".$FP{$chr}{$start}."\n";
			}
		}
		close(WRITE_DIS);
	}
}
print "\nDone. \n";
print strftime "%Y-%m-%d %H:%M:%S", localtime(time);
print "\n";

