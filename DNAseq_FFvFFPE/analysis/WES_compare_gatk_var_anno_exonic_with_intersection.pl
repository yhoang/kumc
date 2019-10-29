#!/usr/bin/perl -w
# Author Yen Hoang, M.Sc.
# 2013

use strict;
use POSIX qw/strftime/;

# Compares annotated SNVs between FF and FFPE in exonic intersected positions.

print strftime "%Y-%m-%d %H:%M:%S", localtime(time);
print "\n";
print STDERR "./WES_compare_gatk_var_anno_exonic_with_intersection.pl \n";

our @method = ("q43cov13");
our @sample = ("2474","2561","2640","2685","2938","3050","3356","4079","4191","CN","CT","DN","DT");
our ($FF, $FFPE, $INTERSECT, $WRITE_CON, $WRITE_DIS, $TITLE);

$TITLE = "chr\tpos\tref\tFF_var\thet\tDepth\tAllele Freq\t1000G\tdbSNP\tFunc\tGene\tExonicFunc\tAAChange\tLJB_SIFT_Pred\tFFPE_var\thet\tDepth\tAllele Freq\n";

print "sample\tthreshold\ttotal_FF\ttotal_FFPE\tTs/Tv_FF\tTs/Tv_FFPE\tintersect_FF\tintersect_FFPE\tTs/Tv_FF\tTs/Tv_FFPE\tsame positions\tcon\tdis\tFN\tFP\t% con\t% con_FF\t% con_FFPE\n";

for (my $o = 0; $o <= $#sample; $o++) {

	foreach my $p (@method) {

		$FF = "</project/results/me/WES/$sample[$o]_FFex/genotyper/anno/$sample[$o]_FFex.bwa.ug.snp_${p}.exome_summary.csv";
		$FFPE = "</project/results/me/WES/$sample[$o]_FFPEex/genotyper/anno/$sample[$o]_FFPEex.bwa.ug.snp_${p}.exome_summary.csv";
		$WRITE_CON = ">/project/results/me/WES/SNV_comparison/genotyper/with_intersection/base_call/$sample[$o]_FFex_FFPEex_${p}_exonic_con";
		$WRITE_DIS = ">/project/results/me/WES/SNV_comparison/genotyper/with_intersection/base_call/$sample[$o]_FFex_FFPEex_${p}_exonic_dis";
		$INTERSECT = "</project/results/me/WES/pileup_distribution/$sample[$o]_FFex_FFPEex_intersect.distribution";
		
		my ( $it_FF, $it_interFF, $it_FFPE, $it_interFFPE, $con, $dis, $FN, $FP, $posit) = qw/0 0 0 0 0 0 0 0 0 0 0/; 
		my (%FF, %FF_ref, %w_FF, %FFPE, %FFPE_ref, %w_FFPE, %w_FFPEHelp, %con, %dis, %FN, %FP);

		open FF_FILE, $FF or die "no $FF";
		while(<FF_FILE>){
			chomp();
			my ( $chr, $pos, $end, $ref, $snp, $het, $cov, $snpfreq, $func, $gene, $exfunc, $AAchange, $sift, $thousand, $dbSNP ) = qw/0 0 0 0 0 0 0 0 0 0 0 0 0/;
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
			$pos = $l[22];
			$end = $l[22];
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
		
			$FF{$chr}{$pos} = $snp;
			$FF_ref{$chr}{$pos} = $ref;
			$w_FF{$chr}{$pos} = join("\t",$ref,$snp,$het,$cov,$snpfreq,$thousand,$dbSNP,$func,$gene,$exfunc,$AAchange,$sift);	# info
			++$it_FF;
		}
		close(FF_FILE);

		open FFPE_FILE, $FFPE or die "no $FFPE";
		while(<FFPE_FILE>){
			chomp();
			my ( $chr, $pos, $end, $ref, $snp, $het, $cov, $snpfreq, $func, $gene, $exfunc, $AAchange, $sift, $thousand, $dbSNP ) = qw/0 0 0 0 0 0 0 0 0 0 0 0 0 0/;
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
			$pos = $l[22];
			$end = $l[22];
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
	
			$FFPE_ref{$chr}{$pos} = $ref;
			$FFPE{$chr}{$pos} = $snp;
			$w_FFPE{$chr}{$pos} = join("\t",$snp,$het,$cov,$snpfreq);	# info
			$w_FFPEHelp{$chr}{$pos} = join("\t",$thousand,$dbSNP,$func,$gene,$exfunc,$AAchange,$sift,$snp,$het,$cov,$snpfreq);
			++$it_FFPE;
		}
		close(FFPE_FILE);

		my (%compareFF, %compareFFPE);
		open PILEUP, $INTERSECT or die "no $INTERSECT";
		while(<PILEUP>){
			chomp();
			my @l = split(/\t/,$_);			# array @l w/o tab
			$l[0] eq "chr" and next;
			my $chr = $l[0];
			my $pos = $l[1];
			my $covFF = $l[3] + $l[4] + $l[5] + $l[6] + $l[7] + $l[8] + $l[9] + $l[10];
			my $covFFPE = $l[12] + $l[13] + $l[14] + $l[15] + $l[16] + $l[17] + $l[18] + $l[19];
				
			if ( exists $FF{$chr}{$pos} ) {
				$compareFFPE{$chr}{$pos} = $covFFPE;
				++$it_interFF;
			} 
		
			if ( exists $FFPE{$chr}{$pos} ) {
				$compareFF{$chr}{$pos} = $covFF;
				++$it_interFFPE;
			} 
		}
		close(PILEUP);

		my ( $Ts_FF, $Ts_FFPE, $Tv_FF, $Tv_FFPE, $Ts2_FF, $Ts2_FFPE, $Tv2_FF, $Tv2_FFPE, $Ts_Tv_ratio_FF, $Ts_Tv_ratio2_FF, $Ts_Tv_ratio_FFPE, $Ts_Tv_ratio2_FFPE) = qw/0 0 0 0 0 0 0 0 0 0 0 0/;

		foreach my $chr ( sort keys %FF){  	# $chr: chr \t pos \t ref	
			if (defined $FFPE{$chr}){
				foreach my $pos (sort keys %{$FF{$chr}}){
					if (defined $FFPE{$chr}{$pos}) {
						++$posit;
						if ($FF{$chr}{$pos} eq $FFPE{$chr}{$pos}){
							$con{$chr}{$pos} = join ("\t",$w_FF{$chr}{$pos},$w_FFPE{$chr}{$pos});
							++$con;
						}
						else {
							$dis{$chr}{$pos} = join ("\t",$w_FF{$chr}{$pos},$w_FFPE{$chr}{$pos});
							++$dis;
						}
						
						my $done = 0;
						if ( $FF_ref{$chr}{$pos} eq "A" || $FF_ref{$chr}{$pos} eq "G") {
							if ( $FF{$chr}{$pos} eq "A" || $FF{$chr}{$pos} eq "G") { ++$Ts2_FF; $done = 1;}
						} elsif ( $FF{$chr}{$pos} eq "C" || $FF{$chr}{$pos} eq "T") { ++$Ts2_FF; $done = 1;}
			
						if ( $done != 1) { ++$Tv2_FF; }
				
						$done = 0;
						if ( $FFPE_ref{$chr}{$pos} eq "A" || $FFPE_ref{$chr}{$pos} eq "G") {
							if ( $FFPE{$chr}{$pos} eq "A" || $FFPE{$chr}{$pos} eq "G") { ++$Ts2_FFPE; $done = 1;}
						} elsif ( $FFPE{$chr}{$pos} eq "C" || $FFPE{$chr}{$pos} eq "T") { ++$Ts2_FFPE; $done = 1;}
				
						if ( $done != 1) { ++$Tv2_FFPE; }
					}
				}
			}
		}
		foreach my $chr ( sort keys %FF){  	# $chr: chr \t pos \t ref	
			foreach my $pos (sort keys %{$FF{$chr}}){
				my $done = 0;
				if ( $FF_ref{$chr}{$pos} eq "A" || $FF_ref{$chr}{$pos} eq "G") {
					if ( $FF{$chr}{$pos} eq "A" || $FF{$chr}{$pos} eq "G") { ++$Ts_FF; $done = 1;}
				} elsif ( $FF{$chr}{$pos} eq "C" || $FF{$chr}{$pos} eq "T") { ++$Ts_FF; $done = 1;}
			
				if ( $done != 1) { ++$Tv_FF; }
					
				if  (!exists $dis{$chr}{$pos} && !exists $con{$chr}{$pos} && exists $compareFFPE{$chr}{$pos} ) {
				
					if (($p eq "q43cov13" || $p eq "both") && $compareFFPE{$chr}{$pos} > 12) {
						$FN{$chr}{$pos} = join ("\t",$w_FF{$chr}{$pos},"\t",$compareFFPE{$chr}{$pos});
						++$FN;
						if ( $done != 1) { ++$Tv2_FF; } else { ++$Ts2_FF; }
					} elsif ($p eq "pass" && $compareFFPE{$chr}{$pos} > 1) {
						$FN{$chr}{$pos} = join ("\t",$w_FF{$chr}{$pos},"\t",$compareFFPE{$chr}{$pos});
						++$FN;
						if ( $done != 1) { ++$Tv2_FF; } else { ++$Ts2_FF; }
					}
				}
			}
		}

		foreach my $chr ( sort keys %FFPE){  	# $chr: chr \t pos \t ref
			foreach my $pos (sort keys %{$FFPE{$chr}}){
				my $done = 0;
				if ( $FFPE_ref{$chr}{$pos} eq "A" || $FFPE_ref{$chr}{$pos} eq "G") {
					if ( $FFPE{$chr}{$pos} eq "A" || $FFPE{$chr}{$pos} eq "G") { ++$Ts_FFPE; $done = 1;}
				} elsif ( $FFPE{$chr}{$pos} eq "C" || $FFPE{$chr}{$pos} eq "T") { ++$Ts_FFPE; $done = 1;}
				
				if ( $done != 1) { ++$Tv_FFPE; }
				
				if  (!exists $dis{$chr}{$pos} && !exists $con{$chr}{$pos} && exists $compareFF{$chr}{$pos}) {
					if ( ($p eq "q43cov13" || $p eq "both") && $compareFF{$chr}{$pos} > 12) {
						$FP{$chr}{$pos} = join ("\t",$FFPE_ref{$chr}{$pos},"\t\t".$compareFF{$chr}{$pos}."\t",$w_FFPEHelp{$chr}{$pos});
						++$FP;
						if ( $done != 1) { ++$Tv2_FFPE; } else { ++$Ts2_FFPE; }
					} elsif ($p eq "pass" && $compareFF{$chr}{$pos} > 1) {
						$FP{$chr}{$pos} = join ("\t",$FFPE_ref{$chr}{$pos},"\t\t".$compareFF{$chr}{$pos}."\t",$w_FFPEHelp{$chr}{$pos});
						++$FP;
						if ( $done != 1) { ++$Tv2_FFPE; } else { ++$Ts2_FFPE; }
					}
				}
			}
		}
		

		my $sum1 = $con + $dis + $FN;
		my $sum2 = $con + $dis + $FP;
		
		my $div  = $con/$posit*100;
		my $div1 = $con/$sum1*100;
		my $div2 = $con/$sum2*100;
		
		if ( $Tv_FF != 0 ) { $Ts_Tv_ratio_FF = $Ts_FF/$Tv_FF; }
		if ( $Tv2_FF != 0 ) { $Ts_Tv_ratio2_FF = $Ts2_FF/$Tv2_FF; }
		if ( $Tv_FFPE != 0 ) { $Ts_Tv_ratio_FFPE = $Ts_FFPE/$Tv_FFPE; }
		if ( $Tv2_FFPE != 0 ) { $Ts_Tv_ratio2_FFPE = $Ts2_FFPE/$Tv2_FFPE; }
				
		print "$sample[$o]\t${p}\t$it_FF\t$it_FFPE\t$Ts_Tv_ratio_FF\t$Ts_Tv_ratio_FFPE\t$it_interFF\t$it_interFFPE\t$Ts_Tv_ratio2_FF\t$Ts_Tv_ratio2_FFPE\t$posit\t$con\t$dis\t$FN\t$FP\t$div\t$div1\t$div2\n";
				

		open WRITE_CON, $WRITE_CON or die;
		print WRITE_CON $TITLE;

		foreach my $chr ( sort keys %con ){
			foreach my $pos (sort keys %{$con{$chr}} ){
				print WRITE_CON $chr."\t".$pos."\t".$con{$chr}{$pos}."\n";
			}
		}
		close (WRITE_CON);

		open WRITE_DIS, $WRITE_DIS or die;
		print WRITE_DIS "$TITLE";

		foreach my $chr ( sort keys %dis ){
			foreach my $pos (sort keys %{$dis{$chr}} ){
				print WRITE_DIS $chr."\t".$pos."\t".$dis{$chr}{$pos}."\n";
			}
		} print WRITE_DIS "\n";
		foreach my $chr ( sort keys %FN ){
			foreach my $pos (sort keys %{$FN{$chr}} ){
				print WRITE_DIS $chr."\t".$pos."\t".$FN{$chr}{$pos}."\n";
			}
		} print WRITE_DIS "\n";
		foreach my $chr ( sort keys %FP ){
			foreach my $pos (sort keys %{$FP{$chr}} ){
				print WRITE_DIS $chr."\t".$pos."\t".$FP{$chr}{$pos}."\n";
			}
		}
		close(WRITE_DIS);
	}
}

print "Done. \n";
print strftime "%Y-%m-%d %H:%M:%S", localtime(time);
print "\n";

