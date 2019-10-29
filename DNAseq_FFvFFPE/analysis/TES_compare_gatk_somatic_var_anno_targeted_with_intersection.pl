#!/usr/bin/perl -w
# Author Yen Hoang, M.Sc.
# 2013

use strict;
use POSIX qw/strftime/;

# Compares somatic SNVs between FF and FFPE in targeted and interesected positions.

print strftime "%Y-%m-%d %H:%M:%S", localtime(time);
print "\n";
print STDERR "./TES_compare_gatk_somatic_var_anno_targeted_with_intersection.pl \n";

my @sample = ("C","D");
my @sample2 = ("FF","FFPE");

our @method = ( "q43cov13");

our ($NORMAL, $TUMOR, $INTERSECT, $WRITE_CON, $WRITE_DIS, $TITLE);
$TITLE = "chr\tpos\tref\tvar_NORMAL\thet\tDepth\tAllele freq\t1000g\tdbSNP\tFunc\tGene\tExonicFunc\tAAChange\tLJB_SIFT_Pred\tref\tvar_TUMOR\thet\tDepth\tAllele freq\n";

print "sample\tthreshold\ttotal_NORMAL\ttotal_TUMOR\tTs/Tv_NORMAL\tTs/Tv_TUMOR\tintersect_NORMAL\tintersect_TUMOR\tTs/Tv_NORMAL\tTs/Tv_TUMOR\tsame positions\tcon\tdis\tFP\t% con\n";

my $targeted = "</project/results/me/TES/targeted_positions.list";
my %bed;
my $it_bed = 0;

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

for (my $o = 0; $o <= $#sample; $o++) {
	for (my $q = 0; $q <= $#sample2; $q++) {
		
		foreach my $p (@method) {
		
			$NORMAL = "</project/results/me/TES/$sample[$o]N_$sample2[$q]/genotyper/anno/$sample[$o]N_$sample2[$q].bwa.ug.snp_${p}.genome_summary.csv";
			$TUMOR = "</project/results/me/TES/$sample[$o]T_$sample2[$q]/genotyper/anno/$sample[$o]T_$sample2[$q].bwa.ug.snp_${p}.genome_summary.csv";
			$INTERSECT = "</project/results/me/TES/pileup_distribution/$sample[$o]N_$sample[$o]T_$sample2[$q]_intersect.distribution";
			$WRITE_CON = ">/project/results/me/TES/SNV_comparison/genotyper/with_intersection/somatic/$sample[$o]N_$sample[$o]T_$sample2[$q]_${p}_exonic_con";
			$WRITE_DIS = ">/project/results/me/TES/SNV_comparison/genotyper/with_intersection/somatic/$sample[$o]N_$sample[$o]T_$sample2[$q]_${p}_exonic_dis";
	
			my ( $it_NORMAL, $it_interNORMAL, $it_TUMOR, $it_interTUMOR, $con, $dis, $FN, $FP, $posit) = qw/0 0 0 0 0 0 0 0 0 0 0/; 
			my (%NORMAL, %NORMAL_snp, %NORMAL_ref, %w_NORMAL, %TUMOR, %TUMOR_snp, %TUMOR_ref, %w_TUMOR, %w_TUMORHelp, %con, %dis, %FN, %FP);

			open NORMAL_FILE, $NORMAL or die "no $NORMAL";
			while(<NORMAL_FILE>){
				chomp();
				my ( $chr, $pos, $ref, $snp, $het, $cov, $snpfreq, $func, $gene, $exfunc, $AAchange, $sift, $thousand, $dbSNP ) = qw/0 0 0 0 0 0 0 0 0 0 0 0 0/;
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
				if (exists $bed{$chr}{$pos}) {
					$NORMAL_snp{$chr}{$pos} = $snp;
					$NORMAL_ref{$chr}{$pos} = $ref;
					$NORMAL{$chr}{$pos} = $snp."\t".$het;
					$w_NORMAL{$chr}{$pos} = join("\t",$ref,$snp,$het,$cov,$snpfreq,$thousand,$dbSNP,$func,$gene,$exfunc,$AAchange,$sift);	# info
					++$it_NORMAL;
				}
			}
			close(NORMAL_FILE);

			open TUMOR_FILE, $TUMOR or die "no $TUMOR";
			while(<TUMOR_FILE>){
				chomp();
				my ( $chr, $pos, $ref, $snp, $het, $cov, $snpfreq, $func, $gene, $exfunc, $AAchange, $sift, $thousand, $dbSNP ) = qw/0 0 0 0 0 0 0 0 0 0 0 0 0 0/;
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
				
				if (exists $bed{$chr}{$pos}) {
					$TUMOR_ref{$chr}{$pos} = $ref;
					$TUMOR_snp{$chr}{$pos} = $snp;
					$TUMOR{$chr}{$pos} = $snp."\t".$het;
					$w_TUMOR{$chr}{$pos} = join("\t",$ref,$snp,$het,$cov,$snpfreq);	# info
					$w_TUMORHelp{$chr}{$pos} = join("\t",$thousand,$dbSNP,$func,$gene,$exfunc,$AAchange,$sift,$ref,$snp,$het,$cov,$snpfreq);
					++$it_TUMOR;
				}
			}
			close(TUMOR_FILE);

			my (%compareNORMAL, %compareTUMOR);
			open PILEUP, $INTERSECT or die "no $INTERSECT";
			while(<PILEUP>){
				chomp();
				my @l = split(/\t/,$_);			# array @l w/o tab
				$l[0] eq "chr" and next;
				my $chr = $l[0];
				my $pos = $l[1];
				my $covNORMAL = $l[3] + $l[4] + $l[5] + $l[6] + $l[7] + $l[8] + $l[9] + $l[10];
				my $covTUMOR = $l[12] + $l[13] + $l[14] + $l[15] + $l[16] + $l[17] + $l[18] + $l[19];
				
				if ( exists $NORMAL_snp{$chr}{$pos} ) {
					$compareTUMOR{$chr}{$pos} = $covTUMOR;
					++$it_interNORMAL;
				} 
		
				if ( exists $TUMOR_snp{$chr}{$pos} ) {
					$compareNORMAL{$chr}{$pos} = $covNORMAL;
					++$it_interTUMOR;
				} 
			}
			close(PILEUP);

			my ( $Ts_NORMAL, $Ts_TUMOR, $Tv_NORMAL, $Tv_TUMOR, $Ts2_NORMAL, $Ts2_TUMOR, $Tv2_NORMAL, $Tv2_TUMOR, $Ts_Tv_ratio_NORMAL, $Ts_Tv_ratio2_NORMAL, $Ts_Tv_ratio_TUMOR, $Ts_Tv_ratio2_TUMOR) = qw/0 0 0 0 0 0 0 0 0 0 0 0/;

			foreach my $chr ( sort keys %NORMAL_snp){  	# $chr: chr \t pos \t ref	
				if (defined $TUMOR_snp{$chr}){
					foreach my $pos (sort keys %{$NORMAL_snp{$chr}}){
						if (defined $TUMOR_snp{$chr}{$pos}) {
							++$posit;
							if ($NORMAL{$chr}{$pos} eq $TUMOR{$chr}{$pos}){
								$con{$chr}{$pos} = join ("\t",$w_NORMAL{$chr}{$pos},$w_TUMOR{$chr}{$pos});
								++$con;
							}
							else {
								$dis{$chr}{$pos} = join ("\t",$w_NORMAL{$chr}{$pos},$w_TUMOR{$chr}{$pos});
								++$dis;
							}
						
							my $done = 0;
							if ( $NORMAL_ref{$chr}{$pos} eq "A" || $NORMAL_ref{$chr}{$pos} eq "G") {
								if ( $NORMAL_snp{$chr}{$pos} eq "A" || $NORMAL_snp{$chr}{$pos} eq "G") { ++$Ts2_NORMAL; $done = 1;}
							} elsif ( $NORMAL_snp{$chr}{$pos} eq "C" || $NORMAL_snp{$chr}{$pos} eq "T") { ++$Ts2_NORMAL; $done = 1;}
			
							if ( $done != 1) { ++$Tv2_NORMAL; }
				
							$done = 0;
							if ( $TUMOR_ref{$chr}{$pos} eq "A" || $TUMOR_ref{$chr}{$pos} eq "G") {
								if ( $TUMOR_snp{$chr}{$pos} eq "A" || $TUMOR_snp{$chr}{$pos} eq "G") { ++$Ts2_TUMOR; $done = 1;}
							} elsif ( $TUMOR_snp{$chr}{$pos} eq "C" || $TUMOR_snp{$chr}{$pos} eq "T") { ++$Ts2_TUMOR; $done = 1;}
				
							if ( $done != 1) { ++$Tv2_TUMOR; }
						}
					}
				}
			}

			foreach my $chr ( sort keys %TUMOR_snp){  	# $chr: chr \t pos \t ref
				foreach my $pos (sort keys %{$TUMOR_snp{$chr}}){
					my $done = 0;
					if ( $TUMOR_ref{$chr}{$pos} eq "A" || $TUMOR_ref{$chr}{$pos} eq "G") {
						if ( $TUMOR_snp{$chr}{$pos} eq "A" || $TUMOR_snp{$chr}{$pos} eq "G") { ++$Ts_TUMOR; $done = 1;}
					} elsif ( $TUMOR_snp{$chr}{$pos} eq "C" || $TUMOR_snp{$chr}{$pos} eq "T") { ++$Ts_TUMOR; $done = 1;}
				
					if ( $done != 1) { ++$Tv_TUMOR; }
				
					if  (!exists $dis{$chr}{$pos} && !exists $con{$chr}{$pos} && exists $compareNORMAL{$chr}{$pos}) {
						if ( ($p eq "both" || $p eq "q43cov13") && $compareNORMAL{$chr}{$pos} > 12) {
							$FP{$chr}{$pos} = join ("\t",$TUMOR_ref{$chr}{$pos},"\t\t".$compareNORMAL{$chr}{$pos},"\t".$w_TUMORHelp{$chr}{$pos});
							++$FP;
							if ( $done != 1) { ++$Tv2_TUMOR; } else { ++$Ts2_TUMOR; }
						} elsif ($p eq "pass" && $compareNORMAL{$chr}{$pos} > 1) {
							$FP{$chr}{$pos} = join ("\t",$TUMOR_ref{$chr}{$pos},"\t\t".$compareNORMAL{$chr}{$pos},"\t".$w_TUMORHelp{$chr}{$pos});
							++$FP;
							if ( $done != 1) { ++$Tv2_TUMOR; } else { ++$Ts2_TUMOR; }
						}
					}
				}
			}
			
			foreach my $chr ( sort keys %NOMAL_snp){  	# $chr: chr \t pos \t ref
				foreach my $pos (sort keys %{$NOMAL_snp{$chr}}){
					my $done = 0;
					if ( $NOMAL_ref{$chr}{$pos} eq "A" || $NOMAL_ref{$chr}{$pos} eq "G") {
						if ( $NOMAL_snp{$chr}{$pos} eq "A" || $NOMAL_snp{$chr}{$pos} eq "G") { ++$Ts_NOMAL; $done = 1;}
					} elsif ( $NOMAL_snp{$chr}{$pos} eq "C" || $NOMAL_snp{$chr}{$pos} eq "T") { ++$Ts_NOMAL; $done = 1;}
				
					if ( $done != 1) { ++$Tv_NOMAL; }
				}
			}
				
			if ( $Tv_NORMAL != 0 ) { $Ts_Tv_ratio_NORMAL = $Ts_NORMAL/$Tv_NORMAL; }
			if ( $Tv2_NORMAL != 0 ) { $Ts_Tv_ratio2_NORMAL = $Ts2_NORMAL/$Tv2_NORMAL; }
			if ( $Tv_TUMOR != 0 ) { $Ts_Tv_ratio_TUMOR = $Ts_TUMOR/$Tv_TUMOR; }
			if ( $Tv2_TUMOR != 0 ) { $Ts_Tv_ratio2_TUMOR = $Ts2_TUMOR/$Tv2_TUMOR; }
			
			my $div = $con/$posit*100;
			
			print "$sample[$o]($sample2[$q])\t${p}\t$it_NORMAL\t$it_TUMOR\t$Ts_Tv_ratio_NORMAL\t$Ts_Tv_ratio_TUMOR\t$it_interNORMAL\t$it_interTUMOR\t$Ts_Tv_ratio2_NORMAL\t$Ts_Tv_ratio2_TUMOR\t$posit\t$con\t$dis\t$FP\t$div\n";
				

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
			foreach my $chr ( sort keys %FP ){
				foreach my $pos (sort keys %{$FP{$chr}} ){
					print WRITE_DIS $chr."\t".$pos."\t".$FP{$chr}{$pos}."\n";
				}
			}
			close(WRITE_DIS);
		}
	}
}
print "Done. \n";
print strftime "%Y-%m-%d %H:%M:%S", localtime(time);
print "\n";
