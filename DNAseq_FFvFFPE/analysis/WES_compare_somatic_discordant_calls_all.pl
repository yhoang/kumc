#!/usr/bin/perl -w
# Author Yen Hoang, M.Sc.
# 2013

use strict;
use POSIX qw/strftime/;

# Compares somatic exonic discordant SNV calls between FF and FFPE.

print strftime "%Y-%m-%d %H:%M:%S", localtime(time);
print "\n";
print STDERR "./compare_gatk_var_anno_with_intersection.pl \n";

my @sample = ("C","D");
our @method = ("q43cov13");
our ($FF, $FFPE, $INTERSECT, $WRITE_CON, $WRITE_DIS, $TITLE);

print "sample\tthreshold\ttotal_FF\ttotal_FFPE\tsame positions\tcon\tdis\n";
$TITLE = "chr\tpos\tref\tvar_FF\thet\mapQ\tcov\tallele freq\t1000g\tdbsnp\tfunc\tgene\tsynonymous\tAAchange\tsift\t\tref\tvar_FF\thet\mapQ\tcov\tallele freq\n";
for (my $o = 0; $o <= $#sample; $o++) {
		
		foreach my $p (@method) {
			$FF = "</project/results/me/WES/SNV_comparison/genotyper/with_intersection/somatic/$sample[$o]N_$sample[$o]T_FFex_${p}_exonic_dis";
			$FFPE = "</project/results/me/WES/SNV_comparison/genotyper/with_intersection/somatic/$sample[$o]N_$sample[$o]T_FFPEex_${p}_exonic_dis";
			$WRITE_DIS = ">/project/results/me/WES/SNV_comparison/genotyper/with_intersection/somatic/somatic_$sample[$o]_FFex_FFPEex_exonic_dis";
			
			my ( $dis_FF, $FN_FF, $FP_FF, $dis_FFPE, $FN_FFPE, $FP_FFPE, $con, $dis, $FN, $FP, $posit, $it_FF, $it_FFPE) = qw/0 0 0 0 0 0 0 0 0 0 0 0 0/; 
			my (%FF, %FF_snp, %FF_ref, %w_FF, %FFPE, %FFPE_snp, %FFPE_ref, %w_FFPE, %w_FFPEHelp, %con, %dis, %FN, %FP);

			open FF_FILE, $FF or die "no $FF";
			while(<FF_FILE>){
				chomp();
				my ( $chr, $pos, $ref, $snp, $het, $mapq, $cov, $snpfreq, $func, $gene, $exfunc, $AAchange, $sift, $thousand, $dbSNP ) = qw/0 0 0 0 0 0 0 0 0 0 0 0 0 0/;
				my @l = split(/\t/,$_);			# array @l w/o tab
				
				if ($l[0] eq "chr" || $l[0] eq "") {next;}
				$chr = $l[0];
				$pos = $l[1];
				$ref = $l[2];
				
				$snp = $l[16];
				$het = $l[17];
				$mapq = $l[18];
				$cov = $l[19];
				$snpfreq = $l[20];
				
				$thousand = $l[8];
				$dbSNP = $l[9];
				$func = $l[10];
				$gene = $l[11];
				$exfunc = $l[12];
				$AAchange = $l[13];
				$sift = $l[14];
		
				$FF_snp{$chr}{$pos} = $snp;
				$FF_ref{$chr}{$pos} = $ref;
				$FF{$chr}{$pos} = $snp."\t".$het;
				$w_FF{$chr}{$pos} = join("\t",$ref,$snp,$het,$mapq,$cov,$snpfreq,$thousand,$dbSNP,$func,$gene,$exfunc,$AAchange,$sift);	# info
				++$it_FF;
			}
			close(FF_FILE);

			open FFPE_FILE, $FFPE or die "no $FFPE";
			while(<FFPE_FILE>){
				chomp();
				my ( $chr, $pos, $ref, $snp, $het, $mapq, $cov, $snpfreq, $func, $gene, $exfunc, $AAchange, $sift, $thousand, $dbSNP ) = qw/0 0 0 0 0 0 0 0 0 0 0 0 0 0/;
				my @l = split(/\t/,$_);			# array @l w/o tab
				
				if ($l[0] eq "chr" || $l[0] eq "") {next;}
				$chr = $l[0];
				$pos = $l[1];
				$ref = $l[2];
				
				$snp = $l[16];
				$het = $l[17];
				$mapq = $l[18];
				$cov = $l[19];
				$snpfreq = $l[20];
				
				$thousand = $l[8];
				$dbSNP = $l[9];
				$func = $l[10];
				$gene = $l[11];
				$exfunc = $l[12];
				$AAchange = $l[13];
				$sift = $l[14];
					
				$FFPE_ref{$chr}{$pos} = $ref;
				$FFPE_snp{$chr}{$pos} = $snp;
				$FFPE{$chr}{$pos} = $snp."\t".$het;
				$w_FFPE{$chr}{$pos} = join("\t",$ref,$snp,$het,$cov,$snpfreq);	# info
				$w_FFPEHelp{$chr}{$pos} = join("\t",$thousand,$dbSNP,$func,$gene,$exfunc,$AAchange,$sift,$ref,$snp,$het,$mapq,$cov,$snpfreq);
				++$it_FFPE;
			}
			close(FFPE_FILE);


			foreach my $chr ( sort keys %FF_snp){  	# $chr: chr \t pos \t ref	
				if (defined $FFPE_snp{$chr}){
					foreach my $pos (sort keys %{$FF_snp{$chr}}){
						if (defined $FFPE_snp{$chr}{$pos}) {
							++$posit;
							if ($FF{$chr}{$pos} eq $FFPE{$chr}{$pos}){
								$con{$chr}{$pos} = join ("\t",$w_FF{$chr}{$pos},$w_FFPE{$chr}{$pos});
								++$con;
							}
							else {
								$dis{$chr}{$pos} = join ("\t",$w_FF{$chr}{$pos},$w_FFPE{$chr}{$pos});
								++$dis;
							}
						}
					}
				}
			}

			print "$sample[$o]\t${p}\t$it_FF\t$it_FFPE\t$posit\t$con\t$dis\n";
			
			open WRITE_DIS, $WRITE_DIS or die "Cannot write int $WRITE_DIS!\n";
			foreach my $chr (sort keys %dis) {
				foreach my $pos (sort keys %{$dis{$chr}}) {
					print "$chr\t$pos\t$dis{$chr}{$pos}\n";
					print WRITE_DIS "$chr\t$pos\t$dis{$chr}{$pos}\n";
				}
			}
	}
}
print "Done. \n";
print strftime "%Y-%m-%d %H:%M:%S", localtime(time);
print "\n";

