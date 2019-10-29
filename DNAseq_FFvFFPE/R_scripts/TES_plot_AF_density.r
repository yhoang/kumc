#!/bin/bash
#for (i in 1:13) {
#print("Set i (1-13)");
SAMPLE=array(0,13)
threshold = "q43cov13";
threshold = "pass";

textcolor = "#2E2E2E";
SAMPLE[1]="2474";
SAMPLE[2]="2561";
SAMPLE[3]="2640";
SAMPLE[4]="2685";
SAMPLE[5]="2938";
SAMPLE[6]="3050";
SAMPLE[7]="3356";
SAMPLE[8]="4079";
SAMPLE[9]="4191";
SAMPLE[10]="CN";
SAMPLE[11]="CT";
SAMPLE[12]="DN";
SAMPLE[13]="DT";

SAMPLE2=array(0,13)
SAMPLE2[10]="14119(N)";
SAMPLE2[11]="14119(T)";
SAMPLE2[12]="22285(N)";
SAMPLE2[13]="22285(T)";

FILE_FF = sprintf("/project/results/me/TES/plots/alleleFreq/genotyper/density_data/%s_AF2_FF_%s_100",SAMPLE[s],threshold);
FILE_FFPE = sprintf("/project/results/me/TES/plots/alleleFreq/genotyper/density_data/%s_AF2_FFPE_%s_100",SAMPLE[s],threshold);

pdffile =  sprintf("/project/results/me/TES/plots/alleleFreq/genotyper/density_plot/%s_AF2_%s_100_density.pdf",SAMPLE[s],threshold);

AF_FF=read.table(FILE_FF,sep="\t",header=TRUE);
AF_FFPE=read.table(FILE_FFPE,sep="\t",header=TRUE);

x2=c((0:5)*20)
densFFPE=density(AF_FFPE$GATK)
densFF=density(AF_FF$GATK)

pdf (file = pdffile, width=300, height=200, pointsize=300)
par (oma=c(0,0,0,0),mar=c(5,6,1,2),col=textcolor)
plot(densFFPE,type="p",pch=3,lwd=30,xlab=("Allele Frequency of SNVs"),cex.lab=2.4,cex.axis=2.2, main = "",xlim = c(4,96.2),cex=1.4,col.lab=textcolor,col.axis=textcolor)
lines(densFF,type="p",pch=20,lwd=30,col="#A2CD5A",cex=1.4)
legend("top",c(sprintf("FFPE (%s)",length(AF_FFPE$GATK)),sprintf("FF (%s)",length(AF_FF$GATK))),col=c("black","#A2CD5A"),text.col=textcolor,pch=c(3,20),cex=1.9,pt.cex=c(1.4,1.4))
dev.off()

say=sprintf("%s created!",pdffile);
print(say);

