#!/bin/bash
#for (i in 1:13) {
#print("Set i (1-13)");
SAMPLE=array(0,13)
threshold = "q43cov13";
#threshold = "pass";

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

FILE_FF = sprintf("/project/results/me/WES/plots/alleleFreq/genotyper/density_data/%s_AF2_FF_%s_100",SAMPLE[s],threshold);
FILE_FFPE = sprintf("/project/results/me/WES/plots/alleleFreq/genotyper/density_data/%s_AF2_FFPE_%s_100",SAMPLE[s],threshold);

pdffile =  sprintf("/project/results/me/WES/plots/alleleFreq/genotyper/density_plot/%s_AF2_%s_100_density.pdf",SAMPLE[s],threshold);

AF_FF=read.table(FILE_FF,sep="\t",header=TRUE);
AF_FFPE=read.table(FILE_FFPE,sep="\t",header=TRUE);

x2=c((0:5)*20)
densFF=density(AF_FF$GATK);
densFFPE=density(AF_FFPE$GATK);

pdf (file = pdffile, width=300, height=200, pointsize=300)
par (oma=c(0,0,0,0),mar=c(5,6,1,2),col=textcolor)
plot(densFFPE,type="p",pch=3,lwd=30,xlab=("Allele Frequency of SNVs"),cex.lab=1.8,cex.axis=1.8, main = "",xlim = c(4,96.2),cex=1.2,col=textcolor)
lines(densFF,type="p",pch=20,lwd=30,col="#9ACD32",cex=1.2)
#lines(density(AF_FF$PILEUP),col="blue",lty=4,lwd=2)
#lines(density(AF_FFPE$PILEUP),col="red",lty=5,lwd=2)
box("plot", col=textcolor)
legend("top",c(sprintf("FFPE (%s)",length(AF_FFPE$GATK)),sprintf("FF (%s)",length(AF_FF$GATK))),col=c(textcolor,"#9ACD32"),pch=c(3,20),cex=1.5,pt.cex=c(1.3,1.3),border = textcolor,text.col=textcolor)
dev.off()


say=sprintf("%s created!",pdffile);
print(say);

