#!/bin/bash
#for (i in 1:13) {
#print("Set i (1-13)");
SAMPLE=array(0,13)
#threshold = "q43cov13";
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

FILE_FF = sprintf("/project/results/me/TES/plots/alleleFreq/genotyper/density_data/%s_AF2_FF_%s_100",SAMPLE[i],threshold);
FILE_FFPE = sprintf("/project/results/me/TES/plots/alleleFreq/genotyper/density_data/%s_AF2_FFPE_%s_100",SAMPLE[i],threshold);

pdffile =  sprintf("/project/results/me/TES/plots/alleleFreq/genotyper/density_plot/%s_AF2_%s_100_density.pdf",SAMPLE[i],threshold);

AF_FF=read.table(FILE_FF,sep="\t",header=TRUE);
AF_FFPE=read.table(FILE_FFPE,sep="\t",header=TRUE);

x2=c((0:5)*20)


pdf (file = pdffile, width=300, height=200, pointsize=300)
par (oma=c(0,0,0,0),mar=c(5,6,1,2),col=textcolor)
plot(densFFPE,type="p",pch=3,lwd=30,xlab=("Allele Frequency of SNVs"),cex.lab=1.8,cex.axis=1.8, main = "",xlim = c(4,96.2),cex=1.2,col.lab=textcolor,col.axis=textcolor)
lines(densFF,type="p",pch=20,lwd=30,col="#A2CD5A",cex=1.2)
#lines(density(AF_FF$PILEUP),col="blue",lty=4,lwd=2)
#lines(density(AF_FFPE$PILEUP),col="red",lty=5,lwd=2)
legend("top",c(sprintf("FFPE (%s)",length(AF_FFPE$GATK)),sprintf("FF (%s)",length(AF_FF$GATK))),col=c("black","#A2CD5A"),text.col=textcolor,pch=c(3,20),cex=1.5,pt.cex=c(1.3,1.3))
dev.off()

say=sprintf("%s created!",pdffile);
print(say);

pdf (file = pdffile, width=300, height=200, pointsize=300)
par (oma=c(0,0,0,0),mar=c(5,6,1,2),col=textcolor)
hist(AF_FFPE$GATK,breaks=100,freq=TRUE, main = "",ylim=c(0,6000),xlim = c(0,105))
par(new=T)
hist(AF_FF$GATK,breaks=100,freq=TRUE,col="green", main = "",xaxt="n",yaxt="n",ylab="",xlab="",ylim=c(0,6000),xlim = c(0,105))
#par(new=T)
#hist(AF_FF$PILEUP,breaks=100,freq=TRUE,col="blue", main = "",xaxt="n",yaxt="n",ylab="",xlab="",ylim=c(0,6000),xlim = c(0,105))
#par(new=T)
#hist(AF_FFPE$PILEUP,breaks=100,freq=TRUE,col="red", main = "",xaxt="n",yaxt="n",ylab="",xlab="",ylim=c(0,6000),xlim = c(0,105))
plot(densFFPE,type="p",pch=3,lwd=30,xlab=("Allele Frequency of SNVs"),cex.lab=2,cex.axis=2, main = "",xlim = c(2,96.2),cex=1.2,col.lab=textcolor,col.axis=textcolor)
lines(densFF,type="p",pch=20,lwd=30,col="#A2CD5A",cex=1.2)
#lines(density(AF_FF$PILEUP),col="blue",lty=4,lwd=2)
#lines(density(AF_FFPE$PILEUP),col="red",lty=5,lwd=2)
legend("top",c(sprintf("FFPE (%s)",length(AF_FFPE$GATK)),sprintf("FF (%s)",length(AF_FF$GATK))),col=c("black","#A2CD5A"),text.col=textcolor,pch=c(3,20),cex=1.5,pt.cex=c(1.2,1.5))
dev.off()



x2=c((0:5)*20)
#densFF=density(AF_FF$PILEUP);
#densFFPE=density(AF_FFPE$PILEUP);
densFF=density(AF_FF$GATK);
densFFPE=density(AF_FFPE$GATK);
svg(file = snvfile, width=150, height=100, pointsize=100)
par (oma=c(0,0,0,0),mar=c(5,6,1,2),col=textcolor)
plot(densFFPE,type="b",lty=1,pch=3,lwd=20,xlab=sprintf("Allele Frequency of SNVs, n=%s",length(AF_FF$GATK)),cex.lab=2,cex.axis=2, main = "",axes=FALSE,xlim = c(0,96.2),cex=2,col.lab=textcolor,col.axis=textcolor)
axis(1,lty=1,lwd=15,cex.axis=1.5,col.axis=textcolor)
axis(2,lty=1,lwd=15,cex.axis=1.5,col.axis=textcolor)
lines(densFF,type="b",lty=1,pch=20,lwd=20,col="#A2CD5A",cex=2)
legend("top",c("FFPE","FF"),col=c("black","#A2CD5A"),lty=c(1,1),pch=c(3,20),lwd=15,cex=3,pt.cex=c(3,3),bty="n")
dev.off()

pdf (file = pdffile, width=300, height=200, pointsize=300)
par (oma=c(0,0,0,0),mar=c(5,6,1,2))
hist(AF_FFPE$GATK,breaks=100,freq=TRUE, main = "",ylim=c(0,6000),xlim = c(0,120))
par(new=T)
hist(AF_FF$GATK,breaks=100,freq=TRUE,col="green", main = "",xaxt="n",yaxt="n",ylab="",xlab="",ylim=c(0,6000),xlim = c(0,120))

plot(densFFPE,type="b",lty=1,pch=3,lwd=30,xlab=sprintf("Allele Frequency of SNVs, n=%s","6105"),cex.lab=2,cex.axis=2, main = "",xlim = c(0,96.2),cex=1)
#par(new=T)
#plot(densFFPE,type="l",col="#A2CD5A",lty="dashed",lwd=2,xaxt="n",yaxt="n",ylab="",xlab="",cex.lab=2,cex.axis=2, main = "",xlim = c(0,200),ylim=c(0,0.08))
lines(densFF,type="b",lty=1,pch=20,lwd=30,col="#A2CD5A",cex=1)
#lines(density(AF_FF$GATK),col="blue",lty=4,lwd=2)
#lines(density(AF_FFPE$GATK),col="red",lty=5,lwd=2)
legend("top",c("FFPE","FF"),col=c("black","#A2CD5A"),lty=c(1,1),pch=c(3,20),lwd=c(30,30),text.col=textcolor,cex=1.5,pt.cex=c(1,1))
dev.off()

say=sprintf("%s created!",pdffile);
print(say);

#par (mfrow=c(1,3),oma=c(4,2,2,1),mar=c(4,2,3,3))
#par(new=T)
#plot(density(AF_FF$countPILEUP),col="#A2CD5A",xaxt="n",yaxt="n",ylab="",xlab="", main = "",xlim=c(0,100))
#par(new=T)
#plot(density(AF_FFPE$countPILEUP),xaxt="n",yaxt="n",ylab="",xlab="", main = "",xlim=c(0,100),lty=3)
#hist(AF_FF$countGATK,col="#A2CD5A",xlab="Allele Frequency of SNVs",xlim=c(0,100),cex.lab=2,cex.axis=2, main = "")
#hist(AF_FFPE$countGATK,xaxt="n",yaxt="n",ylab="",xlab="", main = "",xlim=c(0,100))

#densFF=density(AF_FF$GATK);
#densFFPE=density(AF_FFPE$GATK);

for (j in 1:length(AF_FF[,1])){
if (AF_FF[j,1]>100) {say=sprintf("%s is groesser!%s",AF_FF[j,1],j);print(say);}
}

