#!/bin/bash
file1 = "/project/results/me/WES/plots/mean_cov_PCRdup_barplot/WES_mean_cov_table.list";
file2 = "/project/results/me/WES/plots/mean_cov_PCRdup_barplot/WES_PCRdup.list";
pdffile = "/project/results/me/WES/plots/mean_cov_PCRdup_barplot/WES_mean_cov_PCRdup.pdf";
mean_cov=read.table(file1,sep="\t",header=TRUE);
PCRdup=read.table(file2,sep="\t",header=TRUE);
x2=c((0:10)*10)
x3=c((0:2)*10)
max=50

textcolor = "#2E2E2E";
rownames(mean_cov)=mean_cov[,1]
y1_tick=(0:5)*10
y2_tick=(0:5)*20

pdf (file = pdffile, width=30, height=20, pointsize=30)  
	par (oma=c(0,0,0,0),mar=c(3,6,3,6),col=textcolor)
	barplot(rbind(mean_cov[,2],mean_cov[,3]),beside=TRUE,width=0.008,ylim=c(0,max),col=c("#9ACD32","#1E90FF"),border=NULL,xaxt="n",yaxt="n", main = "",ylab="")
	abline(h=10,col="#A3A3A3",lwd=7)
	abline(h=20,col="#A3A3A3",lwd=7)
	abline(h=30,col="#A3A3A3",lwd=7)
	abline(h=40,col="#A3A3A3",lwd=7)
	barplot(rbind(mean_cov[,2],mean_cov[,3]),beside=TRUE,add=TRUE,width=(0.008),ylim=c(0,max),
	 col=c("#9ACD32","#1E90FF"),cex.lab=2.2, font.lab = 2, main = "",border=c("#9ACD32","#1E90FF"),ylab="Total Mean Coverage",col.lab=textcolor,yaxt="n")
	axis(2,lty=1,lwd=7,at=(y1_tick),labels=y1_tick,las=0, cex.axis=2,col=textcolor,col.lab=textcolor,col.axis=textcolor)
	mtext("     2474    2561  2640   2685   2938  3050   3356   4079   4191   (N)14119(T) (N)22285(T)", side=1,line=1, cex=1.4, col=textcolor,adj=0, font =2)

	par(new=T)
	axis(4,lty=1,lwd=7,at=(y2_tick/2),labels=y2_tick,las=0, cex.axis=2,col=textcolor,col.lab=textcolor,col.axis=textcolor)
	mtext("% PCR Duplicates",side=4,line=3,cex=2.2,col=textcolor,font=2)
	par(new=T)
	plot(c(2.9,10.8,18.7,26.3,34.6,42,50,58,65.9,73.9,81.5,90,97.4),(PCRdup$FF),type="l",lty=1,lwd=12,col="#EE9A00",cex=1.8,xlim=c(0,100),ylim=c(0,100),xaxt="n",yaxt="n",axes=FALSE,ylab="",xlab="")
	lines(c(2.9,10.8,18.7,26.3,34.6,42,50,58,65.9,73.9,81.5,90,97.4),(PCRdup$FFPE),type="l",lty=1,lwd=12,cex=1.8,col="purple")
	legend(1.4,106.5,c("  Total Mean Coverage FF","  Total Mean Coverage FFPE"),bty="n",col=c("#9ACD32","#1E90FF"),fill=c("#9ACD32","#1E90FF"),cex=1.7,border=c("#9ACD32","#1E90FF"),text.col=textcolor)
	legend(0,95.8,c("% PCR Duplicates FF","% PCR Duplicates FFPE"),bty="n",col=c("#EE9A00","purple"),lty=c(1,1),lwd=7,cex=1.7, pt.lwd=c(5,5),text.col=textcolor)
dev.off()

say=sprintf("%s created!",pdffile);
print(say);
