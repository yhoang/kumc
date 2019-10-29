#!/bin/bash
print("Set s (1-13)");
targeted = "targeted_";
#targeted = "";


SAMPLE = array(0,13)
SAMPLE  =  c("2474","2561","2640","2685","2938","3050","3356","4079","4191","CN","CT","DN","DT")

SAMPLE2=array(0,13)
SAMPLE2[10]="14119(N)";
SAMPLE2[11]="14119(T)";
SAMPLE2[12]="22285(N)";
SAMPLE2[13]="22285(T)";

color_FF = array(0,13)
color_FF  =  c()

color_FFPE = array(0,13)
color_FFPE  =  c()

file = sprintf("/project/results/me/TES/cns_concordance/%s_FF_FFPE_%scns_cov",SAMPLE[s],targeted);
pdffile =  sprintf("/project/results/me/TES/plots/cns_concordance/%s%s_FF_FFPE_cns_cov.pdf",targeted,SAMPLE[s]);

con=read.table(file,sep="\t",header=T);

pdf (file = pdffile, width=300, height=200, pointsize=300)
par (oma=c(0,0,0,0),mar=c(5,5,1,1),col=textcolor)
plot(con$cov,(100*con$conFF),lwd=30,xlim=c(0,100),ylim=c(90,100),col="#9ACD32",pch=3,cex=1.5,xlab="minimal coverage",ylab = "% concordant positions",cex.lab=2,cex.axis=2,col.axis=textcolor,col.lab=textcolor)
par(new=T)
plot(con$cov,(100*con$conFFPE),xlim=c(0,100),ylim=c(90,100),pch=3,cex=1.5,xaxt="n",yaxt="n",ylab="",xlab="",col=textcolor)
abline(h=95,col="#A3A3A3",lwd=7)
abline(h=99,col="#A3A3A3",lwd=7)
box("plot", col=textcolor)
legend("bottom",c("FF","FFPE"),col=c("#9ACD32",textcolor), pch=c(3,3), title="preparation",cex=1.5,pt.cex=1.5,border = textcolor,text.col=textcolor)
dev.off()

say=sprintf("%s created!",pdffile);	
print(say);

