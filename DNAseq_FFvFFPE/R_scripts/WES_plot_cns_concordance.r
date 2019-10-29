#!/bin/bash
print("Set s (1-13)");
SAMPLE=array(0,13)
targeted = "targeted_";
#targeted = "";
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


file = sprintf("/project/results/me/WES/cns_concordance/%s_FFex_FFPEex_%scns_cov",SAMPLE[s],targeted);
pdffile =  sprintf("/project/results/me/WES/plots/cns_concordance/%s%s_FFex_FFPEex_cns_cov.pdf",targeted,SAMPLE[s]);

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

