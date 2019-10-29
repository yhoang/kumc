#!/bin/bash

threshold = "q43cov13";
threshold = "pass";


SNV_FILE = sprintf("/project/results/me/WES/plots/intersection_barplot/SNV_targeted_%s",threshold);
pdffile = sprintf("/project/results/me/WES/plots/intersection_barplot/targeted_intersection_SNV_%s.pdf",threshold);
table=read.table(SNV_FILE,sep="\t",header=TRUE);

max=1600
pdf (file = pdffile, width=30, height=20, pointsize=30)  
par (oma=c(0,0,0,0),mar=c(3,5,2,5),col=textcolor)
barplot(rbind(table[,2],table[,3]),beside=T, width=(0.02),ylim=c(0,max),col=c("#9ACD32","#1E90FF"),cex.axis=1.5,cex.lab=1.5, main = "",border="dark grey",ylab="SNV counts",cex=1.4,col.lab=textcolor,col.axis=textcolor)
#names.arg=c("2474","2561","2640","2685","2938","3050","3356","4079","4191","14119(N)","14119(T)","22285(N)","22285(T)"),ylab="SNV counts",cex=1.4,col.lab=textcolor,col.axis=textcolor)
mtext("       2474     2561    2640     2685     2938     3050     3356    4079    4191      (N) 14119 (T)    (N) 22285 (T)", side=1,line=1, cex=1.24, col=textcolor,adj=0)
legend("topright",c("FF","FFPE"),col=c("#9ACD32","#1E90FF"),fill=c("#9ACD32","#1E90FF"),cex=1.5,border="dark grey",text.col=textcolor,bg = "white")
dev.off()

say=sprintf("%s created!",pdffile);
print(say);   





textcolor = "#2E2E2E";
SNV_FILE = sprintf("/project/results/me/WES/plots/intersection_barplot/SNV_%s",threshold);
pdffile = sprintf("/project/results/me/WES/plots/intersection_barplot/intersection_SNV_%s.pdf",threshold);

table=read.table(SNV_FILE,sep="\t",header=TRUE);
max=8000
pdf (file = pdffile, width=30, height=20, pointsize=30)  
par (oma=c(0,0,0,0),mar=c(3,5,2,5),col=textcolor)
barplot(rbind(table[,2],table[,3]), beside=T,width=(0.008),ylim=c(0,max),col=c("#9ACD32","#1E90FF"),cex.axis=1.5,cex.lab=1.5, main = "",border="dark grey",ylab="SNV counts",cex=1.4,col.lab=textcolor,col.axis=textcolor)
#names.arg=c("2474","2561","2640","2685","2938","3050","3356","4079","4191","14119(N)","14119(T)","22285(N)","22285(T)"),ylab="SNV counts",cex=1.4,col.lab=textcolor,col.axis=textcolor)
mtext("       2474     2561    2640     2685     2938     3050     3356    4079    4191      (N) 14119 (T)    (N) 22285 (T)", side=1,line=1, cex=1.24, col=textcolor,adj=0)
legend("topright",c("FF","FFPE"), col=c("#9ACD32","#1E90FF"),fill=c("#9ACD32","#1E90FF"),cex=1.5,border="dark grey",text.col=textcolor,bg = "white")
dev.off()

say=sprintf("%s created!",pdffile);
print(say);   

