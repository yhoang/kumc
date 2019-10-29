#!/bin/bash

threshold = "recal";
#threshold = "q43cov13";
#threshold = "pass";


max=2200
targeted = "targeted_";

textcolor = "#2E2E2E";
SNV_FILE = sprintf("/project/results/me/TES/plots/SNV_condis_barplot/TES_SNV_%s",threshold);
pdffile = sprintf("/project/results/me/TES/plots/SNV_condis_barplot/TES_targeted_SNV_condis_%s.pdf",threshold);

table=read.table(SNV_FILE,sep="\t",header=TRUE);

y_tick = c((0:5)*400)
pdf (file = pdffile, width=30, height=20, pointsize=30)  
par (oma=c(0,0,0,0),mar=c(3,6,2,2),col=textcolor)
barplot(rbind(table[,10],table[,11],table[,13],table[,12]), width=(0.02),ylim=c(0,max),col=c("#999999","#FF4500","#1E90FF","#9ACD32"),cex.lab=2.4, font.lab = 2, main = "",border="dark grey",ylab="SNV counts",col.lab=textcolor,yaxt="n")
axis(2,lty=1,lwd=7,at=(y_tick),labels=y_tick,las=0, cex.axis=2,col=textcolor,col.lab=textcolor,col.axis=textcolor)
mtext("       2474    2561   2640    2685   2938   3050    3356    4079    4191   (N)14119(T)   (N)22285(T)", side=1,line=1, cex=1.4, col=textcolor,adj=0, font =2)
legend("topright",c("concordant","discordant","false positives","false negatives"),col=c("#999999","#FF4500","#1E90FF","#9ACD32"),
fill=c("#999999","#FF4500","#1E90FF","#9ACD32"),cex=1.7,border="dark grey",text.col=textcolor,bg = "white")
dev.off()

say=sprintf("%s created!",pdffile);
print(say);   

########

max=280000
y_tick = c((0:6)*50000)
SNV_FILE = sprintf("/project/results/me/TES/plots/SNV_condis_barplot/TES_SNV_%s",threshold);
pdffile = sprintf("/project/results/me/TES/plots/SNV_condis_barplot/TES_SNV_condis_%s.pdf",threshold);

table=read.table(SNV_FILE,sep="\t",header=TRUE);

pdf (file = pdffile, width=30, height=20, pointsize=30)  
par (oma=c(0,0,0,0),mar=c(3,6,2,2),col=textcolor)
barplot(rbind(table[,4],table[,5],table[,7],table[,6]), width=(0.02),ylim=c(0,max),col=c("#999999","#FF4500","#1E90FF","#9ACD32"),cex.lab=2.4, font.lab = 2, main = "",border="dark grey",ylab="SNV counts",col.lab=textcolor,yaxt="n")
axis(2,lty=1,lwd=7,at=(y_tick),labels=y_tick,las=0, cex.axis=2,col=textcolor,col.lab=textcolor,col.axis=textcolor)
mtext("       2474    2561   2640    2685   2938   3050    3356    4079   4191    (N)14119(T)   (N)22285(T)", side=1,line=1, cex=1.4, col=textcolor,adj=0, font =2)
legend("topright",c("concordant","discordant","false positives","false negatives"),col=c("#999999","#FF4500","#1E90FF","#9ACD32"),
fill=c("#999999","#FF4500","#1E90FF","#9ACD32"),cex=1.7,border="dark grey",text.col=textcolor,bg = "white")
dev.off()

say=sprintf("%s created!",pdffile);
print(say);   

