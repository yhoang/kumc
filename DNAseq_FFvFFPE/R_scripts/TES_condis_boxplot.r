targeted  =  "";
SAMPLE = array(0,13)
SAMPLE  =  c("2474","2561","2640","2685","2938","3050","3356","4079","4191","CN","CT","DN","DT")
targeted = "targeted_";

threshold = "recal";
SAMPLE2 = array(0,13)
SAMPLE2[10] = "14119(N)";
SAMPLE2[11] = "14119(T)";
SAMPLE2[12] = "22285(N)";
SAMPLE2[13] = "22285(T)";

textcolor = "#2E2E2E";

file_con = sprintf("/project/results/me/TES/SNV_comparison/genotyper/targeted/%s_FF_FFPE_%s_%scon",SAMPLE[s],threshold,targeted);
file_dis = sprintf("/project/results/me/TES/SNV_comparison/genotyper/targeted/%s_FF_FFPE_%s_%sdis",SAMPLE[s],threshold,targeted);

table_con=read.table(file_con,sep="\t",header=T,fill=T);
table_dis=read.table(file_dis,sep="\t",header=T,fill=T);

pdffile = sprintf("/project/results/me/TES/plots/SNV_boxplot/TES_%s_FF_FFPE_%s_%sboxplot.pdf",SAMPLE[s],threshold,targeted);
pdf (file = pdffile, width=5, height = 5, pointsize=12)
par (las = 1,mar=c(4,1,1,1),oma=c(0,2,0,0),col=textcolor)
boxplot (table_con$RMS.mapping.quality.1,table_dis$RMS.mapping.quality.1,xlab="concordance (FFPE)      discordance (FFPE)",ylim=c(10,67), ylab="snp quality",col=c(3,2),cex.main=1,cex.lab=1,cex=1,cex.axis=.7,col.lab=textcolor,col.axis=textcolor)
mtext (sprintf("mapping quality"),side=2,line=2,cex=1,las=0,col=textcolor)
boxplot (table_con$QualbyDepth.1,table_dis$QualbyDepth.1,xlab="concordance (FFPE)      discordance (FFPE)",ylim=c(0,40), ylab="snp quality",col=c(3,2),cex.main=1,cex.lab=1,cex=1,cex.axis=.7,col.lab=textcolor,col.axis=textcolor)
mtext (sprintf("snp quality"),side=2,line=2,cex=1,las=0,col=textcolor)
boxplot (table_con$Depth.1,table_dis$Depth.1,ylim=c(0,600),xlab="concordance (FFPE)      discordance (FFPE)", ylab="coverage",col=c(3,2),cex.main=1,cex.lab=1,cex=1,cex.axis=.7,col.lab=textcolor,col.axis=textcolor)
mtext (sprintf("coverage"),side=2,line=2,cex=1,las=0,col=textcolor)
boxplot (table_con$Depth.1,table_dis$Depth.1,ylim=c(0,120),xlab="concordance (FFPE)      discordance (FFPE)", ylab="coverage",col=c(3,2),cex.main=1,cex.lab=1,cex=1,cex.axis=.7,col.lab=textcolor,col.axis=textcolor)
mtext (sprintf("coverage"),side=2,line=2,cex=1,las=0,col=textcolor)
dev.off()


say=sprintf("%s created!",pdffile);
print(say);

#############
targeted = "targeted";
threshold = "recal";

median_dis=array(0,13)
for (s in 1:13) {
	file_dis = sprintf("/project/results/me/TES/SNV_comparison/genotyper/base_call/%s_FF_FFPE_%s_%s_dis",SAMPLE[s],threshold,targeted);
	table_dis=read.table(file_dis,sep="\t",header=T,fill=T);
	median_dis[s]=median(table_dis$Depth.1)
}
median_dis[3]=0
median_dis[9]=0
median_dis[13]=0
y_tick = c((0:6)*30)
max=max(median_dis)+10
barplotfile = sprintf("/project/results/me/TES/plots/SNV_boxplot/TES_FF_FFPE_%s_%s_basecall_barplot.pdf",threshold,targeted);
pdf (file = barplotfile, width=30, height=20, pointsize=30)  
par (oma=c(0,0,0,0),mar=c(4,6,2,2),col=textcolor)
barplot(median_dis, beside=T,width=(0.008),ylim=c(0,200),col=c("#FF4500"),cex.lab=2.2, main = "",border="#FF4500",ylab="median coverage of discordant calls in FFPE",col.lab=textcolor,yaxt="n")
abline(h=13,col="#A3A3A3",lwd=5)
mtext("       2474    2561   2640    2685    2938   3050    3356   4079   4191   (N) 14119 (T) (N) 22285 (T)", side=1,line=1, cex=1.4, col=textcolor,adj=0, font =2)
axis(2,lty=1,lwd=7,at=(y_tick),labels=y_tick,las=0, cex.axis=2,col=textcolor,col.lab=textcolor,col.axis=textcolor)
mtext("          (1)       (4)         (0)        (2)         (3)        (1)        (2)        (4)       (0)        (3)         (2)       (4)         (0)", side=1, line = 2, cex=1.3, col=textcolor, adj=0)
dev.off()


say=sprintf("%s created!",barplotfile);
print(say);



#########################################
targeted = "targeted";
threshold = "recal";

median_dis=array(0,13)
for (s in 1:13) {
	file_dis = sprintf("/project/results/me/TES/SNV_comparison/genotyper/targeted/%s_FF_FFPE_%s_%s_dis",SAMPLE[s],threshold,targeted);
	table_dis=read.table(file_dis,sep="\t",header=T,fill=T);
	median_dis[s]=median(table_dis$Depth.1)
}
median_dis[3]=0

y_tick = c((0:6)*30)
max=max(median_dis)+10
barplotfile = sprintf("/project/results/me/TES/plots/SNV_boxplot/TES_FF_FFPE_%s_%s_barplot.pdf",threshold,targeted);
pdf (file = barplotfile, width=30, height=20, pointsize=30)  
par (oma=c(0,0,0,0),mar=c(4,6,2,2),col=textcolor)
barplot(median_dis, beside=T,width=(0.008),ylim=c(0,200),col=c("#FF4500"),cex.lab=2.2, main = "",border="#FF4500",ylab="median coverage of discordant calls in FFPE",col.lab=textcolor,yaxt="n")
abline(h=13,col="#A3A3A3",lwd=5)
mtext("       2474    2561   2640    2685    2938   3050    3356   4079   4191   (N) 14119 (T) (N) 22285 (T)", side=1,line=1, cex=1.4, col=textcolor,adj=0, font =2)
axis(2,lty=1,lwd=7,at=(y_tick),labels=y_tick,las=0, cex.axis=2,col=textcolor,col.lab=textcolor,col.axis=textcolor)
mtext("          (8)       (16)       (0)        (9)       (75)      (10)      (44)       (9)        (1)       (10)       (8)      (23)        (65)", side=1, line = 2, cex=1.3, col=textcolor, adj=0)
dev.off()


say=sprintf("%s created!",barplotfile);
print(say);



#############
targeted = "targeted";
threshold = "recal";

median_dis=array(0,13)
median_con=array(0,13)
for (s in 1:13) {
	file_con = sprintf("/project/results/me/TES/SNV_comparison/genotyper/targeted/%s_FF_FFPE_%s_%s_con",SAMPLE[s],threshold,targeted);
	file_dis = sprintf("/project/results/me/TES/SNV_comparison/genotyper/targeted/%s_FF_FFPE_%s_%s_dis",SAMPLE[s],threshold,targeted);
	table_con=read.table(file_con,sep="\t",header=T,fill=T);
	table_dis=read.table(file_dis,sep="\t",header=T,fill=T);
	median_dis[s]=median(table_dis$Depth.1)
	median_con[s]=median(table_con$Depth.1)
}

max=max(median_con)+10
barplotfile = sprintf("/project/results/me/TES/plots/SNV_boxplot/TES_FF_FFPE_%s_%s_barplot.pdf",threshold,targeted);
pdf (file = barplotfile, width=30, height=20, pointsize=30)  
par (oma=c(0,0,0,0),mar=c(3,5,2,5),col=textcolor)
barplot(median_con, beside=T,width=(0.008),ylim=c(0,max),col=c("#999999"),border="#999999",xaxt="n",yaxt="n", main = "",ylab="")
par(new=T)
barplot(median_dis, beside=T,width=(0.008),ylim=c(0,max),col=c("#FF4500"),cex.axis=1.5,cex.lab=1.5, main = "",border="#FF4500",ylab="median coverage of variant calls in FFPE",cex=1.4,col.lab=textcolor,col.axis=textcolor)
abline(h=13,col="#9ACD32",lwd=5)
mtext("       2474     2561    2640     2685     2938     3050     3356    4079     4191    (N) 14119 (T)    (N) 22285 (T)", side=1,line=1, cex=1.24, col=textcolor,adj=0)
legend("topright",c("concordant","discordant"), col=c("#999999","#FF4500"),bty="n",fill=c("#999999","#FF4500"), cex=1.5,border = c("#999999","#FF4500"),text.col=textcolor)
dev.off()

say=sprintf("%s created!",barplotfile);
print(say);



