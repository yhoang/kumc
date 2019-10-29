targeted  =  "";
SAMPLE = array(0,13)
SAMPLE  =  c("2474","2561","2640","2685","2938","3050","3356","4079","4191","CN","CT","DN","DT")
#targeted = "exonic";
max_snp = 1;
max_snp = 2;
max_cov = 20;
max_cov = 50;

threshold = "recal";
SAMPLE2 = array(0,13)
SAMPLE2[10] = "14119(N)";
SAMPLE2[11] = "14119(T)";
SAMPLE2[12] = "22285(N)";
SAMPLE2[13] = "22285(T)";

textcolor = "#2E2E2E";

file_con = sprintf("/project/results/me/WES/SNV_comparison/genotyper/total/%s_FFex_FFPEex_%s_%scon",SAMPLE[s],threshold,targeted);
file_dis = sprintf("/project/results/me/WES/SNV_comparison/genotyper/total/%s_FFex_FFPEex_%s_%sdis",SAMPLE[s],threshold,targeted);
#file_con = sprintf("/project/results/me/WES/SNV_comparison/genotyper/%s/%s_FFex_FFPEex_%s_con",targeted,SAMPLE[s],threshold);
#file_dis = sprintf("/project/results/me/WES/SNV_comparison/genotyper/%s/%s_FFex_FFPEex_%s_dis",targeted,SAMPLE[s],threshold);


table_con=read.table(file_con,sep="\t",header=T,fill=T);
table_dis=read.table(file_dis,sep="\t",header=T,fill=T);

pdffile = sprintf("/project/results/me/WES/plots/SNV_boxplot/WES_%s_FFex_FFPEex_%s_%s_boxplot.pdf",SAMPLE[s],threshold,targeted);
pdf (file = pdffile, width=5, height = 5, pointsize=12)
par (las = 1,mar=c(4,1,1,1),oma=c(0,2,0,0),col=textcolor)
boxplot (table_con$RMS.mapping.quality.1,table_dis$RMS.mapping.quality.1,xlab="concordance (FFPE)      discordance (FFPE)",ylim=c(0,70), ylab="snp quality",col=c("#999999","#FF4500"),cex.main=1,cex.lab=1,cex=1,cex.axis=.7,col.lab=textcolor,col.axis=textcolor)
mtext (sprintf("mapping quality"),side=2,line=2,cex=1,las=0,col=textcolor)
boxplot (table_con$QualbyDepth.1,table_dis$QualbyDepth.1,xlab="concordance (FFPE)      discordance (FFPE)",ylim=c(0,max_snp), ylab="snp quality",col=c("#999999","#FF4500"),cex.main=1,cex.lab=1,cex=1,cex.axis=.7,col.lab=textcolor,col.axis=textcolor)
mtext (sprintf("snp quality"),side=2,line=2,cex=1,las=0,col=textcolor)
boxplot (table_con$Depth.1,table_dis$Depth.1,ylim=c(0,300),xlab="concordance (FFPE)      discordance (FFPE)", ylab="coverage",col=c("#999999","#FF4500"),cex.main=1,cex.lab=1,cex=1,cex.axis=.7,col.lab=textcolor,col.axis=textcolor)
mtext (sprintf("coverage"),side=2,line=2,cex=1,las=0,col=textcolor)
boxplot (table_con$Depth.1,table_dis$Depth.1,ylim=c(0,max_cov),xlab="concordance (FFPE)      discordance (FFPE)", ylab="coverage",col=c("#999999","#FF4500"),cex.main=1,cex.lab=1,cex=1,cex.axis=.7,col.lab=textcolor,col.axis=textcolor)
mtext (sprintf("coverage"),side=2,line=2,cex=1,las=0,col=textcolor)
box("plot", col=textcolor)  
dev.off()


say=sprintf("%s created!",pdffile);
print(say); 

##################
targeted = "exonic";
threshold = "recal";

median_dis=array(0,13)
for (s in 1:13) {
	file_dis = sprintf("/project/results/me/WES/SNV_comparison/genotyper/base_call/%s_FFex_FFPEex_%s_%s_dis",SAMPLE[s],threshold,targeted);
	table_dis=read.table(file_dis,sep="\t",header=T,fill=T);
	median_dis[s]=median(table_dis$QualbyDepth.1)
}

y_tick = c((0:8)*30)
max=max(median_dis)+10
barplotfile = sprintf("/project/results/me/WES/plots/SNV_boxplot/WES_FF_FFPE_%s_%s_basecall_barplot.pdf",threshold,targeted);
pdf (file = barplotfile, width=30, height=20, pointsize=30)  
par (oma=c(0,0,0,0),mar=c(4,6,2,2),col=textcolor)
barplot(median_dis, beside=T,width=(0.008),ylim=c(0,250),col=c("#FF4500"),cex.lab=2.2, main = "",border="#FF4500",ylab="median coverage of discordant calls in FFPE",col.lab=textcolor,yaxt="n")
abline(h=13,col="#A3A3A3",lwd=5)
mtext("       2474    2561   2640    2685    2938   3050    3356   4079   4191   (N) 14119 (T) (N) 22285 (T)", side=1,line=1, cex=1.4, col=textcolor,adj=0, font =2)
axis(2,lty=1,lwd=7,at=(y_tick),labels=y_tick,las=0, cex.axis=2,col=textcolor,col.lab=textcolor,col.axis=textcolor)
mtext("          (3)       (4)         (2)        (3)         (7)        (6)        (6)       (10)       (8)        (5)         (3)       (5)         (3)", side=1, line = 2, cex=1.3, col=textcolor, adj=0)
dev.off()


say=sprintf("%s created!",barplotfile);
print(say);



#########################
targeted = "exonic";
threshold = "recal";

median_dis=array(0,13)
for (s in 1:13) {
	file_dis = sprintf("/project/results/me/WES/SNV_comparison/genotyper/%s/%s_FFex_FFPEex_%s_dis.tsv",targeted,SAMPLE[s],threshold);
	table_dis=read.table(file_dis,sep="\t",header=T,fill=T);
	median_dis[s]=median(table_dis$Depth.1)
}

max=max(median_dis)+5
y_tick = c((0:6)*30)
barplotfile = sprintf("/project/results/me/WES/plots/SNV_boxplot/WES_FFex_FFPEex_%s_%s_barplot.pdf",threshold,targeted);
pdf (file = barplotfile, width=30, height=20, pointsize=30)   
par (oma=c(0,0,0,0),mar=c(4,6,2,2),col=textcolor)
barplot(median_dis, beside=T,width=(0.008),ylim=c(0,200),col=c("#FF4500"),cex.lab=2.2, main = "",border="#FF4500",ylab="median coverage of discordant calls in FFPE",col.lab=textcolor,yaxt="n")
abline(h=13,col="#A3A3A3",lwd=5)
mtext("       2474    2561   2640    2685    2938   3050    3356   4079   4191   (N) 14119 (T) (N) 22285 (T)", side=1,line=1, cex=1.4, col=textcolor,adj=0, font =2)
axis(2,lty=1,lwd=7,at=(y_tick),labels=y_tick,las=0, cex.axis=2,col=textcolor,col.lab=textcolor,col.axis=textcolor)
mtext("         (10)      (86)       (7)        (8)      (616)     (15)     (434)     (10)      (37)      (20)      (43)      (19)      (209)", side=1, line = 2, cex=1.3, col=textcolor, adj=0)
dev.off()

say=sprintf("%s created!",barplotfile);
print(say);








#####################
targeted = "exonic";
threshold = "recal";

median_dis=array(0,13)
median_con=array(0,13)
for (s in 1:13) {
	file_con = sprintf("/project/results/me/WES/SNV_comparison/genotyper/%s/%s_FFex_FFPEex_%s_con",targeted,SAMPLE[s],threshold);
	file_dis = sprintf("/project/results/me/WES/SNV_comparison/genotyper/%s/%s_FFex_FFPEex_%s_dis",targeted,SAMPLE[s],threshold);
	table_con=read.table(file_con,sep="\t",header=T,fill=T);
	table_dis=read.table(file_dis,sep="\t",header=T,fill=T);
	median_dis[s]=median(table_dis$Depth.1)
	median_con[s]=median(table_con$Depth.1)
}

max=100
barplotfile = sprintf("/project/results/me/WES/plots/SNV_boxplot/WES_FFex_FFPEex_%s_%s_barplot.pdf",threshold,targeted);
pdf (file = barplotfile, width=30, height=20, pointsize=30)  
par (oma=c(0,0,0,0),mar=c(3,5,2,5),col=textcolor)
barplot(median_con, beside=T,width=(0.008),ylim=c(0,max),col=c("#999999"),border="#999999",xaxt="n",yaxt="n", main = "",ylab="")
par(new=T)
barplot(median_dis, beside=T,width=(0.008),ylim=c(0,max),col=c("#FF4500"),cex.axis=1.5,cex.lab=1.5, main = "",border="#FF4500",ylab="median coverage of variant calls in FFPE",cex=1.4,col.lab=textcolor,col.axis=textcolor)
abline(h=13,col="#9ACD32",lwd=5)
mtext("       2474     2561    2640     2685     2938     3050     3356    4079     4191    (N) 14119 (T)    (N) 22285 (T)", side=1,line=1, cex=1.24, col=textcolor,adj=0)
legend("topright",c("concordant","discordant"), col=c("#999999","#FF4500"),bty="n",fill=c("#999999","#FF4500"), cex=1.5,border = c("#999999","#FF4500"),text.col=textcolor)
legend( 0, 0, "coverage cutoff of 13", col="#9ACD32", lty=1, lwd=5, cex=1.5, text.col=textcolor)
dev.off()


