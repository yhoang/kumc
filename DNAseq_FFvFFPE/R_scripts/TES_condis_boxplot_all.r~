textcolor = "#2E2E2E";
SAMPLE = array(0,13)
SAMPLE  =  c("2474","2561","2640","2685","2938","3050","3356","4079","4191","CN","CT","DN","DT")
max_cov = 100;

threshold = "recal";
target = "targeted";
SAMPLE2 = array(0,13)
SAMPLE2[10] = "14119(N)";
SAMPLE2[11] = "14119(T)";
SAMPLE2[12] = "22285(N)";
SAMPLE2[13] = "22285(T)";

textcolor = "#2E2E2E";

## READ TABLE
file_con = sprintf("/project/results/me/TES/SNV_comparison/genotyper/%s/%s_FF_FFPE_%s_con",target,SAMPLE[1],threshold);
file_dis = sprintf("/project/results/me/TES/SNV_comparison/genotyper/%s/%s_FF_FFPE_%s_dis",target,SAMPLE[1],threshold);
table_con1=read.table(file_con,sep="\t",header=T,fill=T);
table_dis1=read.table(file_dis,sep="\t",header=T,fill=T);
len_con1=length(table_con1[,1])
len_dis1=length(table_dis1[,1])

file_con = sprintf("/project/results/me/TES/SNV_comparison/genotyper/%s/%s_FF_FFPE_%s_con",target,SAMPLE[2],threshold);
file_dis = sprintf("/project/results/me/TES/SNV_comparison/genotyper/%s/%s_FF_FFPE_%s_dis",target,SAMPLE[2],threshold);
table_con2=read.table(file_con,sep="\t",header=T,fill=T);
table_dis2=read.table(file_dis,sep="\t",header=T,fill=T);
len_con2=length(table_con2[,2])
len_dis2=length(table_dis2[,2])

file_con = sprintf("/project/results/me/TES/SNV_comparison/genotyper/%s/%s_FF_FFPE_%s_con",target,SAMPLE[3],threshold);
file_dis = sprintf("/project/results/me/TES/SNV_comparison/genotyper/%s/%s_FF_FFPE_%s_dis",target,SAMPLE[3],threshold);
table_con3=read.table(file_con,sep="\t",header=T,fill=T);
table_dis3=read.table(file_dis,sep="\t",header=T,fill=T);
len_con3=length(table_con3[,3])
len_dis3=length(table_dis3[,3])

file_con = sprintf("/project/results/me/TES/SNV_comparison/genotyper/%s/%s_FF_FFPE_%s_con",target,SAMPLE[4],threshold);
file_dis = sprintf("/project/results/me/TES/SNV_comparison/genotyper/%s/%s_FF_FFPE_%s_dis",target,SAMPLE[4],threshold);
table_con4=read.table(file_con,sep="\t",header=T,fill=T);
table_dis4=read.table(file_dis,sep="\t",header=T,fill=T);
len_con4=length(table_con4[,4])
len_dis4=length(table_dis4[,4])

file_con = sprintf("/project/results/me/TES/SNV_comparison/genotyper/%s/%s_FF_FFPE_%s_con",target,SAMPLE[5],threshold);
file_dis = sprintf("/project/results/me/TES/SNV_comparison/genotyper/%s/%s_FF_FFPE_%s_dis",target,SAMPLE[5],threshold);
table_con5=read.table(file_con,sep="\t",header=T,fill=T);
table_dis5=read.table(file_dis,sep="\t",header=T,fill=T);
len_con5=length(table_con5[,5])
len_dis5=length(table_dis5[,5])

file_con = sprintf("/project/results/me/TES/SNV_comparison/genotyper/%s/%s_FF_FFPE_%s_con",target,SAMPLE[6],threshold);
file_dis = sprintf("/project/results/me/TES/SNV_comparison/genotyper/%s/%s_FF_FFPE_%s_dis",target,SAMPLE[6],threshold);
table_con6=read.table(file_con,sep="\t",header=T,fill=T);
table_dis6=read.table(file_dis,sep="\t",header=T,fill=T);
len_con6=length(table_con6[,6])
len_dis6=length(table_dis6[,6])

file_con = sprintf("/project/results/me/TES/SNV_comparison/genotyper/%s/%s_FF_FFPE_%s_con",target,SAMPLE[7],threshold);
file_dis = sprintf("/project/results/me/TES/SNV_comparison/genotyper/%s/%s_FF_FFPE_%s_dis",target,SAMPLE[7],threshold);
table_con7=read.table(file_con,sep="\t",header=T,fill=T);
table_dis7=read.table(file_dis,sep="\t",header=T,fill=T);
len_con7=length(table_con7[,7])
len_dis7=length(table_dis7[,7])

file_con = sprintf("/project/results/me/TES/SNV_comparison/genotyper/%s/%s_FF_FFPE_%s_con",target,SAMPLE[8],threshold);
file_dis = sprintf("/project/results/me/TES/SNV_comparison/genotyper/%s/%s_FF_FFPE_%s_dis",target,SAMPLE[8],threshold);
table_con8=read.table(file_con,sep="\t",header=T,fill=T);
table_dis8=read.table(file_dis,sep="\t",header=T,fill=T);
len_con8=length(table_con8[,8])
len_dis8=length(table_dis8[,8])

file_con = sprintf("/project/results/me/TES/SNV_comparison/genotyper/%s/%s_FF_FFPE_%s_con",target,SAMPLE[9],threshold);
file_dis = sprintf("/project/results/me/TES/SNV_comparison/genotyper/%s/%s_FF_FFPE_%s_dis",target,SAMPLE[9],threshold);
table_con9=read.table(file_con,sep="\t",header=T,fill=T);
table_dis9=read.table(file_dis,sep="\t",header=T,fill=T);
len_con9=length(table_con9[,9])
len_dis9=length(table_dis9[,9])

file_con = sprintf("/project/results/me/TES/SNV_comparison/genotyper/%s/%s_FF_FFPE_%s_con",target,SAMPLE[10],threshold);
file_dis = sprintf("/project/results/me/TES/SNV_comparison/genotyper/%s/%s_FF_FFPE_%s_dis",target,SAMPLE[10],threshold);
table_con10=read.table(file_con,sep="\t",header=T,fill=T);
table_dis10=read.table(file_dis,sep="\t",header=T,fill=T);
len_con10=length(table_con10[,10])
len_dis10=length(table_dis10[,10])

file_con = sprintf("/project/results/me/TES/SNV_comparison/genotyper/%s/%s_FF_FFPE_%s_con",target,SAMPLE[11],threshold);
file_dis = sprintf("/project/results/me/TES/SNV_comparison/genotyper/%s/%s_FF_FFPE_%s_dis",target,SAMPLE[11],threshold);
table_con11=read.table(file_con,sep="\t",header=T,fill=T);
table_dis11=read.table(file_dis,sep="\t",header=T,fill=T);
len_con11=length(table_con11[,11])
len_dis11=length(table_dis11[,11])

file_con = sprintf("/project/results/me/TES/SNV_comparison/genotyper/%s/%s_FF_FFPE_%s_con",target,SAMPLE[12],threshold);
file_dis = sprintf("/project/results/me/TES/SNV_comparison/genotyper/%s/%s_FF_FFPE_%s_dis",target,SAMPLE[12],threshold);
table_con12=read.table(file_con,sep="\t",header=T,fill=T);
table_dis12=read.table(file_dis,sep="\t",header=T,fill=T);
len_con12=length(table_con12[,12])
len_dis12=length(table_dis12[,12])

file_con = sprintf("/project/results/me/TES/SNV_comparison/genotyper/%s/%s_FF_FFPE_%s_con",target,SAMPLE[13],threshold);
file_dis = sprintf("/project/results/me/TES/SNV_comparison/genotyper/%s/%s_FF_FFPE_%s_dis",target,SAMPLE[13],threshold);
table_con13=read.table(file_con,sep="\t",header=T,fill=T);
table_dis13=read.table(file_dis,sep="\t",header=T,fill=T);
len_con13=length(table_con13[,13])
len_dis13=length(table_dis13[,13])



####### PLOT
pdffile = sprintf("/project/results/me/TES/plots/SNV_boxplot/TES_FF_FFPE_%s_boxplot_all.pdf",threshold);
pdf (file = pdffile, width=15, height = 4, pointsize=1)
	par (las = 1, mar=c(7,6,1,2),oma=c(0,2,0,0),col=textcolor)
	
	boxplot ( table_con1$Depth.1, table_dis1$Depth.1, table_con2$Depth.1, table_dis2$Depth.1, table_con3$Depth.1, table_dis3$Depth.1, 
	 table_con4$Depth.1, table_dis4$Depth.1, table_con5$Depth.1, table_dis5$Depth.1, table_con6$Depth.1, table_dis6$Depth.1, 
	 table_con7$Depth.1, table_dis7$Depth.1, table_con8$Depth.1, table_dis8$Depth.1, table_con9$Depth.1, table_dis9$Depth.1, 
	 table_con10$Depth.1, table_dis10$Depth.1, table_con11$Depth.1, table_dis11$Depth.1, table_con12$Depth.1, table_dis12$Depth.1, 
	 table_con13$Depth.1, table_dis13$Depth.1,
	 ylim=c(0,max_cov) , xaxt = "n", space = 2, boxwex=0.7,
	 col=c("#999999","#FF4500"), cex.axis =2 )
	mtext (sprintf("coverage"),side=2,line=4,cex=2.5,las=0,col=textcolor, font = 2)
	legend("top", c("concordant","discordant"), col=c("#999999","#FF4500"),fill=c("#999999","#FF4500"), cex=2.1, text.col=textcolor, bg="white")
	abline(h=13,col=textcolor,lwd=3)
	mtext("              2474         2561         2640         2685          2938         3050         3356         4079         4191          (N)  14119  (T)            (N)  22285  (T)", side=1,line=1, cex=2.5, col=textcolor,adj=0, font =2)
	mtext(sprintf("            (%s,%s)   (%s,%s)   (%s,%s)   (%s,%s)  (%s,%s)   (%s,%s)  (%s,%s)  (%s,%s)  (%s,%s)   (%s,%s) (%s,%s)   (%s,%s) (%s,%s)",
	 len_con1, len_dis1, len_con2, len_dis2, len_con3, len_dis3, len_con4, len_dis4, len_con5, len_dis5, len_con6, len_dis6, len_con7, len_dis7,
	 len_con8, len_dis8, len_con9, len_dis9, len_con10, len_dis10, len_con11, len_dis11, len_con12, len_dis12, len_con13, len_dis13 ), 
	 side=1, line = 3, cex=2.08, col=textcolor, adj=0 )
	box("plot", col=textcolor)
dev.off()


say=sprintf("%s created!",pdffile);
print(say); 

