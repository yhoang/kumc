SAMPLE = array(0,13)
SAMPLE  =  c("2474","2561","2640","2685","2938","3050","3356","4079","4191","CN","CT","DN","DT")

SAMPLE2 = array(0,13)
SAMPLE2[10] = "14119(N)";
SAMPLE2[11] = "14119(T)";
SAMPLE2[12] = "22285(N)";
SAMPLE2[13] = "22285(T)";

textcolor = "#2E2E2E";

fileFF1 = sprintf("/project/results/me/TES/mapQ/%s_FF2_cut",SAMPLE[1]);
fileFFPE1 = sprintf("/project/results/me/TES/mapQ/%s_FFPE2_cut",SAMPLE[1]);
tableFF1=read.table(fileFF1,sep="\t",header=T,fill=T);
tableFFPE1=read.table(fileFFPE1,sep="\t",header=T,fill=T);

fileFF2 = sprintf("/project/results/me/TES/mapQ/%s_FF2_cut",SAMPLE[2]);
fileFFPE2 = sprintf("/project/results/me/TES/mapQ/%s_FFPE2_cut",SAMPLE[2]);
tableFF2=read.table(fileFF2,sep="\t",header=T,fill=T);
tableFFPE2=read.table(fileFFPE2,sep="\t",header=T,fill=T);

fileFF3 = sprintf("/project/results/me/TES/mapQ/%s_FF2_cut",SAMPLE[3]);
fileFFPE3 = sprintf("/project/results/me/TES/mapQ/%s_FFPE2_cut",SAMPLE[3]);
tableFF3=read.table(fileFF3,sep="\t",header=T,fill=T);
tableFFPE3=read.table(fileFFPE3,sep="\t",header=T,fill=T);

fileFF4 = sprintf("/project/results/me/TES/mapQ/%s_FF2_cut",SAMPLE[4]);
fileFFPE4 = sprintf("/project/results/me/TES/mapQ/%s_FFPE2_cut",SAMPLE[4]);
tableFF4=read.table(fileFF4,sep="\t",header=T,fill=T);
tableFFPE4=read.table(fileFFPE4,sep="\t",header=T,fill=T);

fileFF5 = sprintf("/project/results/me/TES/mapQ/%s_FF2_cut",SAMPLE[5]);
fileFFPE5 = sprintf("/project/results/me/TES/mapQ/%s_FFPE2_cut",SAMPLE[5]);
tableFF5=read.table(fileFF5,sep="\t",header=T,fill=T);
tableFFPE5=read.table(fileFFPE5,sep="\t",header=T,fill=T);

fileFF6 = sprintf("/project/results/me/TES/mapQ/%s_FF2_cut",SAMPLE[6]);
fileFFPE6 = sprintf("/project/results/me/TES/mapQ/%s_FFPE2_cut",SAMPLE[6]);
tableFF6=read.table(fileFF6,sep="\t",header=T,fill=T);
tableFFPE6=read.table(fileFFPE6,sep="\t",header=T,fill=T);

fileFF7 = sprintf("/project/results/me/TES/mapQ/%s_FF2_cut",SAMPLE[7]);
fileFFPE7 = sprintf("/project/results/me/TES/mapQ/%s_FFPE2_cut",SAMPLE[7]);
tableFF7=read.table(fileFF7,sep="\t",header=T,fill=T);
tableFFPE7=read.table(fileFFPE7,sep="\t",header=T,fill=T);

fileFF8 = sprintf("/project/results/me/TES/mapQ/%s_FF2_cut",SAMPLE[8]);
fileFFPE8 = sprintf("/project/results/me/TES/mapQ/%s_FFPE2_cut",SAMPLE[8]);
tableFF8=read.table(fileFF8,sep="\t",header=T,fill=T);
tableFFPE8=read.table(fileFFPE8,sep="\t",header=T,fill=T);

fileFF9 = sprintf("/project/results/me/TES/mapQ/%s_FF2_cut",SAMPLE[9]);
fileFFPE9 = sprintf("/project/results/me/TES/mapQ/%s_FFPE2_cut",SAMPLE[9]);
tableFF9=read.table(fileFF9,sep="\t",header=T,fill=T);
tableFFPE9=read.table(fileFFPE9,sep="\t",header=T,fill=T);

fileFF10 = sprintf("/project/results/me/TES/mapQ/%s_FF2_cut",SAMPLE[10]);
fileFFPE10 = sprintf("/project/results/me/TES/mapQ/%s_FFPE2_cut",SAMPLE[10]);
tableFF10=read.table(fileFF10,sep="\t",header=T,fill=T);
tableFFPE10=read.table(fileFFPE10,sep="\t",header=T,fill=T);

fileFF11 = sprintf("/project/results/me/TES/mapQ/%s_FF2_cut",SAMPLE[11]);
fileFFPE11 = sprintf("/project/results/me/TES/mapQ/%s_FFPE2_cut",SAMPLE[11]);
tableFF11=read.table(fileFF11,sep="\t",header=T,fill=T);
tableFFPE11=read.table(fileFFPE11,sep="\t",header=T,fill=T);

fileFF12 = sprintf("/project/results/me/TES/mapQ/%s_FF2_cut",SAMPLE[12]);
fileFFPE12 = sprintf("/project/results/me/TES/mapQ/%s_FFPE2_cut",SAMPLE[12]);
tableFF12=read.table(fileFF12,sep="\t",header=T,fill=T);
tableFFPE12=read.table(fileFFPE12,sep="\t",header=T,fill=T);

fileFF13 = sprintf("/project/results/me/TES/mapQ/%s_FF2_cut",SAMPLE[13]);
fileFFPE13 = sprintf("/project/results/me/TES/mapQ/%s_FFPE2_cut",SAMPLE[13]);
tableFF13=read.table(fileFF13,sep="\t",header=T,fill=T);
tableFFPE13=read.table(fileFFPE13,sep="\t",header=T,fill=T);

# border=c("#9ACD32","#1E90FF"),
pdffile = sprintf("/project/results/me/TES/plots/mapQ_boxplot/TES_FF_FFPE_boxplot.pdf");
pdf (file = pdffile, width=10, height = 4, pointsize=1)
	par (las = 1,mar=c(5,5,1,1),oma=c(0,2,0,0),col=textcolor)
	boxplot ( tableFF1$AverageMapQuality, tableFFPE1$AverageMapQuality, tableFF2$AverageMapQuality, tableFFPE2$AverageMapQuality,  tableFF3$AverageMapQuality, tableFFPE3$AverageMapQuality, 
	 tableFF4$AverageMapQuality, tableFFPE4$AverageMapQuality, tableFF5$AverageMapQuality, tableFFPE5$AverageMapQuality, tableFF6$AverageMapQuality, tableFFPE6$AverageMapQuality, 
	 tableFF7$AverageMapQuality, tableFFPE7$AverageMapQuality,  tableFF8$AverageMapQuality, tableFFPE8$AverageMapQuality,  tableFF9$AverageMapQuality, tableFFPE9$AverageMapQuality, 
	 tableFF10$AverageMapQuality, tableFFPE10$AverageMapQuality,  tableFF11$AverageMapQuality, tableFFPE11$AverageMapQuality, tableFF12$AverageMapQuality, tableFFPE12$AverageMapQuality, 
	 tableFF13$AverageMapQuality, tableFFPE13$AverageMapQuality, outline=F, xlab="", ylab="", ylim = c(33,62),
	 col=c("#9ACD32","#1E90FF"), fill=c("#9ACD32","#1E90FF"), cex.main=1, cex=1,cex.axis=2, col.lab=textcolor, col.axis=textcolor, adj=1)
	legend("bottomright", c("FF","FFPE"), col=c("#9ACD32","#1E90FF"),fill=c("#9ACD32","#1E90FF"), cex=2.1, text.col=textcolor, bty="n")
	mtext (sprintf("mapping quality"),side=2,line=4,cex=2.5,las=0,col=textcolor,font=2)
	mtext("         2474     2561     2640     2685     2938     3050     3356     4079      4191     (N)14119(T)       (N)22285(T)", side=1,line=2, cex=2.2, col=textcolor,adj=0, font =2)
dev.off()


say=sprintf("%s created!",pdffile);
print(say); 


