SAMPLE = array(0,13)
SAMPLE  =  c("2474","2561","2640","2685","2938","3050","3356","4079","4191","CN","CT","DN","DT")

SAMPLE2 = array(0,13)
SAMPLE2[10] = "14119(N)";
SAMPLE2[11] = "14119(T)";
SAMPLE2[12] = "22285(N)";
SAMPLE2[13] = "22285(T)";

textcolor = "#2E2E2E";

fileFF1 = sprintf("/project/results/me/WES/mapQ/%s_FFex2_cut10",SAMPLE[1]);
fileFFPE1 = sprintf("/project/results/me/WES/mapQ/%s_FFPEex2_cut10",SAMPLE[1]);
tableFF1=read.table(fileFF1,sep="\t",header=F,fill=T);
tableFFPE1=read.table(fileFFPE1,sep="\t",header=F,fill=T);

fileFF2 = sprintf("/project/results/me/WES/mapQ/%s_FFex2_cut10",SAMPLE[2]);
fileFFPE2 = sprintf("/project/results/me/WES/mapQ/%s_FFPEex2_cut10",SAMPLE[2]);
tableFF2=read.table(fileFF2,sep="\t",header=F,fill=T);
tableFFPE2=read.table(fileFFPE2,sep="\t",header=F,fill=T);

fileFF3 = sprintf("/project/results/me/WES/mapQ/%s_FFex2_cut10",SAMPLE[3]);
fileFFPE3 = sprintf("/project/results/me/WES/mapQ/%s_FFPEex2_cut10",SAMPLE[3]);
tableFF3=read.table(fileFF3,sep="\t",header=F,fill=T);
tableFFPE3=read.table(fileFFPE3,sep="\t",header=F,fill=T);

fileFF4 = sprintf("/project/results/me/WES/mapQ/%s_FFex2_cut10",SAMPLE[4]);
fileFFPE4 = sprintf("/project/results/me/WES/mapQ/%s_FFPEex2_cut10",SAMPLE[4]);
tableFF4=read.table(fileFF4,sep="\t",header=F,fill=T);
tableFFPE4=read.table(fileFFPE4,sep="\t",header=F,fill=T);

fileFF5 = sprintf("/project/results/me/WES/mapQ/%s_FFex2_cut10",SAMPLE[5]);
fileFFPE5 = sprintf("/project/results/me/WES/mapQ/%s_FFPEex2_cut10",SAMPLE[5]);
tableFF5=read.table(fileFF5,sep="\t",header=F,fill=T);
tableFFPE5=read.table(fileFFPE5,sep="\t",header=F,fill=T);

fileFF6 = sprintf("/project/results/me/WES/mapQ/%s_FFex2_cut10",SAMPLE[6]);
fileFFPE6 = sprintf("/project/results/me/WES/mapQ/%s_FFPEex2_cut10",SAMPLE[6]);
tableFF6=read.table(fileFF6,sep="\t",header=F,fill=T);
tableFFPE6=read.table(fileFFPE6,sep="\t",header=F,fill=T);

fileFF7 = sprintf("/project/results/me/WES/mapQ/%s_FFex2_cut10",SAMPLE[7]);
fileFFPE7 = sprintf("/project/results/me/WES/mapQ/%s_FFPEex2_cut10",SAMPLE[7]);
tableFF7=read.table(fileFF7,sep="\t",header=F,fill=T);
tableFFPE7=read.table(fileFFPE7,sep="\t",header=F,fill=T);

fileFF8 = sprintf("/project/results/me/WES/mapQ/%s_FFex2_cut10",SAMPLE[8]);
fileFFPE8 = sprintf("/project/results/me/WES/mapQ/%s_FFPEex2_cut10",SAMPLE[8]);
tableFF8=read.table(fileFF8,sep="\t",header=F,fill=T);
tableFFPE8=read.table(fileFFPE8,sep="\t",header=F,fill=T);

fileFF9 = sprintf("/project/results/me/WES/mapQ/%s_FFex2_cut10",SAMPLE[9]);
fileFFPE9 = sprintf("/project/results/me/WES/mapQ/%s_FFPEex2_cut10",SAMPLE[9]);
tableFF9=read.table(fileFF9,sep="\t",header=F,fill=T);
tableFFPE9=read.table(fileFFPE9,sep="\t",header=F,fill=T);

fileFF10 = sprintf("/project/results/me/WES/mapQ/%s_FFex2_cut10",SAMPLE[10]);
fileFFPE10 = sprintf("/project/results/me/WES/mapQ/%s_FFPEex2_cut10",SAMPLE[10]);
tableFF10=read.table(fileFF10,sep="\t",header=F,fill=T);
tableFFPE10=read.table(fileFFPE10,sep="\t",header=F,fill=T);

fileFF11 = sprintf("/project/results/me/WES/mapQ/%s_FFex2_cut10",SAMPLE[11]);
fileFFPE11 = sprintf("/project/results/me/WES/mapQ/%s_FFPEex2_cut10",SAMPLE[11]);
tableFF11=read.table(fileFF11,sep="\t",header=F,fill=T);
tableFFPE11=read.table(fileFFPE11,sep="\t",header=F,fill=T);

fileFF12 = sprintf("/project/results/me/WES/mapQ/%s_FFex2_cut10",SAMPLE[12]);
fileFFPE12 = sprintf("/project/results/me/WES/mapQ/%s_FFPEex2_cut10",SAMPLE[12]);
tableFF12=read.table(fileFF12,sep="\t",header=F,fill=T);
tableFFPE12=read.table(fileFFPE12,sep="\t",header=F,fill=T);

fileFF13 = sprintf("/project/results/me/WES/mapQ/%s_FFex2_cut10",SAMPLE[13]);
fileFFPE13 = sprintf("/project/results/me/WES/mapQ/%s_FFPEex2_cut10",SAMPLE[13]);
tableFF13=read.table(fileFF13,sep="\t",header=F,fill=T);
tableFFPE13=read.table(fileFFPE13,sep="\t",header=F,fill=T);

pdffile = sprintf("/project/results/me/WES/plots/mapQ_boxplot/WES_FFex_FFPEex_boxplot.pdf");
pdf (file = pdffile, width=10, height = 4, pointsize=1)
	par (las = 1,mar=c(5,5,1,1),oma=c(0,2,0,0),col=textcolor)
	boxplot ( tableFF1$V1, tableFFPE1$V1, tableFF2$V1, tableFFPE2$V1,  tableFF3$V1, tableFFPE3$V1, tableFF4$V1, tableFFPE4$V1, tableFF5$V1, tableFFPE5$V1, tableFF6$V1, tableFFPE6$V1, 
	 tableFF7$V1, tableFFPE7$V1,  tableFF8$V1, tableFFPE8$V1,  tableFF9$V1, tableFFPE9$V1,  tableFF10$V1, tableFFPE10$V1,  tableFF11$V1, tableFFPE11$V1, 
	 tableFF12$V1, tableFFPE12$V1,  tableFF13$V1, tableFFPE13$V1, outline=F, xlab="", ylab="", ylim = c(33,62),
	 col=c("#9ACD32","#1E90FF"),fill=c("#9ACD32","#1E90FF"),cex.main=1,cex=1,cex.axis=2,col.lab=textcolor,col.axis=textcolor,adj=1)
	legend("bottomright",c("FF","FFPE"),col=c("#9ACD32","#1E90FF"),fill=c("#9ACD32","#1E90FF"), cex=2.1, text.col=textcolor, bty="n")
	mtext (sprintf("mapping quality"),side=2,line=4,cex=2.5,las=0,col=textcolor,font=2)
	mtext("         2474     2561     2640     2685     2938     3050     3356     4079      4191     (N)14119(T)       (N)22285(T)", side=1,line=2, cex=2.2, col=textcolor,adj=0, font =2)
dev.off()


say=sprintf("%s created!",pdffile);
print(say); 













########
table1b=boxplot(tableFFPE1)
table1b$n=1000
table2b=boxplot(tableFFPE2)
table2b$n=1000
table3b=boxplot(tableFFPE3)
table3b$n=1000
table4b=boxplot(tableFFPE4)
table4b$n=1000
table5b=boxplot(tableFFPE5)
table5b$n=1000
table6b=boxplot(tableFFPE6)
table6b$n=1000
table7b=boxplot(tableFFPE7)
table7b$n=1000
table8b=boxplot(tableFFPE8)
table8b$n=1000
table9b=boxplot(tableFFPE9)
table9b$n=1000
table10b=boxplot(tableFFPE10)
table10b$n=1000
table11b=boxplot(tableFFPE11)
table11b$n=1000
table12b=boxplot(tableFFPE12)
table12b$n=1000
table13b=boxplot(tableFFPE13)
table13b$n=1000

table1a=boxplot(tableFF1)
table1a$n=1000
table2a=boxplot(tableFF2)
table2a$n=1000
table3a=boxplot(tableFF3)
table3a$n=1000
table4a=boxplot(tableFF4)
table4a$n=1000
table5a=boxplot(tableFF5)
table5a$n=1000
table6a=boxplot(tableFF6)
table6a$n=1000
table7a=boxplot(tableFF7)
table7a$n=1000
table8a=boxplot(tableFF8)
table8a$n=1000
table9a=boxplot(tableFF9)
table9a$n=1000
table10a=boxplot(tableFF10)
table10a$n=1000
table11a=boxplot(tableFF11)
table11a$n=1000
table12a=boxplot(tableFF12)
table12a$n=1000
table13a=boxplot(tableFF13)
table13a$n=1000

pdffile = sprintf("/project/results/me/WES/plots/mapQ_boxplot/WES_FFex_FFPEex_boxplot.pdf");
pdf (file = pdffile, width=10, height = 4, pointsize=1)
par (las = 1,mar=c(4,6,1,1),oma=c(0,2,0,0),col=textcolor)
bxp (table1a, table1b, add=T) table2a, table2b, table3a, table3b, table4a, table4b, table5a, table5b, table6a, table6b, table7a, table7b, table8a, table8b, table9a, table9b, table10a, table10b, table11a, table11b, table12a, table12b, table13a, table13b)
box("plot", col=textcolor)  
dev.off()



table.box=boxplot(tableFF1$V1, tableFFPE1$V1, tableFF2$V1, tableFFPE2$V1, tableFF3$V1, tableFFPE3,
 tableFF4$V1, tableFFPE4$V1, tableFF5$V1, tableFFPE5$V1, tableFF6$V1, tableFFPE6$V1, tableFF7,
 tableFFPE7$V1, tableFF8$V1, tableFFPE8$V1, tableFF9$V1, tableFFPE9$V1, tableFF10$V1, tableFFPE10$V1, 
 tableFF11$V1, tableFFPE11$V1, tableFF12$V1, tableFFPE12$V1, tableFF13$V1, tableFFPE13$V1,plot=F)
 
pdffile = sprintf("/project/results/me/WES/plots/mapQ_boxplot/WES_FFex_FFPEex_boxplot.pdf");
pdf (file = pdffile, width=10, height = 4, pointsize=1)
par (las = 1,mar=c(4,6,1,1),oma=c(0,2,0,0),col=textcolor)
bxp (table.box)
dev.off()





test=cbind(table_FF[1],table_FF[2],table_FF[3])

pdffile = sprintf("/project/results/me/WES/plots/mapQ_boxplot/WES_FFex_FFPEex_boxplot.pdf");
pdf (file = pdffile, width=10, height = 4, pointsize=1)
par (las = 1,mar=c(4,6,1,1),oma=c(0,2,0,0),col=textcolor)
boxplot ( test, xlab="         FF                   FFPE",
 ylab="mapping quality",col=c(3,2),cex.main=1,cex.lab=1.8,cex=1,cex.axis=1.5,col.lab=textcolor,col.axis=textcolor,adj=1)
#mtext (sprintf("mapping quality"),side=2,line=2,cex=1,las=0,col=textcolor)
box("plot", col=textcolor)  
dev.off()

