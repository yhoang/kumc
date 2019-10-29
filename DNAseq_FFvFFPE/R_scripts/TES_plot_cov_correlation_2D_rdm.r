INTERSECTION = c(7240035,6769069,6756966,7101083,7257881,6625630,6922354,6694114,6837727,5693503,5183853,3125273,2709650)
INTERSECTION = c(1565296,1563561,1562043,1564751,1566659,1564287,1563359,1561726,1564989,1558394,1550765,1483348,1406885)

targeted = "";
targeted = "targeted_";
max = 4000;
random = 2000;
textcolor = "#2E2E2E";
SAMPLE=array(0,13)
SAMPLE = c("2474","2561","2640","2685","2938","3050","3356","4079","4191","CN","CT","DN","DT")
SAMPLE2=array(0,13)
SAMPLE2[10]="14119(N)";
SAMPLE2[11]="14119(T)";
SAMPLE2[12]="22285(N)";
SAMPLE2[13]="22285(T)";

file = sprintf("/project/results/me/TES/pileup_distribution/%s_FF_FFPE_%sintersect.distribution",SAMPLE[s],targeted);
pdffile = sprintf("/project/results/me/TES/plots/cov_correlation/2D/intersection/%s/TES_%s%s_FF_FFPE_pearson_correlation.pdf",random,targeted,SAMPLE[s]);
pdffile_fw = sprintf("/project/results/me/TES/plots/cov_correlation/2D/intersection/%s/TES_%s%s_FF_FFPE_forward_pearson_correlation.pdf",random,targeted,SAMPLE[s]);
pdffile_rv = sprintf("/project/results/me/TES/plots/cov_correlation/2D/intersection/%s/TES_%s%s_FF_FFPE_reverse_pearson_correlation.pdf",random,targeted,SAMPLE[s]);
table=read.table(file,sep="\t",header=T);

rdm = sample(1:INTERSECTION[s],random)
table_rdm=rbind(table[rdm[1],],table[rdm[2],],table[rdm[3],],table[rdm[4],],table[rdm[5],])
for (k in seq(6,length(rdm),5)){
	table_rdm=rbind(table_rdm,table[rdm[k],],table[rdm[k+1],],table[rdm[k+2],],table[rdm[k+3],],table[rdm[k+4],])
}

cor=cor((table_rdm$A+table_rdm$T+table_rdm$G+table_rdm$C+table_rdm$a+table_rdm$t+table_rdm$g+table_rdm$c),(table_rdm$A.1+table_rdm$T.1+table_rdm$G.1+table_rdm$C.1+table_rdm$a.1+table_rdm$t.1+table_rdm$g.1+table_rdm$c.1),method="pearson",use="na.or.complete")
y1_tick=(0:4)*1000
x1_tick=(0:4)*1000
######## PDFs
pdf (file = pdffile, width=100, height=100, pointsize=100)
par (oma=c(0,0,0,0),mar=c(6,6,1,1),col=textcolor)
if (s<10) {
	plot((table_rdm$A+table_rdm$T+table_rdm$G+table_rdm$C+table_rdm$a+table_rdm$t+table_rdm$g+table_rdm$c),(table_rdm$A.1+table_rdm$T.1+table_rdm$G.1+table_rdm$C.1+table_rdm$a.1+table_rdm$t.1+table_rdm$g.1+table_rdm$c.1),
	 pch=3,lwd=10,xlab = sprintf("coverage of %s_FF",SAMPLE[s]),ylab = sprintf("coverage of %s_FFPE",SAMPLE[s]),xlim=c(0,max),ylim=c(0,max),cex.lab=2.2,font.lab=2,cex=2,
	 ,col="#999999",col.lab=textcolor,yaxt="n",xaxt="n")
	axis(1,lty=1,lwd=15,at=(x1_tick),labels=x1_tick,las=0, cex.axis=2,col=textcolor,col.lab=textcolor,col.axis=textcolor)
	axis(2,lty=1,lwd=15,at=(y1_tick),labels=y1_tick,las=0, cex.axis=2,col=textcolor,col.lab=textcolor,col.axis=textcolor)
} else {
	plot((table_rdm$A+table_rdm$T+table_rdm$G+table_rdm$C+table_rdm$a+table_rdm$t+table_rdm$g+table_rdm$c),(table_rdm$A.1+table_rdm$T.1+table_rdm$G.1+table_rdm$C.1+table_rdm$a.1+table_rdm$t.1+table_rdm$g.1+table_rdm$c.1),
	 pch=3,lwd=10,xlab = sprintf("coverage of %s_FF",SAMPLE2[s]),ylab = sprintf("coverage of %s_FFPE",SAMPLE2[s]),xlim=c(0,max),ylim=c(0,max),cex.lab=2.2,font.lab=2,cex=2,
	 ,col="#999999",col.lab=textcolor,yaxt="n",xaxt="n")
	axis(1,lty=1,lwd=15,at=(x1_tick),labels=x1_tick,las=0, cex.axis=2,col=textcolor,col.lab=textcolor,col.axis=textcolor)
	axis(2,lty=1,lwd=15,at=(y1_tick),labels=y1_tick,las=0, cex.axis=2,col=textcolor,col.lab=textcolor,col.axis=textcolor)
}
mtext (sprintf("correlation coefficient = %.3f",cor),side=1,line=5,cex=2,adj=0.5,col=textcolor)
box("plot", col=textcolor) 
dev.off()

say=sprintf("%s created!",pdffile);
print(say); 

#####
cor_fw=cor((table_rdm$A+table_rdm$T+table_rdm$G+table_rdm$C),(table_rdm$A.1+table_rdm$T.1+table_rdm$G.1+table_rdm$C.1),method="pearson",use="na.or.complete")

pdf (file = pdffile_fw, width=100, height=100, pointsize=100)
par (oma=c(0,0,0,0),mar=c(6,6,1,1),col=textcolor)
if (s<10) {
	plot((table_rdm$A+table_rdm$T+table_rdm$G+table_rdm$C),(table_rdm$A.1+table_rdm$T.1+table_rdm$G.1+table_rdm$C.1),
	 pch=3,lwd=10,xlab = sprintf("coverage of %s_FF (forward strand)",SAMPLE[s]),ylab = sprintf("coverage of %s_FFPE (forward strand)",SAMPLE[s]),xlim=c(0,max),ylim=c(0,max),cex.lab=2.2,font.lab=2,cex=2,
	 ,col="#999999",col.lab=textcolor,yaxt="n",xaxt="n")
	axis(1,lty=1,lwd=15,at=(x1_tick),labels=x1_tick,las=0, cex.axis=2,col=textcolor,col.lab=textcolor,col.axis=textcolor)
	axis(2,lty=1,lwd=15,at=(y1_tick),labels=y1_tick,las=0, cex.axis=2,col=textcolor,col.lab=textcolor,col.axis=textcolor)
} else {
	plot((table_rdm$A+table_rdm$T+table_rdm$G+table_rdm$C),(table_rdm$A.1+table_rdm$T.1+table_rdm$G.1+table_rdm$C.1),
	 pch=3,lwd=10,xlab = sprintf("coverage of %s_FF (forward strand)",SAMPLE2[s]),ylab = sprintf("coverage of %s_FFPE (forward strand)",SAMPLE2[s]),xlim=c(0,max),ylim=c(0,max),cex.lab=2.2,font.lab=2,cex=2,
	 ,col="#999999",col.lab=textcolor,yaxt="n",xaxt="n")
	axis(1,lty=1,lwd=15,at=(x1_tick),labels=x1_tick,las=0, cex.axis=2,col=textcolor,col.lab=textcolor,col.axis=textcolor)
	axis(2,lty=1,lwd=15,at=(y1_tick),labels=y1_tick,las=0, cex.axis=2,col=textcolor,col.lab=textcolor,col.axis=textcolor)
}
mtext (sprintf("correlation coefficient = %.3f",cor_fw),side=1,line=5,cex=2,adj=0.5,col=textcolor)
box("plot", col=textcolor) 
dev.off()

say=sprintf("%s created!",pdffile);
print(say); 

#####
cor_rv=cor((table_rdm$a+table_rdm$t+table_rdm$g+table_rdm$c),(table_rdm$a.1+table_rdm$t.1+table_rdm$g.1+table_rdm$c.1),method="pearson",use="na.or.complete")

pdf (file = pdffile_rv, width=100, height=100, pointsize=100)
par (oma=c(0,0,0,0),mar=c(6,6,1,1),col=textcolor)
if (s<10) {
	plot((table_rdm$a+table_rdm$t+table_rdm$g+table_rdm$c),(table_rdm$a.1+table_rdm$t.1+table_rdm$g.1+table_rdm$c.1),
	 pch=3,lwd=10,xlab = sprintf("coverage of %s_FF (reverse strand)",SAMPLE[s]),ylab = sprintf("coverage of %s_FFPE (reverse strand)",SAMPLE[s]),xlim=c(0,max),ylim=c(0,max),cex.lab=2.2,font.lab=2,cex=2,
	 ,col="#999999",col.lab=textcolor,yaxt="n",xaxt="n")
	axis(1,lty=1,lwd=15,at=(x1_tick),labels=x1_tick,las=0, cex.axis=2,col=textcolor,col.lab=textcolor,col.axis=textcolor)
	axis(2,lty=1,lwd=15,at=(y1_tick),labels=y1_tick,las=0, cex.axis=2,col=textcolor,col.lab=textcolor,col.axis=textcolor)
} else {
	plot((table_rdm$a+table_rdm$t+table_rdm$g+table_rdm$c),(table_rdm$a.1+table_rdm$t.1+table_rdm$g.1+table_rdm$c.1),
	 pch=3,lwd=10,xlab = sprintf("coverage of %s_FF (reverse strand)",SAMPLE2[s]),ylab = sprintf("coverage of %s_FFPE (reverse strand)",SAMPLE2[s]),xlim=c(0,max),ylim=c(0,max),cex.lab=2.2,font.lab=2,cex=2,
	 ,col="#999999",col.lab=textcolor,yaxt="n",xaxt="n")
	axis(1,lty=1,lwd=15,at=(x1_tick),labels=x1_tick,las=0, cex.axis=2,col=textcolor,col.lab=textcolor,col.axis=textcolor)
	axis(2,lty=1,lwd=15,at=(y1_tick),labels=y1_tick,las=0, cex.axis=2,col=textcolor,col.lab=textcolor,col.axis=textcolor)
}
mtext (sprintf("correlation coefficient = %.3f",cor_rv),side=1,line=5,cex=2,adj=0.5,col=textcolor)
box("plot", col=textcolor) 
dev.off()

say=sprintf("%s created!",pdffile);
print(say); 

