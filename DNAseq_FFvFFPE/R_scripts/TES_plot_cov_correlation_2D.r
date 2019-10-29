targeted = "";
targeted = "targeted_";

textcolor = "#2E2E2E";
SAMPLE = array(0,13)
SAMPLE  =  c("2474","2561","2640","2685","2938","3050","3356","4079","4191","CN","CT","DN","DT")

SAMPLE2=array(0,13)
SAMPLE2[10]="14119(N)";
SAMPLE2[11]="14119(T)";
SAMPLE2[12]="22285(N)";
SAMPLE2[13]="22285(T)";

#i=6
max = 6000;
file = sprintf("/project/results/me/TES/pileup_distribution/%s_FF_FFPE_%sintersect.distribution",SAMPLE[s],targeted);
pdffile = sprintf("/project/results/me/TES/plots/cov_correlation/2D/intersection/TES_%s%s_FF_FFPE_pearson_correlation.pdf",targeted,SAMPLE[s]);
pdffile_fw = sprintf("/project/results/me/TES/plots/cov_correlation/2D/intersection/TES_%s%s_FF_FFPE_forward_pearson_correlation.pdf",targeted,SAMPLE[s]);
pdffile_rv = sprintf("/project/results/me/TES/plots/cov_correlation/2D/intersection/TES_%s%s_FF_FFPE_reverse_pearson_correlation.pdf",targeted,SAMPLE[s]);
table=read.table(file,sep="\t",header=T);

cor=cor((table$A+table$T+table$G+table$C+table$a+table$t+table$g+table$c),(table$A.1+table$T.1+table$G.1+table$C.1+table$a.1+table$t.1+table$g.1+table$c.1),method="pearson")

######## PDFs
pdf (file = pdffile, width=100, height=100, pointsize=100)
par (oma=c(0,0,0,0),mar=c(6,6,1,1),col=textcolor)
if (s<10) {
	plot((table$A+table$T+table$G+table$C+table$a+table$t+table$g+table$c),(table$A.1+table$T.1+table$G.1+table$C.1+table$a.1+table$t.1+table$g.1+table$c.1),
	 pch=3,lwd=10,xlab = sprintf("coverage of %s_FF",SAMPLE[s]),ylab = sprintf("coverage of %s_FFPE",SAMPLE[s]),xlim=c(0,max),ylim=c(0,max),cex.axis=1.8,cex=1.8,cex.lab=1.8,
	 ,col="#999999",col.lab=textcolor,col.axis=textcolor)
} else {
	plot((table$A+table$T+table$G+table$C+table$a+table$t+table$g+table$c),(table$A.1+table$T.1+table$G.1+table$C.1+table$a.1+table$t.1+table$g.1+table$c.1),
	 pch=3,lwd=10,xlab = sprintf("coverage of %s_FF",SAMPLE2[s]),ylab = sprintf("coverage of %s_FFPE",SAMPLE2[s]),xlim=c(0,max),ylim=c(0,max),cex.axis=1.8,cex=1.8,cex.lab=1.8,
	 ,col="#999999",col.lab=textcolor,col.axis=textcolor)
}
mtext (sprintf("correlation coefficient = %.3f",cor),side=1,line=5,cex=1.6,adj=0.5,col=textcolor)
box("plot", col=textcolor) 
dev.off()

say=sprintf("%s created!",pdffile);
print(say); 

#####
cor_fw=cor((table$A+table$T+table$G+table$C),(table$A.1+table$T.1+table$G.1+table$C.1),method="pearson")

pdf (file = pdffile_fw, width=100, height=100, pointsize=100)
par (oma=c(0,0,0,0),mar=c(6,6,1,1),col=textcolor)
if (s<10) {
	plot((table$A+table$T+table$G+table$C),(table$A.1+table$T.1+table$G.1+table$C.1),
	 pch=3,lwd=10,xlab = sprintf("coverage of %s_FF",SAMPLE[s]),ylab = sprintf("coverage of %s_FFPE",SAMPLE[s]),xlim=c(0,max),ylim=c(0,max),cex.axis=1.8,cex=1.8,cex.lab=1.8,
	 ,col="#999999",col.lab=textcolor,col.axis=textcolor)
} else {
	plot((table$A+table$T+table$G+table$C),(table$A.1+table$T.1+table$G.1+table$C.1),
	 pch=3,lwd=10,xlab = sprintf("coverage of %s_FF",SAMPLE2[s]),ylab = sprintf("coverage of %s_FFPE",SAMPLE2[s]),xlim=c(0,max),ylim=c(0,max),cex.axis=1.8,cex=1.8,cex.lab=1.8,
	 ,col="#999999",col.lab=textcolor,col.axis=textcolor)
}
mtext (sprintf("correlation coefficient = %.3f",cor_fw),side=1,line=5,cex=1.6,adj=0.5,col=textcolor)
box("plot", col=textcolor) 
dev.off()

say=sprintf("%s created!",pdffile);
print(say); 

#####
cor_rv=cor((table$a+table$t+table$g+table$c),(table$a.1+table$t.1+table$g.1+table$c.1),method="pearson")

pdf (file = pdffile_rv, width=100, height=100, pointsize=100)
par (oma=c(0,0,0,0),mar=c(6,6,1,1),col=textcolor)
if (s<10) {
	plot((table$a+table$t+table$g+table$c),(table$a.1+table$t.1+table$g.1+table$c.1),
	 pch=3,lwd=10,xlab = sprintf("coverage of %s_FF",SAMPLE[s]),ylab = sprintf("coverage of %s_FFPE",SAMPLE[s]),xlim=c(0,max),ylim=c(0,max),cex.axis=1.8,cex=1.8,cex.lab=1.8,
	 ,col="#999999",col.lab=textcolor,col.axis=textcolor)
} else {
	plot((table$a+table$t+table$g+table$c),(table$a.1+table$t.1+table$g.1+table$c.1),
	 pch=3,lwd=10,xlab = sprintf("coverage of %s_FF",SAMPLE2[s]),ylab = sprintf("coverage of %s_FFPE",SAMPLE2[s]),xlim=c(0,max),ylim=c(0,max),cex.axis=1.8,cex=1.8,cex.lab=1.8,
	 ,col="#999999",col.lab=textcolor,col.axis=textcolor)
}
mtext (sprintf("correlation coefficient = %.3f",cor_rv),side=1,line=5,cex=1.6,adj=0.5,col=textcolor)
box("plot", col=textcolor) 
dev.off()

say=sprintf("%s created!",pdffile);
print(say); 

