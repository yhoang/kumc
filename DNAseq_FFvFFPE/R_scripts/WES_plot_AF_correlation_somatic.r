#print("Set s (1-13)");
SAMPLE=array(0,13)
targeted = "somatic";
#targeted = "exonic";
threshold = "q43cov13";
#threshold = "pass";

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
SAMPLE[10]="C";
SAMPLE[11]="C";
SAMPLE[12]="D";
SAMPLE[13]="D";


SAMPLE2=array(0,13)
SAMPLE2[10]="14119";
SAMPLE2[11]="14119";
SAMPLE2[12]="22285";
SAMPLE2[13]="22285";

s=10

file1 = sprintf("/project/results/me/WES/SNV_comparison/genotyper/with_intersection/%s/%sN_%sT_FFex_%s_exonic_dis",targeted,SAMPLE[s],SAMPLE[s],threshold);
file2 = sprintf("/project/results/me/WES/SNV_comparison/genotyper/with_intersection/%s/%sN_%sT_FFPEex_%s_exonic_dis",targeted,SAMPLE[s],SAMPLE[s],threshold);

table_FF=read.table(file1,sep="\t",header=T,fill=TRUE);
table_FFPE=read.table(file2,sep="\t",header=T,fill=TRUE);

pdffile = sprintf("/project/results/me/WES/plots/mapQ_correlation/%s_WES_%s_FFex_FFPEex_%s_AF1.pdf",targeted,threshold,SAMPLE[s],threshold);

len_FF=length(table_FF[,4])
len_FFPE=length(table_FFPE[,1])
it_FP=0
it_FN=0
it_FFPE=0
corFF=cor(table_FF[,8],table_FF[,21],method="pearson")
corFFPE=cor(table_FFPE[,8],table_FFPE[,21],method="pearson")

#maxFFplot = 2000
#plot_FF = table_FF[0:2000,]

y1_tick=(0:10)*0.2
x1_tick=(0:10)*0.2
pdf (file = pdffile, width=300, height=300, pointsize=300)
par (oma=c(0,0,0,0),mar=c(6,6,1,1),col=textcolor)

if (s <10) {
	plot(table_FF[,8],table_FF[,21], type = "p", xlab = sprintf("allele frequency of %s_NORMAL",SAMPLE[s]), ylab = sprintf("allele frequency of %s_TUMOR",SAMPLE[s]), pch=3,lwd=30,xlim=c(0,1),ylim=c(0,1)
	 ,col="#9ACD32",cex.axis=2,cex.lab=2.2,cex=2,col.lab=textcolor,yaxt="n",xaxt="n")
	axis(1,lty=1,lwd=7,at=(x1_tick),labels=x1_tick,las=0, cex.axis=2,col=textcolor,col.lab=textcolor,col.axis=textcolor)
	axis(2,lty=1,lwd=7,at=(y1_tick),labels=y1_tick,las=0, cex.axis=2,col=textcolor,col.lab=textcolor,col.axis=textcolor)
} else {
	plot(table_FF[,8],table_FF[,21], type = "p", xlab = sprintf("allele frequency of %s_NORMAL",SAMPLE2[s]), ylab = sprintf("allele frequency of %s_TUMOR",SAMPLE2[s]), pch=3,lwd=30,xlim=c(0,1),ylim=c(0,1)
	 ,col="#9ACD32",cex.axis=2,cex.lab=2.2,cex=2,col.lab=textcolor,yaxt="n",xaxt="n")
	axis(1,lty=1,lwd=7,at=(x1_tick),labels=x1_tick,las=0, cex.axis=2,col=textcolor,col.lab=textcolor,col.axis=textcolor)
	axis(2,lty=1,lwd=7,at=(y1_tick),labels=y1_tick,las=0, cex.axis=2,col=textcolor,col.lab=textcolor,col.axis=textcolor)
}

par (new=T,col=textcolor)
plot(table_FFPE[,4],table_FFPE[,6], type = "p", pch=3,xlim=c(0,1),ylim=c(0,1),col="#FF4500",lwd=30,xaxt="n",yaxt="n",ylab="",xlab="",cex=2,col.lab=textcolor)

box("plot", col=textcolor)  
legend("topleft",c(sprintf("FF, n=%s",len_FF),sprintf("FFPE, n=%s",len_FFPE)),
 col=c("#1E90FF","#9ACD32"),pch=c(3,3),lwd=c(30), lty=c(NA,NA),cex=2,text.col=textcolor,border = textcolor,bg = "white")
mtext (sprintf("coefficient: FF = %.3f, FFPE = %.3f",corFF,corFFPE),side=1,line=5,cex=1.7,adj=0.5)
dev.off()

say=sprintf("%s created!",pdffile);
print(say); 

