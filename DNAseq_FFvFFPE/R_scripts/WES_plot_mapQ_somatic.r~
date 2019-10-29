#print("Set s (10 or 13)");
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

file1 = sprintf("/project/results/me/WES/SNV_comparison/genotyper/with_intersection/%s/%sN_%sT_FFex_%s_exonic_dis",targeted,SAMPLE[s],SAMPLE[s],threshold);
file2 = sprintf("/project/results/me/WES/SNV_comparison/genotyper/with_intersection/%s/%sN_%sT_FFPEex_%s_exonic_dis",targeted,SAMPLE[s],SAMPLE[s],threshold);

table1=read.table(file1,sep="\t",header=T,fill=TRUE);
table2=read.table(file2,sep="\t",header=T,fill=TRUE);

pdffile = sprintf("/project/results/me/WES/plots/mapQ_correlation/%s_WES_%s_FFex_FFPEex_%s_mapQ.pdf",targeted,threshold,SAMPLE[s],threshold);

minQ = 100
maxQ = 0
# FFcordances
it_FF=0
len_FF=length(table1[,6])
stepFF=array(0,c(len_FF,101))
for (i in 1:len_FF){
        it_FF=it_FF+1
	stepFF[i,1]=table1[i,6]
	diff=table1[i,19]-table1[i,6]
	for (j in 2:101){
		stepFF[i,j]=stepFF[i,j-1]+diff/100
	}
}

# FFPEcordances
it_FFPE=0
len_FFPE=length(table2[,18])
stepFFPE=array(0,c(len_FFPE,101))
for (i in 1:len_FFPE){
	it_FFPE=it_FFPE+1
	stepFFPE[i,1]=round(table2[i,6])
	diff=table2[i,19]-table2[i,6]
	for (j in 2:101){
		stepFFPE[i,j]=stepFFPE[i,j-1]+diff/100
	}
}

x=seq(25,60, by = 2)
xlines=c(0:100)

it_DIS = 0
it_FP = 0
it_FN = 0

pdf (pdffile, width=300, height=350, pointsize=300)
par (mar=c(2,4,2,4)+0.1,oma=c(0,3,0,3),col=textcolor)
plot(0,0,ylab="",xlab="",ylim=c(25,60),xlim=c(4,96),type="n",xaxt="n",yaxt="n")
# draw an axis on the left and right
axis(2, at=x,lwd=7,labels=x, las=2, cex.axis=2,col=textcolor,col.lab=textcolor,col.axis=textcolor)
axis(4, at=x,lwd=7,labels=x, las=2, cex.axis=2,col=textcolor,col.lab=textcolor,col.axis=textcolor)
# add a title for the right and left axis
if (s<10) {
	mtext(sprintf("mapping quality of %s_TUMOR",SAMPLE[s]), side=4, line=5, cex=2.2,col=textcolor,font=2)
	mtext(sprintf("mapping quality of %s_NORMAL",SAMPLE[s]), side=2, line=5, cex=2.2,col=textcolor,font=2)
} else {
	mtext(sprintf("mapping quality of %s_TUMOR",SAMPLE2[s]), side=4, line=5, cex=2.2,col=textcolor,font=2)
	mtext(sprintf("mapping quality of %s_NORMAL",SAMPLE2[s]), side=2, line=5, cex=2.2,col=textcolor,font=2)
}
for (i in 1:len_FFPE){
	lines(xlines,stepFFPE[i,],col="#9ACD32",lwd=60)
}
for (i in 1:len_FF){
	lines(xlines,stepFF[i,],col="#1E90FF",lwd=40)
}

legend ("bottom", c(paste("FF, n=",len_FF),paste("FFPE, n=",len_FFPE)),
 col=c("#1E90FF","#9ACD32"),lty=c(1,1),lwd=(55),cex=1.9,text.col=textcolor,border = textcolor,bg = "white")
dev.off()


say=sprintf("%s created!",pdffile);
print(say);   

