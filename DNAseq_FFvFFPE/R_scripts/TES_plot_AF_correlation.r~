#!/bin/bash
#for (i in 1:13) {
#print("Set s (1-13)");
SAMPLE=array(0,13)

targeted = "";
targeted = "targeted_";
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
SAMPLE[10]="CN";
SAMPLE[11]="CT";
SAMPLE[12]="DN";
SAMPLE[13]="DT";

SAMPLE2=array(0,13)
SAMPLE2[10]="14119(N)";
SAMPLE2[11]="14119(T)";
SAMPLE2[12]="22285(N)";
SAMPLE2[13]="22285(T)";

file_con = sprintf("/project/results/me/TES/plots/alleleFreq/genotyper/correlation_data/AF_%s_FF_FFPE_%s_%scon_correlation",SAMPLE[s],threshold,targeted);
file_dis = sprintf("/project/results/me/TES/plots/alleleFreq/genotyper/correlation_data/AF_%s_FF_FFPE_%s_%sdis_correlation",SAMPLE[s],threshold,targeted);
pdffile =  sprintf("/project/results/me/TES/plots/alleleFreq/genotyper/correlation_plot_1D/%s/TES_%sAF1D_%s_FF_FFPE_%s_condis_correlation.pdf",threshold,targeted,SAMPLE[s],threshold);

table_con=read.table(file_con,sep="\t",header=T,fill=F);
table_dis=read.table(file_dis,sep="\t",header=T,fill=F);


minQ = 100
maxQ = 0
# concordances
len_con=length(table_con[,4])
stepcon=array(0,c(len_con,101))
for (i in 1:len_con){
	stepcon[i,1]=table_con[i,4]
	diff=table_con[i,6]-table_con[i,4]
	for (j in 2:101){
		stepcon[i,j]=stepcon[i,j-1]+diff/100
	}
}

# discordances
len_dis=length(table_dis[,6])
stepdis=array(0,c(len_dis,101))
for (i in 1:len_dis){
	stepdis[i,1]=table_dis[i,4]
	diff=table_dis[i,6]-table_dis[i,4]
	for (j in 2:101){
		stepdis[i,j]=stepdis[i,j-1]+diff/100
	}
}

x=seq(0,1, by = 0.2)
xlines=c(0:100)

it_CON = 0
it_FP = 0
it_FN = 0
it_DIS = 0

pdf (pdffile, width=300, height=350, pointsize=300)
par (mar=c(2,4,2,4)+0.1,oma=c(0,3,0,3),col=textcolor)
plot(0,0,ylab="",xlab="",ylim=c(0,1),xlim=c(4,96),type="n",xaxt="n",yaxt="n")
# draw an axis on the left and right

axis(2, at=x,lwd=7,labels=x, las=2, cex.axis=2,col=textcolor,col.lab=textcolor,col.axis=textcolor)
axis(4, at=x,lwd=7,labels=x, las=2, cex.axis=2,col=textcolor,col.lab=textcolor,col.axis=textcolor)
# add a title for the right and left axis
if (s<10){
	mtext(sprintf("allele frequency of %s_FFPE",SAMPLE[s]), side=4, line=5, cex=2.2,col=textcolor,font=2)
	mtext(sprintf("allele frequency of %s_FF",SAMPLE[s]), side=2, line=5, cex=2.2,col=textcolor,font=2)
} else {
	mtext(sprintf("allele frequency of %s_FFPE",SAMPLE2[s]), side=4, line=5, cex=2.2,col=textcolor,font=2)
	mtext(sprintf("allele frequency of %s_FF",SAMPLE2[s]), side=2, line=5, cex=2.2,col=textcolor,font=2)
}
# add lines
for (i in 1:len_con){
	lines(xlines,stepcon[i,],col="#999999",lwd=50)
	it_CON=it_CON+1
}
for (i in 1:len_dis){
	if (table_dis[i,3]=="") { # False Positive
		lines(xlines,stepdis[i,],col="#1E90FF",lwd=50)
		it_FP=it_FP+1
	} else if (table_dis[i,5]=="") { # False Negative
		lines(xlines,stepdis[i,],col="#9ACD32",lwd=50)
		it_FN=it_FN+1
	} 
}
for (i in 1:len_dis){
	if (table_dis[i,3]!="" && table_dis[i,5]!="") {
		lines(xlines,stepdis[i,],col="#FF4500",lwd=55)
		it_DIS=it_DIS+1
	}
}
it_CON=1338
it_FP=7
legend ("bottom", c(paste("concordance, n=",it_CON),paste("discordant calls, n=",it_DIS),paste("false positives, n=",it_FP),paste("false negatives, n=",it_FN)), 
 col=c("#999999","#FF4500","#1E90FF","#9ACD32"),lty=c(1,1),lwd=50,cex=1.9,text.col=textcolor,border = textcolor,bg = "white")
dev.off()


say=sprintf("%s created!",pdffile);
print(say);   

