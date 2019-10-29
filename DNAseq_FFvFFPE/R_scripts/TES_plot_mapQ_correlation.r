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



file1 = sprintf("/project/results/me/TES/SNV_comparison/genotyper/with_intersection/%s_FF_FFPE_%s_%scon",SAMPLE[s],threshold,targeted);
file2 = sprintf("/project/results/me/TES/SNV_comparison/genotyper/with_intersection/%s_FF_FFPE_%s_%sdis",SAMPLE[s],threshold,targeted);

table1=read.table(file1,sep="\t",header=T,fill=TRUE);
table2=read.table(file2,sep="\t",header=T,fill=TRUE);

pdffile = sprintf("/project/results/me/TES/plots/mapQ_correlation/q43_or_cov13/%s%s_FF_FFPE_%s_mapQ.pdf",targeted,SAMPLE[s],threshold);

minQ = 100
maxQ = 0
# concordances
it_con=0
len_con=length(table1[,9])
stepcon=array(0,c(len_con,101))
for (i in 1:len_con){
        it_con=it_con+1
	stepcon[i,1]=table1[i,9]
	diff=table1[i,20]-table1[i,9]
	for (j in 2:101){
		stepcon[i,j]=stepcon[i,j-1]+diff/100
	}
}

# discordances
it_dis=0
len_dis=length(table2[,20])
stepdis=array(0,c(len_dis,101))
for (i in 1:len_dis){
	it_dis=it_dis+1
	stepdis[i,1]=round(table2[i,9])
	diff=round(table2[i,20])-round(table2[i,9])
	for (j in 2:101){
		stepdis[i,j]=stepdis[i,j-1]+diff/100
	}
}

x=seq(28,64, by = 2)
xlines=c(0:100)

it_dis1 = 0
it_dis2 = 0
it_dis3 = 0

pdf (pdffile, width=300, height=350, pointsize=300)
par (mar=c(2,4,2,4)+0.1,oma=c(0,2,0,2),col=textcolor)
plot(0,0,ylab="",xlab="",ylim=c(28,64),xlim=c(4,96),type="n",xaxt="n",yaxt="n")
# draw an axis on the left and right
axis(2, at=x,labels=x, las=2, cex.axis=1.5,col=textcolor,col.lab=textcolor,col.axis=textcolor)
axis(4, at=x,labels=x, las=2, cex.axis=1.5,col=textcolor,col.lab=textcolor,col.axis=textcolor)
# add a title for the right and left axis
if (s<10) {
	mtext(sprintf("mapping quality of %s_FFPE",SAMPLE[s]), side=4, line=4, cex=1.7,col=textcolor)
	mtext(sprintf("mapping quality of %s_FF",SAMPLE[s]), side=2, line=4, cex=1.7,col=textcolor)
} else {
	mtext(sprintf("mapping quality of %s_FFPE",SAMPLE2[s]), side=4, line=4, cex=1.7,col=textcolor)
	mtext(sprintf("mapping quality of %s_FF",SAMPLE2[s]), side=2, line=4, cex=1.7,col=textcolor)
}
# add concordance lines
for (i in 1:len_con){
	lines(xlines,stepcon[i,],col="#999999",lwd=40)
}
for (i in 1:len_dis){
	if (stepdis[i,1]<43 || stepdis[i,101]<43) { 
		lines(xlines,stepdis[i,],col="#1E90FF",lwd=50)
		it_dis1=it_dis1+1
	} else if ( (table2[i,7]<13) || (table2[i,18]<13)){ 
		lines(xlines,stepdis[i,],col="#9ACD32",lwd=50)
		it_dis2=it_dis2+1
	} else {
		lines(xlines,stepdis[i,],col="#FF4500",lwd=50)
		it_dis3=it_dis3+1
	}
	
}
box("plot", col=textcolor)  
legend ("bottom", c(paste("concordance, n=",it_con),paste("discordance with mapQ<43, n=",it_dis1),paste("discordance with cov < 13, n=",it_dis2),paste("discordance, n=",it_dis3)),
 col=c("#999999","#1E90FF","#9ACD32","#FF4500"),lty=c(1,1),lwd=50,cex=1.5,text.col=textcolor,border = textcolor,bg = "white")
dev.off()

say=sprintf("%s created!",pdffile);
print(say);   

