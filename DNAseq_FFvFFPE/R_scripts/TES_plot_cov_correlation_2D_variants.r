print("Set s (1-13) and max (5000)");
max=5000
targeted = "";
#targeted = "targeted_";

threshold = "q43cov13";
#threshold = "pass";

textcolor = "#2E2E2E";
SAMPLE=array(0,13)
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

file_con = sprintf("/project/results/me/TES/SNV_comparison/genotyper/with_intersection/%s_FF_FFPE_%s_%scon",SAMPLE[s],threshold,targeted);
file_dis = sprintf("/project/results/me/TES/SNV_comparison/genotyper/with_intersection/%s_FF_FFPE_%s_%sdis",SAMPLE[s],threshold,targeted);

pdffile = sprintf("/project/results/me/TES/plots/cov_correlation/2D/SNV/%s%s_FF_FFPE_%s_pearson_correlation.pdf",targeted,SAMPLE[s],threshold);

table_con=read.table(file_con,sep="\t",header=T);
table_dis=read.table(file_dis,sep="\t",header=T);
len_con=length(table_con[,20])
len_dis=length(table_dis[,20])
it_FP=0
it_FN=0
it_DIS=0
FP=array(0,c(len_dis,2))
FN=array(0,c(len_dis,2))
DIS=array(0,c(len_dis,2))

corcon=cor(table_con$Depth,table_con$Depth.1,method="pearson")
cordis=cor(table_dis$Depth,table_dis$Depth.1,method="pearson")

pdf (file = pdffile, width=300, height=300, pointsize=300)
par (oma=c(0,0,0,0),mar=c(6,6,1,1),col=textcolor)
if (s<10) {
	plot(table_con$Depth,table_con$Depth.1, type = "p", xlab = sprintf("coverage of %s_FF",SAMPLE[s]), ylab = sprintf("coverage of %s_FFPE",SAMPLE[s]), pch=3,lwd=30,xlim=c(0,max),ylim=c(0,max)
	 ,col="#999999",cex.axis=1.8,cex=1.8,cex.lab=2,col.lab=textcolor,col.axis=textcolor)
} else {
	plot(table_con$Depth,table_con$Depth.1, type = "p", xlab = sprintf("coverage of %s_FF",SAMPLE2[s]), ylab = sprintf("coverage of %s_FFPE",SAMPLE2[s]), pch=3,lwd=30,xlim=c(0,max),ylim=c(0,max)
	 ,col="#999999",cex.axis=1.8,cex=1.8,cex.lab=2,col.lab=textcolor,col.axis=textcolor)
}
par (new=T,col.lab=textcolor)

box("outer",col=textcolor) 
box("figure", col=textcolor)
box("inner",col=textcolor) 

for (i in 1:len_dis){
	if (is.na(table_dis[i,6])==TRUE) { #False Positives
		par (new=T,col.lab=textcolor)
		plot(table_dis[i,7],table_dis[i,18], type = "p", pch=3,xlim=c(0,max),ylim=c(0,max),col="#1E90FF",lwd=30,xaxt="n",yaxt="n",ylab="",xlab="",cex=1.8,col.lab=textcolor)
		it_FP=it_FP+1
	} else if (is.na(table_dis[i,19])==TRUE) { #False Negatives
		par (new=T,col.lab=textcolor)
		plot(table_dis[i,7],table_dis[i,18], type = "p", pch=3,xlim=c(0,max),ylim=c(0,max),col="#9ACD32",lwd=30,xaxt="n",yaxt="n",ylab="",xlab="",cex=1.8,col.lab=textcolor)
		it_FN=it_FN+1
	}
}
for (i in 1:len_dis){
	if ( is.na(table_dis[i,6])==FALSE && (is.na(table_dis[i,19])==FALSE) ){
		par (new=T,col.lab=textcolor)
		plot(table_dis[i,7],table_dis[i,18], type = "p", pch=3,xlim=c(0,max),ylim=c(0,max),col="#FF4500",lwd=30,xaxt="n",yaxt="n",ylab="",xlab="",cex=1.8,col.lab=textcolor)
		it_DIS=it_DIS+1
	}
}

box("plot", col=textcolor)  

legend("topleft",c(sprintf("concordance, n=%s",len_con),sprintf("discordance, n=%s",it_DIS),sprintf("false positives, n=%s",it_FP),sprintf("false negatives, n=%s",it_FN)),
 col=c("#999999","#FF4500","#1E90FF","#9ACD32"),pch=3,lwd=20, lty=NA,cex=1.7,border = textcolor,text.col=textcolor)
mtext (sprintf("coefficient: concordance = %.3f, discordance = %.3f",corcon,cordis),side=1,line=5,cex=1.7,col=textcolor,adj=0.5)

dev.off()

say=sprintf("%s created!",pdffile);
print(say); 
