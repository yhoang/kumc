print("Set s (1-9)");
SAMPLE=array(0,9)

targeted = "";
#targeted = "targeted_";

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
pdffile = sprintf("/project/results/me/TES/plots/alleleFreq/genotyper/correlation_plot_2D/%s/TES_%sAF2D_%s_FF_FFPE_%s_condis_correlation.pdf",threshold,targeted,SAMPLE[s],threshold);

table_con=read.table(file_con,sep="\t",header=T,fill=T);
table_dis=read.table(file_dis,sep="\t",header=T,fill=T);

len_con=length(table_con[,4])
len_dis=length(table_dis[,1])
it_FP=0
it_FN=0
it_DIS=0
corcon=cor(table_con[,4],table_con[,6],method="pearson")
cordis=cor(table_dis[,4],table_dis[,6],method="pearson")
y1_tick=(0:10)*0.2
x1_tick=(0:10)*0.2

pdf (file = pdffile, width=300, height=300, pointsize=300)
par (oma=c(0,0,0,0),mar=c(6,6,1,1),col=textcolor)
if (s<10){
	plot(table_con[,4],table_con[,6], type = "p", xlab = sprintf("allele frequency of %s_FF",SAMPLE[s]), ylab = sprintf("allele frequency of %s_FFPE",SAMPLE[s]), font.lab=2,
	 pch=3,lwd=30,xlim=c(0,1),ylim=c(0,1),col="#999999",cex.lab=2.2,cex=2,col.lab=textcolor,yaxt="n",xaxt="n")
	axis(1,lty=1,lwd=7,at=(x1_tick),labels=x1_tick,las=0, cex.axis=2,col=textcolor,col.lab=textcolor,col.axis=textcolor)
	axis(2,lty=1,lwd=7,at=(y1_tick),labels=y1_tick,las=0, cex.axis=2,col=textcolor,col.lab=textcolor,col.axis=textcolor)
} else {
	plot(table_con[,4],table_con[,6], type = "p", xlab = sprintf("allele frequency of %s_FF",SAMPLE2[s]), ylab = sprintf("allele frequency of %s_FFPE",SAMPLE2[s]),font.lab=2,
	 pch=3,lwd=30,xlim=c(0,1),ylim=c(0,1),col="#999999",cex.lab=2.2,cex=2,col.lab=textcolor,yaxt="n",xaxt="n")
	axis(1,lty=1,lwd=7,at=(x1_tick),labels=x1_tick,las=0, cex.axis=2,col=textcolor,col.lab=textcolor,col.axis=textcolor)
	axis(2,lty=1,lwd=7,at=(y1_tick),labels=y1_tick,las=0, cex.axis=2,col=textcolor,col.lab=textcolor,col.axis=textcolor)
}
for (i in 1:len_dis){
	if (table_dis[i,3]=="") { # False Positive
		par (new=T,col=textcolor)
		plot(table_dis[i,4],table_dis[i,6], type = "p", pch=3,xlim=c(0,1),ylim=c(0,1),col="#1E90FF",lwd=30,xaxt="n",yaxt="n",ylab="",xlab="",cex=2,col.lab=textcolor)
		it_FP=it_FP+1
	} else if (table_dis[i,5]=="") { # False Negative
		par (new=T,col=textcolor)
		plot(table_dis[i,4],table_dis[i,6], type = "p", pch=3,xlim=c(0,1),ylim=c(0,1),col="#9ACD32",lwd=30,xaxt="n",yaxt="n",ylab="",xlab="",cex=2,col.lab=textcolor)
		it_FN=it_FN+1
	} 
}

for (i in 1:len_dis){
	if (table_dis[i,3]!="" && table_dis[i,5]!="") {
		par (new=T,col=textcolor)
		plot(table_dis[i,4],table_dis[i,6], type = "p", pch=3,xlim=c(0,1),ylim=c(0,1),col="#FF4500",lwd=30,xaxt="n",yaxt="n",ylab="",xlab="",cex=2,col.lab=textcolor)
		it_DIS=it_DIS+1
	}
}
len_con=5026
it_FP=12
box("plot", col=textcolor)  
legend("topleft",c(sprintf("concordance, n=%s",len_con),sprintf("discordance, n=%s",it_DIS),sprintf("false positives, n=%s",it_FP),sprintf("false negatives, n=%s",it_FN)),
 col=c("#999999","#FF4500","#1E90FF","#9ACD32"),pch=c(3,3),lwd=c(30), lty=c(NA,NA),cex=2,text.col=textcolor,border = textcolor,bg = "white")
mtext (sprintf("coefficient: concordance = %.3f, discordance = %.3f",corcon,cordis),side=1,line=5,cex=1.7,adj=0.5)
dev.off()

say=sprintf("%s created!",pdffile);
print(say); 

