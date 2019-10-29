print("Set s (1-13) and max (5000) and max=1100 for targeted");
max=5000
targeted = "total";
#targeted = "exonic";
targeted = "targeted";
threshold = "q43cov13";
threshold = "pass";

SAMPLE = array(0,13)
SAMPLE  =  c("2474","2561","2640","2685","2938","3050","3356","4079","4191","CN","CT","DN","DT")

SAMPLE2=array(0,13)
SAMPLE2[10]="14119(N)";
SAMPLE2[11]="14119(T)";
SAMPLE2[12]="22285(N)";
SAMPLE2[13]="22285(T)";

textcolor = "#2E2E2E";

file_con = sprintf("/project/results/me/WES/SNV_comparison/genotyper/with_intersection/stats/%s/%s_FFex_FFPEex_%s_%s_con.tsv",targeted,SAMPLE[s],threshold,targeted);
file_dis = sprintf("/project/results/me/WES/SNV_comparison/genotyper/with_intersection/stats/%s/%s_FFex_FFPEex_%s_%s_dis.tsv",targeted,SAMPLE[s],threshold,targeted);
#file_con = sprintf("/project/results/me/WES/SNV_comparison/genotyper/with_intersection/stats/%s/%s_FFex_FFPEex_%s_con.tsv",targeted,SAMPLE[s],threshold);
#file_dis = sprintf("/project/results/me/WES/SNV_comparison/genotyper/with_intersection/stats/%s/%s_FFex_FFPEex_%s_dis.tsv",targeted,SAMPLE[s],threshold);


pdffile = sprintf("/project/results/me/WES/plots/cov_correlation/2D/SNV_FR/%s/WES_%s_%s_FF_FFPE_%s_pearson_correlation.pdf",threshold,targeted,SAMPLE[s],threshold);
pdffile_fw = sprintf("/project/results/me/WES/plots/cov_correlation/2D/SNV_FR/%s/WES_%s_%s_FF_FFPE_%s_forward_pearson_correlation.pdf",threshold,targeted,SAMPLE[s],threshold);
pdffile_rv = sprintf("/project/results/me/WES/plots/cov_correlation/2D/SNV_FR/%s/WES_%s_%s_FF_FFPE_%s_reverse_pearson_correlation.pdf",threshold,targeted,SAMPLE[s],threshold);

table_con=read.table(file_con,sep="\t",header=T,fill=T,na.strings="");
table_dis=read.table(file_dis,sep="\t",header=T,fill=T,na.strings="");
len_con=length(table_con[,2])
len_dis=length(table_dis[,2]) 
it_FP=0
it_FN=0
it_DIS=0

y1_tick=(0:5)*200;
x1_tick=(0:5)*200;

corcon=cor((table_con$A_FF+table_con$T_FF+table_con$G_FF+table_con$C_FF+table_con$a_FF+table_con$t_FF+table_con$g_FF+table_con$c_FF),
 (table_con$A_FFPE+table_con$T_FFPE+table_con$G_FFPE+table_con$C_FFPE+table_con$a_FFPE+table_con$t_FFPE+table_con$g_FFPE+table_con$c_FFPE),method="pearson")
cordis=cor((table_dis$A_FF+table_dis$T_FF+table_dis$G_FF+table_dis$C_FF+table_dis$a_FF+table_dis$t_FF+table_dis$g_FF+table_dis$c_FF),
 (table_dis$A_FFPE+table_dis$T_FFPE+table_dis$G_FFPE+table_dis$C_FFPE+table_dis$a_FFPE+table_dis$t_FFPE+table_dis$g_FFPE+table_dis$c_FFPE),method="pearson")

pdf (file = pdffile, width=300, height=300, pointsize=300)
par (oma=c(0,0,0,0),mar=c(6,6,1,1),col=textcolor)
if (s<10) {
	plot( (table_con$A_FF+table_con$T_FF+table_con$G_FF+table_con$C_FF+table_con$a_FF+table_con$t_FF+table_con$g_FF+table_con$c_FF),
 	 (table_con$A_FFPE+table_con$T_FFPE+table_con$G_FFPE+table_con$C_FFPE+table_con$a_FFPE+table_con$t_FFPE+table_con$g_FFPE+table_con$c_FFPE), 
 	 type = "p", xlab = sprintf("coverage of %s_FF",SAMPLE[s]), ylab = sprintf("coverage of %s_FFPE",SAMPLE[s]), pch=3,lwd=30,xlim=c(0,max),ylim=c(0,max),
	 col="#999999",font.lab=2,cex=2,cex.lab=2.2,col.lab=textcolor,yaxt="n",xaxt="n")
	axis(1,lty=1,lwd=15,at=(x1_tick),labels=x1_tick,las=0, cex.axis=2,col=textcolor,col.lab=textcolor,col.axis=textcolor)
	axis(2,lty=1,lwd=15,at=(y1_tick),labels=y1_tick,las=0, cex.axis=2,col=textcolor,col.lab=textcolor,col.axis=textcolor)
} else {
	plot( (table_con$A_FF+table_con$T_FF+table_con$G_FF+table_con$C_FF+table_con$a_FF+table_con$t_FF+table_con$g_FF+table_con$c_FF),
 	 (table_con$A_FFPE+table_con$T_FFPE+table_con$G_FFPE+table_con$C_FFPE+table_con$a_FFPE+table_con$t_FFPE+table_con$g_FFPE+table_con$c_FFPE), 
 	 type = "p", xlab = sprintf("coverage of %s_FF",SAMPLE2[s]), ylab = sprintf("coverage of %s_FFPE",SAMPLE2[s]), pch=3,lwd=30,xlim=c(0,max),ylim=c(0,max),
	 col="#999999",font.lab=2,cex=2,cex.lab=2.2,col.lab=textcolor,yaxt="n",xaxt="n")
	axis(1,lty=1,lwd=15,at=(x1_tick),labels=x1_tick,las=0, cex.axis=2,col=textcolor,col.lab=textcolor,col.axis=textcolor)
	axis(2,lty=1,lwd=15,at=(y1_tick),labels=y1_tick,las=0, cex.axis=2,col=textcolor,col.lab=textcolor,col.axis=textcolor)
}


for (i in 1:len_dis){
	if (is.na(table_dis[i,4])==TRUE) { #False Positives
		par (new=T,col.lab=textcolor)
		plot( (table_dis[i,13]+table_dis[i,14]+table_dis[i,15]+table_dis[i,16]+table_dis[i,17]+table_dis[i,18]+table_dis[i,19]+table_dis[i,20]),
		 (table_dis[i,22]+table_dis[i,23]+table_dis[i,24]+table_dis[i,25]+table_dis[i,26]+table_dis[i,27]+table_dis[i,28]+table_dis[i,29]),
		 type = "p", pch=3,xlim=c(0,max),ylim=c(0,max),col="#1E90FF",lwd=30,xaxt="n",yaxt="n",ylab="",xlab="",cex=1.8,col.lab=textcolor)
		it_FP=it_FP+1
	} else if (is.na(table_dis[i,11])==TRUE) { #False Negatives
		par (new=T,col.lab=textcolor)
		plot( (table_dis[i,13]+table_dis[i,14]+table_dis[i,15]+table_dis[i,16]+table_dis[i,17]+table_dis[i,18]+table_dis[i,19]+table_dis[i,20]),
		 (table_dis[i,22]+table_dis[i,23]+table_dis[i,24]+table_dis[i,25]+table_dis[i,26]+table_dis[i,27]+table_dis[i,28]+table_dis[i,29]),
		 type = "p", pch=3,xlim=c(0,max),ylim=c(0,max),col="#9ACD32",lwd=30,xaxt="n",yaxt="n",ylab="",xlab="",cex=1.8,col.lab=textcolor)
		it_FN=it_FN+1
	}
}
for (i in 1:len_dis){
	if ( is.na(table_dis[i,4])==FALSE && (is.na(table_dis[i,11])==FALSE) ){
		par (new=T,col.lab=textcolor)
		plot( (table_dis[i,13]+table_dis[i,14]+table_dis[i,15]+table_dis[i,16]+table_dis[i,17]+table_dis[i,18]+table_dis[i,19]+table_dis[i,20]),
		 (table_dis[i,22]+table_dis[i,23]+table_dis[i,24]+table_dis[i,25]+table_dis[i,26]+table_dis[i,27]+table_dis[i,28]+table_dis[i,29]),
		 type = "p", pch=3,xlim=c(0,max),ylim=c(0,max),col="#FF4500",lwd=30,xaxt="n",yaxt="n",ylab="",xlab="",cex=1.8,col.lab=textcolor)
		it_DIS=it_DIS+1
	}
}

box("plot", col=textcolor)  
legend("topleft",c(sprintf("concordance, n=%s",len_con),sprintf("discordance, n=%s",it_DIS),sprintf("false positives, n=%s",it_FP),sprintf("false negatives, n=%s",it_FN)),
 col=c("#999999","#FF4500","#1E90FF","#9ACD32"),pch=3,lwd=20, lty=NA,cex=1.9,border = textcolor,text.col=textcolor)
mtext (sprintf("coefficient: concordance = %.3f, discordance = %.3f",corcon,cordis),side=1,line=5,cex=2,col=textcolor,adj=0.5)
dev.off()

say=sprintf("%s created!",pdffile);
print(say); 

######### forward

corcon_fw=cor((table_con$A_FF+table_con$T_FF+table_con$G_FF+table_con$C_FF),
 (table_con$A_FFPE+table_con$T_FFPE+table_con$G_FFPE+table_con$C_FFPE),method="pearson")
cordis_fw=cor((table_dis$A_FF+table_dis$T_FF+table_dis$G_FF+table_dis$C_FF),
 (table_dis$A_FFPE+table_dis$T_FFPE+table_dis$G_FFPE+table_dis$C_FFPE),method="pearson")

pdf (file = pdffile_fw, width=300, height=300, pointsize=300)
par (oma=c(0,0,0,0),mar=c(6,6,1,1),col=textcolor)
if (s<10) {
	plot( (table_con$A_FF+table_con$T_FF+table_con$G_FF+table_con$C_FF),
 	 (table_con$A_FFPE+table_con$T_FFPE+table_con$G_FFPE+table_con$C_FFPE), 
 	 type = "p", xlab = sprintf("coverage of %s_FF (forward strand)",SAMPLE[s]), ylab = sprintf("coverage of %s_FFPE (forward strand)",SAMPLE[s]), pch=3,lwd=30,xlim=c(0,max),ylim=c(0,max),
	 col="#999999",font.lab=2,cex=2,cex.lab=2.2,col.lab=textcolor,yaxt="n",xaxt="n")
	axis(1,lty=1,lwd=15,at=(x1_tick),labels=x1_tick,las=0, cex.axis=2,col=textcolor,col.lab=textcolor,col.axis=textcolor)
	axis(2,lty=1,lwd=15,at=(y1_tick),labels=y1_tick,las=0, cex.axis=2,col=textcolor,col.lab=textcolor,col.axis=textcolor)
} else {
	plot( (table_con$A_FF+table_con$T_FF+table_con$G_FF+table_con$C_FF),
 	 (table_con$A_FFPE+table_con$T_FFPE+table_con$G_FFPE+table_con$C_FFPE), 
 	 type = "p", xlab = sprintf("coverage of %s_FF (forward strand)",SAMPLE2[s]), ylab = sprintf("coverage of %s_FFPE (forward strand)",SAMPLE2[s]), pch=3,lwd=30,xlim=c(0,max),ylim=c(0,max),
	 col="#999999",font.lab=2,cex=2,cex.lab=2.2,col.lab=textcolor,yaxt="n",xaxt="n")
	axis(1,lty=1,lwd=15,at=(x1_tick),labels=x1_tick,las=0, cex.axis=2,col=textcolor,col.lab=textcolor,col.axis=textcolor)
	axis(2,lty=1,lwd=15,at=(y1_tick),labels=y1_tick,las=0, cex.axis=2,col=textcolor,col.lab=textcolor,col.axis=textcolor)
}


for (i in 1:len_dis){
	if (is.na(table_dis[i,4])==TRUE) { #False Positives
		par (new=T,col.lab=textcolor)
		plot( (table_dis[i,13]+table_dis[i,14]+table_dis[i,15]+table_dis[i,16]),
		 (table_dis[i,22]+table_dis[i,23]+table_dis[i,24]+table_dis[i,25]),
		 type = "p", pch=3,xlim=c(0,max),ylim=c(0,max),col="#1E90FF",lwd=30,xaxt="n",yaxt="n",ylab="",xlab="",cex=1.8,col.lab=textcolor)
	} else if (is.na(table_dis[i,11])==TRUE) { #False Negatives
		par (new=T,col.lab=textcolor)
		plot( (table_dis[i,13]+table_dis[i,14]+table_dis[i,15]+table_dis[i,16]),
		 (table_dis[i,22]+table_dis[i,23]+table_dis[i,24]+table_dis[i,25]),
		 type = "p", pch=3,xlim=c(0,max),ylim=c(0,max),col="#9ACD32",lwd=30,xaxt="n",yaxt="n",ylab="",xlab="",cex=1.8,col.lab=textcolor)
	}
}
for (i in 1:len_dis){
	if ( is.na(table_dis[i,4])==FALSE && (is.na(table_dis[i,11])==FALSE) ){
		par (new=T,col.lab=textcolor)
		plot( (table_dis[i,13]+table_dis[i,14]+table_dis[i,15]+table_dis[i,16]),
		 (table_dis[i,22]+table_dis[i,23]+table_dis[i,24]+table_dis[i,25]),
		 type = "p", pch=3,xlim=c(0,max),ylim=c(0,max),col="#FF4500",lwd=30,xaxt="n",yaxt="n",ylab="",xlab="",cex=1.8,col.lab=textcolor)
	}
}

box("plot", col=textcolor)  
legend("topleft",c(sprintf("concordance, n=%s",len_con),sprintf("discordance, n=%s",it_DIS),sprintf("false positives, n=%s",it_FP),sprintf("false negatives, n=%s",it_FN)),
 col=c("#999999","#FF4500","#1E90FF","#9ACD32"),pch=3,lwd=20, lty=NA,cex=1.9,border = textcolor,text.col=textcolor)
mtext (sprintf("coefficient: concordance = %.3f, discordance = %.3f",corcon_fw,cordis_fw),side=1,line=5,cex=2,col=textcolor,adj=0.5)
dev.off()

say=sprintf("%s created!",pdffile_fw);
print(say); 

#### reverse

corcon_rv=cor((table_con$a_FF+table_con$t_FF+table_con$g_FF+table_con$c_FF),
 (table_con$a_FFPE+table_con$t_FFPE+table_con$g_FFPE+table_con$c_FFPE),method="pearson")
cordis_rv=cor((table_dis$a_FF+table_dis$t_FF+table_dis$g_FF+table_dis$c_FF),
 (table_dis$C_FFPE+table_dis$a_FFPE+table_dis$t_FFPE+table_dis$g_FFPE+table_dis$c_FFPE),method="pearson")

pdf (file = pdffile_rv, width=300, height=300, pointsize=300)
par (oma=c(0,0,0,0),mar=c(6,6,1,1),col=textcolor)
if (s<10) {
	plot( (table_con$a_FF+table_con$t_FF+table_con$g_FF+table_con$c_FF),
 	 (table_con$a_FFPE+table_con$t_FFPE+table_con$g_FFPE+table_con$c_FFPE), 
 	 type = "p", xlab = sprintf("coverage of %s_FF (reverse strand)",SAMPLE[s]), ylab = sprintf("coverage of %s_FFPE (reverse strand)",SAMPLE[s]), pch=3,lwd=30,xlim=c(0,max),ylim=c(0,max),
	 col="#999999",font.lab=2,cex=2,cex.lab=2.2,col.lab=textcolor,yaxt="n",xaxt="n")
	axis(1,lty=1,lwd=15,at=(x1_tick),labels=x1_tick,las=0, cex.axis=2,col=textcolor,col.lab=textcolor,col.axis=textcolor)
	axis(2,lty=1,lwd=15,at=(y1_tick),labels=y1_tick,las=0, cex.axis=2,col=textcolor,col.lab=textcolor,col.axis=textcolor)
} else {
	plot( (table_con$a_FF+table_con$t_FF+table_con$g_FF+table_con$c_FF),
 	 (table_con$a_FFPE+table_con$t_FFPE+table_con$g_FFPE+table_con$c_FFPE), 
 	 type = "p", xlab = sprintf("coverage of %s_FF (reverse strand)",SAMPLE2[s]), ylab = sprintf("coverage of %s_FFPE (reverse strand)",SAMPLE2[s]), pch=3,lwd=30,xlim=c(0,max),ylim=c(0,max),
	 col="#999999",font.lab=2,cex=2,cex.lab=2.2,col.lab=textcolor,yaxt="n",xaxt="n")
	axis(1,lty=1,lwd=15,at=(x1_tick),labels=x1_tick,las=0, cex.axis=2,col=textcolor,col.lab=textcolor,col.axis=textcolor)
	axis(2,lty=1,lwd=15,at=(y1_tick),labels=y1_tick,las=0, cex.axis=2,col=textcolor,col.lab=textcolor,col.axis=textcolor)
}


for (i in 1:len_dis){
	if (is.na(table_dis[i,4])==TRUE) { #False Positives
		par (new=T,col.lab=textcolor)
		plot( (table_dis[i,17]+table_dis[i,18]+table_dis[i,19]+table_dis[i,20]),
		 (table_dis[i,26]+table_dis[i,27]+table_dis[i,28]+table_dis[i,29]),
		 type = "p", pch=3,xlim=c(0,max),ylim=c(0,max),col="#1E90FF",lwd=30,xaxt="n",yaxt="n",ylab="",xlab="",cex=1.8,col.lab=textcolor)
	} else if (is.na(table_dis[i,11])==TRUE) { #False Negatives
		par (new=T,col.lab=textcolor)
		plot( (table_dis[i,17]+table_dis[i,18]+table_dis[i,19]+table_dis[i,20]),
		 (table_dis[i,26]+table_dis[i,27]+table_dis[i,28]+table_dis[i,29]),
		 type = "p", pch=3,xlim=c(0,max),ylim=c(0,max),col="#9ACD32",lwd=30,xaxt="n",yaxt="n",ylab="",xlab="",cex=1.8,col.lab=textcolor)
	}
}
for (i in 1:len_dis){
	if ( is.na(table_dis[i,4])==FALSE && (is.na(table_dis[i,11])==FALSE) ){
		par (new=T,col.lab=textcolor)
		plot( (table_dis[i,17]+table_dis[i,18]+table_dis[i,19]+table_dis[i,20]),
		 (table_dis[i,26]+table_dis[i,27]+table_dis[i,28]+table_dis[i,29]),
		 type = "p", pch=3,xlim=c(0,max),ylim=c(0,max),col="#FF4500",lwd=30,xaxt="n",yaxt="n",ylab="",xlab="",cex=1.8,col.lab=textcolor)
	}
}

box("plot", col=textcolor)  
legend("topleft",c(sprintf("concordance, n=%s",len_con),sprintf("discordance, n=%s",it_DIS),sprintf("false positives, n=%s",it_FP),sprintf("false negatives, n=%s",it_FN)),
 col=c("#999999","#FF4500","#1E90FF","#9ACD32"),pch=3,lwd=20, lty=NA,cex=1.9,border = textcolor,text.col=textcolor)
mtext (sprintf("coefficient: concordance = %.3f, discordance = %.3f",corcon_rv,cordis_rv),side=1,line=5,cex=2,col=textcolor,adj=0.5)
dev.off()

say=sprintf("%s created!",pdffile_rv);
print(say); 


