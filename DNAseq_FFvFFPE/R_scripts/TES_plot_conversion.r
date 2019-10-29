#!/bin/bash

method="basecall";
y_tick = c((0:5)*0.2)
textcolor = "#2E2E2E";


sample="CT";
SNV_FILE = sprintf("~/Dropbox/workspace_me/R_scripts/TES_%s_conversion_%s",sample,method);

pretable=read.table(SNV_FILE,sep="\t",header=TRUE);
rownames(pretable)=pretable[,1]
table=subset(pretable, select=c(2:30))


pdffile = sprintf("/project/results/me/TES/plots/conversion/TES_%s_conversion_%s.pdf",sample,method);
pdf (file = pdffile, width=30, height=20, pointsize=30)  
par (oma=c(0,0,0,0),mar=c(3,5,2,2),col=textcolor)
barplot(cbind(table[,1],table[,2],table[,3],table[,4],table[,5],table[,6],table[,7],table[,8],table[,9],table[,10],table[,11],table[,12],table[,13],table[,14],table[,15],
	table[,16],table[,17],table[,18],table[,19],table[,20],table[,21],table[,22],table[,23],table[,24],table[,25],table[,26],table[,27],table[,28],table[,29]), yaxt="n",
	 col=c("#00376c","#0055a7","#0073e2","#0a86ff","#45a3ff","#80c0ff"),
	 space = 0, border=NA)
axis(2,lty=1,lwd=7,at=(y_tick),labels=y_tick,las=0, cex.axis=2,col=textcolor,col.lab=textcolor,col.axis=textcolor)
legend("topright",rownames(table), cex=1.4,border="dark grey",text.col=textcolor,bg = "white", col=c("#00376c","#0055a7","#0073e2","#0a86ff","#45a3ff","#80c0ff"), fill=c("#00376c","#0055a7","#0073e2","#0a86ff","#45a3ff","#80c0ff"))
mtext("       2474   2561  2640   2685  2938   3050  3356   4079  4191  (N)14119(T) (N)22285(T)", side=1,line=1, cex=1.4, col=textcolor,adj=0, font =2)
dev.off()

## "#80c0ff","#45a3ff","#0a86ff","#0073e2","#0055a7","#00376c"  blue light > dark
## "#00376c","#0055a7","#0073e2","#0a86ff","#45a3ff","#80c0ff"  blue dark > light
## "#FF0000","#FF6600","#FFCC00","#ccff00","#66ff00","#00FF00"  red yellow green
## "#c9e491","#b1d961","#9acd32","#769e26","#536e1b","#2f3f0f"	green light > dark
## "#2f3f0f","#536e1b","#769e26","#9acd32","#b1d961","#c9e491"	green dark > light

test1=cbind(table[,1],table[,2],table[,3],table[,4],table[,5],table[,6],table[,7],table[,8],table[,9],table[,10],table[,11],table[,12],table[,13],table[,14],table[,15])
test2=cbind(table[,16],table[,17],table[,18],table[,19],table[,20],table[,21],table[,22],table[,23],table[,24],table[,25],table[,26],table[,27],table[,28],table[,29])
pdffile = sprintf("/project/results/me/TES/plots/conversion/TESt_%s_conversion_%s.pdf",sample,method);
pdf (file = pdffile, width=30, height=20, pointsize=30)  
par (oma=c(0,0,0,0),mar=c(3,5,2,2),col=textcolor)
barplot(test1,test2)
dev.off()


#################

sample="AG";
SNV_FILE = sprintf("~/Dropbox/workspace_me/R_scripts/TES_conversion_%s_%s",sample,method);
pdffile = sprintf("/project/results/me/TES/plots/conversion/TES_conversion_%s_%s.pdf",sample,method);

pretable=read.table(SNV_FILE,sep="\t",header=TRUE);
rownames(pretable)=pretable[,1]
table=subset(pretable, select=c(2:30))


pdffile = sprintf("/project/results/me/TES/plots/conversion/TES_%s_conversion_%s.pdf",sample,method);
pdf (file = pdffile, width=30, height=20, pointsize=30)  
par (oma=c(0,0,0,0),mar=c(3,5,2,2),col=textcolor)
barplot(cbind(table[,1],table[,2],table[,3],table[,4],table[,5],table[,6],table[,7],table[,8],table[,9],table[,10],table[,11],table[,12],table[,13],table[,14],table[,15],
	table[,16],table[,17],table[,18],table[,19],table[,20],table[,21],table[,22],table[,23],table[,24],table[,25],table[,26],table[,27],table[,28],table[,29]), yaxt="n",
	col=c("#c9e491","#b1d961","#9acd32","#769e26","#536e1b","#2f3f0f"), space = 0, border=NA)
axis(2,lty=1,lwd=7,at=(y_tick),labels=y_tick,las=0, cex.axis=2,col=textcolor,col.lab=textcolor,col.axis=textcolor)
legend("topright",rownames(table), cex=1.4,border="dark grey",text.col=textcolor,bg = "white", 
	col=c("#c9e491","#b1d961","#9acd32","#769e26","#536e1b","#2f3f0f"), fill=c("#c9e491","#b1d961","#9acd32","#769e26","#536e1b","#2f3f0f"))
mtext("       2474   2561  2640   2685  2938   3050  3356   4079  4191  (N)14119(T) (N)22285(T)", side=1,line=1, cex=1.4, col=textcolor,adj=0, font =2)
dev.off()


say=sprintf("%s created!",pdffile); 
print(say);   


