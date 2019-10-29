SAMPLE = array(0,13)
SAMPLE  =  c("2474","2561","2640","2685","2938","3050","3356","4079","4191","CN","CT","DN","DT")

SAMPLE2 = array(0,13)
SAMPLE2[10] = "14119(N)";
SAMPLE2[11] = "14119(T)";
SAMPLE2[12] = "22285(N)";
SAMPLE2[13] = "22285(T)";

textcolor = "#2E2E2E";
 
file = sprintf("/home/yhoang/Documents/me/plot_data/TES/pileup_distribution/transition_all.txt");
pretable=read.table(file, sep="\t", header=TRUE, fill=TRUE);

rownames(pretable)=pretable[,1]
table=subset(pretable, select=c(2:25))
group1=table[1:9,]
group2=table[10:13,]

pdffile = sprintf("/home/yhoang/Documents/me/plot_data/TES/plots/conversion/TES_transition_all.pdf");
pdf (file = pdffile, width=18, height = 4, pointsize=1)
	par (las = 1, mfrow = c(1,2), mar=c(5,8,1,5),oma=c(0,5,0,0),col=textcolor)
	
	boxplot ( group1,  xaxt = "n", staplewex = 1, boxwex=0.6, ylim = c( 0.0003,0.0045), #width = c(rep(0.5,24)),
	 col=c("#9ACD32","#1E90FF"), fill=c("#9ACD32","#1E90FF"), cex.main=1, cex=1, cex.axis=3.5, col.lab=textcolor, col.axis=textcolor, adj=1, font.lab=2, cex.lab=2.2)
	legend("topright", c("FF","FFPE"), col=c("#9ACD32","#1E90FF"),fill=c("#9ACD32","#1E90FF"), cex=4, text.col=textcolor, bty="n")
	mtext (sprintf("Group 1"),side=2,line=9,cex=5,las=0,col=textcolor,font=1)
	mtext("     A>T A>G A>C T>A T>G T>C G>A G>T G>C C>A C>T C>G", side=1,line=2, cex=3, col=textcolor,adj=0, font =2)
	
	boxplot ( group2,  xaxt = "n", staplewex = 1, boxwex=0.6, ylim = c( 0.0003,0.0045),
	 col=c("#9ACD32","#1E90FF"), fill=c("#9ACD32","#1E90FF"), cex.main=1, cex=1,cex.axis=3.5, col.lab=textcolor, col.axis=textcolor, adj=1, font.lab=2, cex.lab=2.2)
	mtext (sprintf("Group 2"),side=2,line=9,cex=5,las=0,col=textcolor,font=1)
	mtext("     A>T A>G A>C T>A T>G T>C G>A G>T G>C C>A C>T C>G", side=1,line=2, cex=3, col=textcolor,adj=0, font =2)
dev.off()


say=sprintf("%s created!",pdffile);
print(say); 



