# creating /home/yhoang/Documents/me/plot_data/TES/plots/conversion/TES_transition_2motifs_main.pdf

SAMPLE = array(0,13)
SAMPLE  =  c("2474","2561","2640","2685","2938","3050","3356","4079","4191","CN","CT","DN","DT")

SAMPLE2 = array(0,13)
SAMPLE2[10] = "14119(N)";
SAMPLE2[11] = "14119(T)";
SAMPLE2[12] = "22285(N)";
SAMPLE2[13] = "22285(T)";

textcolor = "#2E2E2E";

##### IMPORT
# IMPORT C Transitions subsequent
fileCsub = "/home/yhoang/Documents/me/plot_data/TES/pileup_distribution/C_transition_2motifs_sub.txt";
pretableCsub=read.table(fileCsub, sep="\t", header=TRUE, fill=TRUE);
rownames(pretableCsub)=pretableCsub[,1]
tableCsub=subset(pretableCsub, select=c(2:9))
groupC1sub=tableCsub[1:9,]
groupC2sub=tableCsub[10:13,]

# IMPORT C Transitions antecedent
fileCant = "/home/yhoang/Documents/me/plot_data/TES/pileup_distribution/C_transition_2motifs_ant.txt";
pretableCant=read.table(fileCant, sep="\t", header=TRUE, fill=TRUE);
rownames(pretableCant)=pretableCant[,1]
tableCant=subset(pretableCant, select=c(2:9))
groupC1ant=tableCant[1:9,]
groupC2ant=tableCant[10:13,]

# IMPORT G Transitions subsequent
fileGsub = "/home/yhoang/Documents/me/plot_data/TES/pileup_distribution/G_transition_2motifs_sub.txt";
pretableGsub=read.table(fileGsub, sep="\t", header=TRUE, fill=TRUE);
rownames(pretableGsub)=pretableGsub[,1]
tableGsub=subset(pretableGsub, select=c(2:9))
groupG1sub=tableGsub[1:9,]
groupG2sub=tableGsub[10:13,]

# IMPORT G Transitions antecedent
fileGant = "/home/yhoang/Documents/me/plot_data/TES/pileup_distribution/G_transition_2motifs_ant.txt";
pretableGant=read.table(fileGant, sep="\t", header=TRUE, fill=TRUE);
rownames(pretableGant)=pretableGant[,1]
tableGant=subset(pretableGant, select=c(2:9))
groupG1ant=tableGant[1:9,]
groupG2ant=tableGant[10:13,]

##### MAIN FILE
pdffile = sprintf("/home/yhoang/Documents/me/plot_data/TES/plots/conversion/TES_transition_2motifs_main.pdf",target);
pdf (file = pdffile, width=20, height = 4, pointsize=1)
	par (las = 1, mfrow = c(1,4), mar=c(5,9,1,5),oma=c(0,0,0,0),col=textcolor)
	
### Group 1 C>T subsequent
	boxplot ( groupC1sub,  xaxt = "n", yaxt = "n", staplewex = 1, boxwex=0.4, ylim = c( 0.0008,0.0076),
	 col=c("#9ACD32","#1E90FF"), fill=c("#9ACD32","#1E90FF"), cex.main=1, cex=1, cex.axis=5, col.lab=textcolor, col.axis=textcolor, adj=1, font.lab=2, cex.lab=9,ylab="Group 1")
	legend("topleft", c("FF","FFPE"), col=c("#9ACD32","#1E90FF"),fill=c("#9ACD32","#1E90FF"), cex=5, text.col=textcolor, bty="n")
	mtext("   CA>TA  CT>TT CG>TG  CC>TC", side=1,line=2.5, cex=3, col=textcolor,adj=0, font =2)
	text (8,0.0072,"A",cex=12, font =2)

### Group 1 G>A antecedent
	boxplot ( groupG1ant,  xaxt = "n", staplewex = 1, boxwex=0.4, ylim = c( 0.0008,0.0076),
	 col=c("#9ACD32","#1E90FF"), fill=c("#9ACD32","#1E90FF"), cex.main=1, cex=1, cex.axis=5, col.lab=textcolor, col.axis=textcolor, adj=1, font.lab=2, cex.lab=2.2)
	mtext("  AG>AA  TG>TA  GG>GA CG>CA", side=1,line=2.5, cex=3, col=textcolor,adj=0, font =2)
	text (8,0.0072,"B",cex=12, font =2)
		
### Group 2 C>T subsequent
	boxplot ( groupC2sub,  xaxt = "n", yaxt = "n", staplewex = 1, boxwex=0.4, ylim = c( 0.0008,0.0076),
	 col=c("#9ACD32","#1E90FF"), fill=c("#9ACD32","#1E90FF"), cex.main=1, cex=1,cex.axis=3.5, col.lab=textcolor, col.axis=textcolor, adj=1, font.lab=2, cex.lab=9,ylab="Group 2")
	mtext("   CA>TA  CT>TT CG>TG  CC>TC", side=1,line=2.5, cex=3, col=textcolor,adj=0, font =2)
	text (8,0.0072,"C",cex=12, font =2)
		
### Group 2 G>A antecedent
	boxplot ( groupG2ant,  xaxt = "n", staplewex = 1, boxwex=0.4, ylim = c( 0.0008,0.0076),
	 col=c("#9ACD32","#1E90FF"), fill=c("#9ACD32","#1E90FF"), cex.main=1, cex=1,cex.axis=5, col.lab=textcolor, col.axis=textcolor, adj=1, font.lab=2, cex.lab=2.2)
	mtext("  AG>AA  TG>TA  GG>GA CG>CA", side=1,line=2.5, cex=3, col=textcolor,adj=0, font =2)	
	text (8,0.0072,"D",cex=12, font =2)

dev.off()
	
say=sprintf("%s created!",pdffile);
print(say);


