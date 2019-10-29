# creating /home/yhoang/Documents/me/plot_data/TES/plots/conversion/TES_WES_transition_2motifs_ungrouped_main.pdf

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

# IMPORT C Transitions antecedent
fileCant = "/home/yhoang/Documents/me/plot_data/TES/pileup_distribution/C_transition_2motifs_ant.txt";
pretableCant=read.table(fileCant, sep="\t", header=TRUE, fill=TRUE);
rownames(pretableCant)=pretableCant[,1]
tableCant=subset(pretableCant, select=c(2:9))

# IMPORT G Transitions subsequent
fileGsub = "/home/yhoang/Documents/me/plot_data/TES/pileup_distribution/G_transition_2motifs_sub.txt";
pretableGsub=read.table(fileGsub, sep="\t", header=TRUE, fill=TRUE);
rownames(pretableGsub)=pretableGsub[,1]
tableGsub=subset(pretableGsub, select=c(2:9))

# IMPORT G Transitions antecedent
fileGant = "/home/yhoang/Documents/me/plot_data/TES/pileup_distribution/G_transition_2motifs_ant.txt";
pretableGant=read.table(fileGant, sep="\t", header=TRUE, fill=TRUE);
rownames(pretableGant)=pretableGant[,1]
tableGant=subset(pretableGant, select=c(2:9))


##### MAIN FILE
pdffile = sprintf("/home/yhoang/Documents/me/plot_data/TES/plots/conversion/TES_transition_2motifs_ungrouped_main.pdf",target);
pdf (file = pdffile, width=10, height = 4, pointsize=1)
	par (las = 1, mfrow = c(1,2), mar=c(5,9,1,5),oma=c(0,0,0,0),col=textcolor)
	
### C>T subsequent
	boxplot ( tableCsub,  xaxt = "n", yaxt = "n", staplewex = 1, boxwex=0.4, ylim = c( 0.0008,0.0076),
	 col=c("#9ACD32","#1E90FF"), fill=c("#9ACD32","#1E90FF"), cex.main=1, cex=1, cex.axis=3, col.lab=textcolor, col.axis=textcolor, adj=1, font.lab=2, cex.lab=7,ylab="TES")
	legend("topleft", c("FF","FFPE"), col=c("#9ACD32","#1E90FF"),fill=c("#9ACD32","#1E90FF"), cex=4, text.col=textcolor, bty="n")
	mtext(" CA>TA CT>TT CG>TG CC>TC", side=1,line=2.5, cex=3, col=textcolor,adj=0, font =2)
	text (8,0.0072,"A",cex=8, font =2)

### G>A antecedent
	boxplot ( tableGant,  xaxt = "n", staplewex = 1, boxwex=0.4, ylim = c( 0.0008,0.0076),
	 col=c("#9ACD32","#1E90FF"), fill=c("#9ACD32","#1E90FF"), cex.main=1, cex=1, cex.axis=4, col.lab=textcolor, col.axis=textcolor, adj=1, font.lab=2, cex.lab=2.2)
	mtext("AG>AA TG>TA GG>GA CG>CA", side=1,line=2.5, cex=3, col=textcolor,adj=0, font =2)
	text (8,0.0072,"B",cex=8, font =2)
		

dev.off()
	
say=sprintf("%s created!",pdffile);
print(say);


