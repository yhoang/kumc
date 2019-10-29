SAMPLE = array(0,13)
SAMPLE  =  c("2474","2561","2640","2685","2938","3050","3356","4079","4191","CN","CT","DN","DT")

SAMPLE2 = array(0,13)
SAMPLE2[10] = "14119(N)";
SAMPLE2[11] = "14119(T)";
SAMPLE2[12] = "22285(N)";
SAMPLE2[13] = "22285(T)";

textcolor = "#2E2E2E";


################## C>T
target ="C";
 
file2 = sprintf("/project/results/me/TES/pileup_distribution/%s_transition_2motifs_ant.txt",target);
pretable_ant=read.table(file2, sep="\t", header=TRUE, fill=TRUE);

rownames(pretable_ant)=pretable_ant[,1]
table_ant=subset(pretable_ant, select=c(2:9))
group1_ant=table_ant[1:9,]
group2_ant=table_ant[10:13,]

file1 = sprintf("/project/results/me/TES/pileup_distribution/%s_transition_2motifs_sub.txt",target);
pretable_sub=read.table(file1, sep="\t", header=TRUE, fill=TRUE);
rownames(pretable_sub)=pretable_sub[,1]
table_sub=subset(pretable_sub, select=c(2:9))
group1_sub=table_sub[1:9,]
group2_sub=table_sub[10:13,]


###### Group1

pdffile = sprintf("/project/results/me/TES/plots/conversion/TES_%s_transition_2motifs_group1.pdf",target);
pdf (file = pdffile, width=10, height = 4, pointsize=1)
	par (las = 1, mfrow = c(1,2), mar=c(5,8,1,2),oma=c(0,2,0,0),col=textcolor)
	boxplot ( group1_ant,  xaxt = "n", staplewex = 0.5, boxwex=0.4, ylim = c( 0.0008,0.0076),
	 col=c("#9ACD32","#1E90FF"), fill=c("#9ACD32","#1E90FF"), cex.main=1, cex=1, cex.axis=2, col.lab=textcolor, col.axis=textcolor, adj=1, font.lab=2, cex.lab=2.2)
	legend("top", c("FF","FFPE"), col=c("#9ACD32","#1E90FF"),fill=c("#9ACD32","#1E90FF"), cex=2.1, text.col=textcolor, bty="n")
	mtext (sprintf("Group 1"),side=2,line=7,cex=2.5,las=0,col=textcolor,font=2)
	mtext("      AC>AT       TC>TT      GC>GT      CC>CT", side=1,line=2, cex=2.2, col=textcolor,adj=0, font =2)
	
	boxplot ( group1_sub,  xaxt = "n", staplewex = 0.5, boxwex=0.4, ylim = c( 0.0008,0.0076),
	 col=c("#9ACD32","#1E90FF"), fill=c("#9ACD32","#1E90FF"), cex.main=1, cex=1, cex.axis=2, col.lab=textcolor, col.axis=textcolor, adj=1, font.lab=2, cex.lab=2.2)
	 legend("top", c("FF","FFPE"), col=c("#9ACD32","#1E90FF"),fill=c("#9ACD32","#1E90FF"), cex=2.1, text.col=textcolor, bty="n")
	 mtext (sprintf("Group 1"),side=2,line=7,cex=2.5,las=0,col=textcolor,font=2)
	 mtext("      CA>TA      CT>TT      CG>TG      CC>TC", side=1,line=2, cex=2.2, col=textcolor,adj=0, font =2)

dev.off()

say=sprintf("%s created!",pdffile);
print(say); 

###### Group2

pdffile = sprintf("/project/results/me/TES/plots/conversion/TES_%s_transition_2motifs_group2.pdf",target);
pdf (file = pdffile, width=10, height = 4, pointsize=1)
	par (las = 1, mfrow = c(1,2), mar=c(5,8,1,2),oma=c(0,2,0,0),col=textcolor)

	boxplot ( group2_ant,  xaxt = "n", staplewex = 0.5, boxwex=0.4, ylim = c( 0.0008,0.0076),
	 col=c("#9ACD32","#1E90FF"), fill=c("#9ACD32","#1E90FF"), cex.main=1, cex=1,cex.axis=2, col.lab=textcolor, col.axis=textcolor, adj=1, font.lab=2, cex.lab=2.2)
	 legend("top", c("FF","FFPE"), col=c("#9ACD32","#1E90FF"),fill=c("#9ACD32","#1E90FF"), cex=2.1, text.col=textcolor, bty="n")
	 mtext (sprintf("Group 2"),side=2,line=7,cex=2.5,las=0,col=textcolor,font=2)
	 mtext("      AC>AT       TC>TT      GC>GT      CC>CT", side=1,line=2, cex=2.2, col=textcolor,adj=0, font =2)

	boxplot ( group2_sub,  xaxt = "n", staplewex = 0.5, boxwex=0.4, ylim = c( 0.0008,0.0076),
	 col=c("#9ACD32","#1E90FF"), fill=c("#9ACD32","#1E90FF"), cex.main=1, cex=1,cex.axis=2, col.lab=textcolor, col.axis=textcolor, adj=1, font.lab=2, cex.lab=2.2)
	legend("top", c("FF","FFPE"), col=c("#9ACD32","#1E90FF"),fill=c("#9ACD32","#1E90FF"), cex=2.1, text.col=textcolor, bty="n")
	mtext (sprintf("Group 2"),side=2,line=7,cex=2.5,las=0,col=textcolor,font=2)
	mtext("      CA>TA      CT>TT      CG>TG      CC>TC", side=1,line=2, cex=2.2, col=textcolor,adj=0, font =2)
dev.off()


say=sprintf("%s created!",pdffile);
print(say);



###################### G>A
target = "G";

file2 = sprintf("/project/results/me/TES/pileup_distribution/%s_transition_2motifs_ant.txt",target);
pretable_ant=read.table(file2, sep="\t", header=TRUE, fill=TRUE);

rownames(pretable_ant)=pretable_ant[,1]
table_ant=subset(pretable_ant, select=c(2:9))
group1_ant=table_ant[1:9,]
group2_ant=table_ant[10:13,]

file1 = sprintf("/project/results/me/TES/pileup_distribution/%s_transition_2motifs_sub.txt",target);
pretable_sub=read.table(file1, sep="\t", header=TRUE, fill=TRUE);

rownames(pretable_sub)=pretable_sub[,1]
table_sub=subset(pretable_sub, select=c(2:9))
group1_sub=table_sub[1:9,]
group2_sub=table_sub[10:13,]


###### Group 1

pdffile = sprintf("/project/results/me/TES/plots/conversion/TES_%s_transition_2motifs_group1.pdf",target);
pdf (file = pdffile, width=10, height = 4, pointsize=1)
	par (las = 1, mfrow = c(1,2), mar=c(5,8,1,2),oma=c(0,2,0,0),col=textcolor)
	
	boxplot ( group1_ant,  xaxt = "n", staplewex = 0.5, boxwex=0.4, ylim = c( 0.0008,0.0076),
	 col=c("#9ACD32","#1E90FF"), fill=c("#9ACD32","#1E90FF"), cex.main=1, cex=1, cex.axis=2, col.lab=textcolor, col.axis=textcolor, adj=1, font.lab=2, cex.lab=2.2)
	 legend("top", c("FF","FFPE"), col=c("#9ACD32","#1E90FF"),fill=c("#9ACD32","#1E90FF"), cex=2.1, text.col=textcolor, bty="n")
	 mtext (sprintf("Group 1"),side=2,line=7,cex=2.5,las=0,col=textcolor,font=2)
	 mtext("     AG>AA       TG>TA     GG>GA     CG>CA", side=1,line=2, cex=2.2, col=textcolor,adj=0, font =2)
	
	boxplot ( group1_sub,  xaxt = "n", staplewex = 0.5, boxwex=0.4, ylim = c( 0.0008,0.0076),
	 col=c("#9ACD32","#1E90FF"), fill=c("#9ACD32","#1E90FF"), cex.main=1, cex=1, cex.axis=2, col.lab=textcolor, col.axis=textcolor, adj=1, font.lab=2, cex.lab=2.2)
	 legend("top", c("FF","FFPE"), col=c("#9ACD32","#1E90FF"),fill=c("#9ACD32","#1E90FF"), cex=2.1, text.col=textcolor, bty="n")
	 mtext (sprintf("Group 1"),side=2,line=7,cex=2.5,las=0,col=textcolor,font=2)
	 mtext("      GA>AA     GT>AT      GG>AG      GC>AC", side=1,line=2, cex=2.2, col=textcolor,adj=0, font =2)
	
	
dev.off()


say=sprintf("%s created!",pdffile);
print(say); 


###### Group 2

pdffile = sprintf("/project/results/me/TES/plots/conversion/TES_%s_transition_2motifs_group2.pdf",target);
pdf (file = pdffile, width=10, height = 4, pointsize=1)
	par (las = 1, mfrow = c(1,2), mar=c(5,8,1,2),oma=c(0,2,0,0),col=textcolor)
	
	boxplot ( group2_ant,  xaxt = "n", staplewex = 0.5, boxwex=0.4, ylim = c( 0.0008,0.0076),
	 col=c("#9ACD32","#1E90FF"), fill=c("#9ACD32","#1E90FF"), cex.main=1, cex=1,cex.axis=2, col.lab=textcolor, col.axis=textcolor, adj=1, font.lab=2, cex.lab=2.2)
	 legend("top", c("FF","FFPE"), col=c("#9ACD32","#1E90FF"),fill=c("#9ACD32","#1E90FF"), cex=2.1, text.col=textcolor, bty="n")
	 mtext (sprintf("Group 2"),side=2,line=7,cex=2.5,las=0,col=textcolor,font=2)
	 mtext("     AG>AA       TG>TA     GG>GA     CG>CA", side=1,line=2, cex=2.2, col=textcolor,adj=0, font =2)
	
	boxplot ( group2_sub,  xaxt = "n", staplewex = 0.5, boxwex=0.4, ylim = c( 0.0008,0.0076),
	 col=c("#9ACD32","#1E90FF"), fill=c("#9ACD32","#1E90FF"), cex.main=1, cex=1,cex.axis=2, col.lab=textcolor, col.axis=textcolor, adj=1, font.lab=2, cex.lab=2.2)
	 legend("top", c("FF","FFPE"), col=c("#9ACD32","#1E90FF"),fill=c("#9ACD32","#1E90FF"), cex=2.1, text.col=textcolor, bty="n")
	 mtext (sprintf("Group 2"),side=2,line=7,cex=2.5,las=0,col=textcolor,font=2)
	 mtext("      GA>AA     GT>AT      GG>AG      GC>AC", side=1,line=2, cex=2.2, col=textcolor,adj=0, font =2)
dev.off()

say=sprintf("%s created!",pdffile);
print(say); 

