
textcolor = "#2E2E2E";
SAMPLE="14119";


file_FF = "/project/results/me/WGS/CN_FF.coverage.chr.hist_cut";
file_FFPE = "/project/results/me/WGS/CN_FFPE.coverage.chr.hist_cut";
pdffile = "/project/results/me/WGS/plots/WGS_cov_distribution.pdf";

table_FF=read.table(file_FF,colClasses=c("character","integer","integer"));
table_FFPE=read.table(file_FFPE,colClasses=c("character","integer","integer"));

x_FF = 1: max(table_FF[,2]) 
x_FFPE = 1: max(table_FFPE[,2]) 
histogram_FF = array( 0, max(table_FF[,2]) )
histogram_FFPE = array( 0, max(table_FFPE[,2]) )

for (k in 1:length(table_FFPE[,1]) ) {
	histogram_FFPE[table_FFPE[k,2]] = histogram_FFPE[table_FFPE[k,2]] + table_FFPE[k,3];
}
for (k in 1:length(table_FF[,1]) ) {
	histogram_FF[table_FF[k,2]] = histogram_FF[table_FF[k,2]] + table_FF[k,3];
}

histogram_FF_cut = array (histogram_FF[1:170],170)
histogram_FFPE_cut = array (histogram_FFPE[1:170],170)
x_cut = 1:170
x_tick = c((0:6)*30)
y_tick = c((0:6)*30000000)

pdf (file = pdffile, width=320, height=230, pointsize=300)
	par (oma=c(0,0,0,0),mar=c(7,7,4,1),col=textcolor)
	plot(histogram_FFPE_cut, ylim = c(0, 150000000), xaxt="n", yaxt="n",xlab="",ylab="", frame.plot=F,
	 col = "#1E90FF", cex = 2.4, type = "b", pch = 3, lty =1, lwd=60)
	axis(1,lty=1,lwd=50,at=(x_tick),labels=x_tick,las=0, cex.axis=2.5,col=textcolor,col.lab=textcolor,col.axis=textcolor)
	axis(2,lty=1,lwd=50,at=(y_tick),labels=y_tick,las=0, cex.axis=2.5,col=textcolor,col.lab=textcolor,col.axis=textcolor)
	par(new=T)
	lines(histogram_FF_cut, col = "#9ACD32", type = "b", pch = 3, lty =1, cex = 2.5, lwd=60)	
	mtext("Positions",side=2,line=4,cex=2.8,col=textcolor,font=2)
	mtext("Coverage",side=1,line=5,cex=2.8,col=textcolor,font=2)
	legend("right",c("FF","FFPE"), col=c("#9ACD32","#1E90FF"), cex=2.5, text.col=textcolor, pch = 3, lty =1, lwd=60)
	text( 5, 148000000 , sprintf("positions(cov=1) = %s", max(histogram_FF_cut) ), cex=2.2, col = "#9ACD32",adj=0)
	text( 5, 138000000 , sprintf("positions(cov=1) = %s", max(histogram_FFPE_cut) ), cex=2.2, col = "#1E90FF",adj=0)
dev.off()

say=sprintf("%s created!",pdffile);
print(say); 

