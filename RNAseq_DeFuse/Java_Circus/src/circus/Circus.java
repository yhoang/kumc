package circus;


import java.awt.Color;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics2D;
import java.awt.BasicStroke;
import java.awt.RenderingHints;
import java.awt.geom.Arc2D;
import java.awt.geom.Ellipse2D;
import java.awt.geom.Line2D;
import java.awt.geom.Path2D;
import java.awt.geom.Rectangle2D;
import java.awt.image.BufferedImage;
import java.io.File;

import javax.imageio.ImageIO;



public class Circus {

	private GenomeHandler gnmHandler;

	
	private BufferedImage bImage;
	private Graphics2D graphics;
	
	private Color color_background = Color.BLACK;
	private Color color_lines = new Color(100,100,100,100);
	
	private int width;
	private int height;
	private int centerX;
	private int centerY;
	private int weightX;
	private int weightY;
	private double radius;

	private double bottomRadius_CS;
	private double bottomRadius_CL;
	private double bottomRadius_TI;
	private double bottomRadius_TX;

	private double weightRadius_TI;
//	private double weightRadius_TX;
	private double weightRadius_CL;
	private double weightRadius_CS;
	
	private double thickness = 0.05;		//thickness of circle
	private double basesPerDegree = 0.0;
	
	private long ideoPositions[];
	
	private String title = "";
	private int title_x = 0;
	private int title_y = 0;
	private int title_s = 0;

	public static final int INNERLINE_TX = 1; 
	public static final int INNERLINE_TI = 2; 
	public static final int OUTERLINE_OUTSIDE = 3; 
	public static final int OUTERLINE_INSIDE = 4;
	
	public static final int STRAIGHTLINE = 1; 
	public static final int QUADLINE = 2; 
	
	private Font titleFont;
	
	private double degreeTurn = 90;
	private int degreeDirection = -1;

	

	
	public Circus(int genomeID) {

		gnmHandler = new GenomeHandler(genomeID);
		
	}
	
	
	public void setTitle(String title, int size, int x, int y) {
		this.title = title;
		this.title_s = size;
		this.title_x = x;
		this.title_y = y;
	}

	
	public void setSize(int width, int height) {
		this.width = width;
		this.height = height;
	}
	
	
	public void setCenter(int x, int y) {
		this.centerX = x;
		this.centerY = y;
	}
	
	
	
	public void setRadius(int rad) {
		this.radius = rad;
	}

	
	
	public void init() {
		weightX = centerX;
		weightY = centerY;

		titleFont = new Font("SansSerif", Font.PLAIN, 55);
		bottomRadius_TX = radius - radius*thickness - 10;

		bottomRadius_TI = radius - radius*thickness - 10;
		weightRadius_TI = centerX * 0.50;
		
		bottomRadius_CL = radius + 10;
		weightRadius_CL = centerX * 0.75 + centerX * 0.25;
		
		bottomRadius_CS = radius - radius*thickness + 5;
		weightRadius_CS = radius - 5;
		
		bImage = new BufferedImage(width, height, BufferedImage.TYPE_INT_ARGB);
		graphics = bImage.createGraphics();

		graphics.setBackground( color_background );
		graphics.clearRect(0, 0, width, height);
		
		graphics.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
		graphics.setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_ON);
		
		long[] chr = new long[ gnmHandler.getChromosomeAmount() ];
		
		for (int c = 0; c < gnmHandler.getChromosomeAmount(); c++ ) {
			chr[c] = gnmHandler.getChromosomeLength(c);
		}
		
		initIdeogram(chr, 10000000);
		initLines(0.0f, 1.0f, 0.0f, 0.1f);

	}
	
	
	
	public void initIdeogram(long[] arcs, int space) {
	
		long totalLength = 0;
		ideoPositions = new long[arcs.length];

		float hueInc = 1.0f/arcs.length;
		//System.out.println(hueInc);
		
		ideoPositions[0] = space/2;

		for (int a = 0; a < arcs.length; a++) {		// changed ifrom int a = 1, now with chr#=0
			ideoPositions[a] = degreeDirection * (space + totalLength);
			totalLength += degreeDirection * (space + arcs[a]);
			//System.out.println(ideoPositions[a]);
		}
		
		
		basesPerDegree = (totalLength) / (360.0);
		//System.out.println(arcPerSegment);
		
		//arcs.length
		for (int a = 0; a < arcs.length; a++) {
			graphics.setColor( Color.getHSBColor(a*hueInc, 0.8f, 1.0f) );
			Arc2D.Double arc2D = new Arc2D.Double(centerX - radius, centerY - radius, radius*2, radius*2, degreeTurn + (ideoPositions[a]/basesPerDegree), arcs[a]/basesPerDegree, Arc2D.PIE); 
			graphics.fill( arc2D );
				
			graphics.setColor( Color.WHITE );
			
			drawTitle();
			
			double radS = ((ideoPositions[a])/basesPerDegree) + degreeTurn;
			double radE = ((ideoPositions[a]+arcs[a])/basesPerDegree) + degreeTurn;
			
			double xS = centerX + ( radius * Math.cos(Math.toRadians(radS)) );
			double yS = centerY - ( radius * Math.sin(Math.toRadians(radS)) );
			double xE = centerX + ( radius * Math.cos(Math.toRadians(radE)) );
			double yE = centerY - ( radius * Math.sin(Math.toRadians(radE)) );
			double xC = centerX + ( (radius*1.05) * Math.cos(Math.toRadians((radS+radE)/2)) );
			double yC = centerY - ( (radius*1.05) * Math.sin(Math.toRadians((radS+radE)/2)) );
			
			graphics.setColor( Color.WHITE );
			Line2D.Double line = new Line2D.Double(xS, yS, centerX, centerY);
			graphics.draw(line);
			line = new Line2D.Double(xE, yE, centerX, centerY);
			graphics.draw(line);
			
			graphics.setFont(titleFont.deriveFont(38f));

			FontMetrics fm = graphics.getFontMetrics();
			Rectangle2D textsize = fm.getStringBounds(gnmHandler.getChromosomeName(a), graphics);
			double xO = textsize.getWidth() / 2.0;
			double yO = textsize.getHeight() / 2.0;
			
			graphics.drawString(gnmHandler.getChromosomeName(a), (int) (xC - xO), (int) (yC + yO));
		}
		
		graphics.setColor( color_background );
		Ellipse2D.Double oval = new Ellipse2D.Double(centerX - radius + radius*thickness, centerY - radius + radius*thickness, (radius - radius*thickness)*2, (radius - radius*thickness)*2);
		graphics.fill(oval);
		
		graphics.setColor( Color.WHITE );
		for (int a = 1; a < arcs.length; a++) {
			graphics.draw( new Arc2D.Double(centerX - radius, centerY - radius, radius*2, radius*2, degreeTurn + (ideoPositions[a]/basesPerDegree), arcs[a]/basesPerDegree, Arc2D.OPEN) );
			graphics.draw( new Arc2D.Double(centerX - radius + radius*thickness, centerY - radius + radius*thickness, radius*2 - radius*thickness*2, radius*2 - radius*thickness*2, degreeTurn + (ideoPositions[a]/basesPerDegree), arcs[a]/basesPerDegree, Arc2D.OPEN) );
		}
		
	}
	
	
	
	
	public void addLink(int chr1, int pos1, int chr2, int pos2, float line_width, int type) {
		
		double angledPos1 = ((ideoPositions[chr1]+pos1)/basesPerDegree) + degreeTurn;
		double angledPos2 = ((ideoPositions[chr2]+pos2)/basesPerDegree) + degreeTurn;
		
		double x1 = 0.0, x2 = 0.0, y1 = 0.0, y2 = 0.0, wX = 0.0, wY = 0.0;
		
		if ( type == INNERLINE_TX ) {
			x1 = centerX + ( bottomRadius_TX * Math.cos(Math.toRadians(angledPos1)) );
			y1 = centerY - ( bottomRadius_TX * Math.sin(Math.toRadians(angledPos1)) );
			
			x2 = centerX + ( bottomRadius_TX * Math.cos(Math.toRadians(angledPos2)) );
			y2 = centerY - ( bottomRadius_TX * Math.sin(Math.toRadians(angledPos2)) );

			wX = weightX;
			wY = weightY;

			Path2D.Double pathB = new Path2D.Double(Path2D.WIND_EVEN_ODD);
			pathB.moveTo( x1, y1 );
			pathB.quadTo( wX, wY, x2, y2);
			
			graphics.setStroke(new BasicStroke(line_width));
			graphics.draw(pathB);
		} else if ( type == INNERLINE_TI ) {
			x1 = centerX + ( (bottomRadius_TI) * Math.cos(Math.toRadians(angledPos1)) );
			y1 = centerY - ( (bottomRadius_TI) * Math.sin(Math.toRadians(angledPos1)) );

			x2 = centerX + ( (bottomRadius_TI) * Math.cos(Math.toRadians(angledPos2)) );
			y2 = centerY - ( (bottomRadius_TI) * Math.sin(Math.toRadians(angledPos2)) );
			
			double midAngle = (angledPos1 + angledPos2) / 2;
			
			wX = centerX + ( (weightRadius_TI) * Math.cos(Math.toRadians(midAngle)) );
			wY = centerY - ( (weightRadius_TI) * Math.sin(Math.toRadians(midAngle)) );

			Path2D.Double pathB = new Path2D.Double(Path2D.WIND_EVEN_ODD);
			pathB.moveTo( x1, y1 );
			pathB.quadTo( wX, wY, x2, y2);
			
			graphics.setStroke(new BasicStroke(line_width));
			graphics.draw(pathB);
		} else if ( type == OUTERLINE_OUTSIDE ) {
			x1 = centerX + ( (bottomRadius_CL) * Math.cos(Math.toRadians(angledPos1)) );
			y1 = centerY - ( (bottomRadius_CL) * Math.sin(Math.toRadians(angledPos1)) );

			x2 = centerX + ( (bottomRadius_CL) * Math.cos(Math.toRadians(angledPos2)) );
			y2 = centerY - ( (bottomRadius_CL) * Math.sin(Math.toRadians(angledPos2)) );
			
			double midAngle = (angledPos1 + angledPos2) / 2;
			
			wX = centerX + ( (weightRadius_CL) * Math.cos(Math.toRadians(midAngle)) );
			wY = centerY - ( (weightRadius_CL) * Math.sin(Math.toRadians(midAngle)) );

			Path2D.Double pathB = new Path2D.Double(Path2D.WIND_EVEN_ODD);
			pathB.moveTo( x1, y1 );
			pathB.quadTo( wX, wY, x2, y2);
			
			graphics.setStroke(new BasicStroke(line_width));
			graphics.draw(pathB);
		} else if ( type == OUTERLINE_INSIDE ) {
			x1 = centerX + ( (bottomRadius_CS) * Math.cos(Math.toRadians(angledPos1)) );
			y1 = centerY - ( (bottomRadius_CS) * Math.sin(Math.toRadians(angledPos1)) );
	
			x2 = centerX + ( (bottomRadius_CS) * Math.cos(Math.toRadians(angledPos2)) );
			y2 = centerY - ( (bottomRadius_CS) * Math.sin(Math.toRadians(angledPos2)) );
			
			double midAngle = (angledPos1 + angledPos2) / 2;
			
			wX = centerX + ( (weightRadius_CS) * Math.cos(Math.toRadians(midAngle)) );
			wY = centerY - ( (weightRadius_CS) * Math.sin(Math.toRadians(midAngle)) );

			Path2D.Double pathB = new Path2D.Double(Path2D.WIND_EVEN_ODD);
			pathB.moveTo( x1, y1 );
			pathB.lineTo( wX, wY );
			
			graphics.setStroke(new BasicStroke(line_width));
			graphics.draw(pathB);
		}
		
	}
	
	
	
	public void drawLine(double x1, double y1, double x2, double y2) {
		
		Line2D.Double line = new Line2D.Double(x1, y1, x2, y2);
		graphics.draw(line);
		
	}
	
	

	public void saveImage(String fileName) throws Exception {
		File f = new File(fileName);
		f.mkdirs();
		ImageIO.write(bImage, "png", f);
		
	}
	


	public void initLines(int r, int g, int b, int a) {
		color_lines = new Color(r, g, b, a);
		graphics.setColor( color_lines );
	}

	
	
	public void initLines(float h, float s, float b, float a) {
		Color col = Color.getHSBColor(h, s, b);
		color_lines = new Color(col.getRed(), col.getGreen(), col.getBlue(), Math.min((int) (a*255), 255));
		graphics.setColor( color_lines );
	}



	public void drawString(String str, int x, int y, float size) {
		graphics.setFont( titleFont.deriveFont(size) );
		graphics.drawString(str, x, y);

	}



	public BufferedImage getBufferedImage() {
		return bImage;
	}



	public void drawTitle() {
		graphics.setFont( titleFont.deriveFont(title_s) );
		graphics.drawString(title, title_x, title_y);
	}



	public class GenomeHandler {

		/* Rattus Norvegicus RGSC 3.4 */
		private static final int RATTUS_NORVEGICUS_CHROMLENGTH_1  = 267910886;		//chromosome:RGSC3.4:1:1:267,910,886:1
		private static final int RATTUS_NORVEGICUS_CHROMLENGTH_2  = 258207540;		//chromosome:RGSC3.4:2:1:258,207,540:1
		private static final int RATTUS_NORVEGICUS_CHROMLENGTH_3  = 171063335;		//chromosome:RGSC3.4:3:1:171,063,335:1
		private static final int RATTUS_NORVEGICUS_CHROMLENGTH_4  = 187126005;		//chromosome:RGSC3.4:4:1:187,126,005:1
		private static final int RATTUS_NORVEGICUS_CHROMLENGTH_5  = 173096209;		//chromosome:RGSC3.4:5:1:173,096,209:1
		private static final int RATTUS_NORVEGICUS_CHROMLENGTH_6  = 147636619;		//chromosome:RGSC3.4:6:1:147,636,619:1
		private static final int RATTUS_NORVEGICUS_CHROMLENGTH_7  = 143002779;		//chromosome:RGSC3.4:7:1:143,002,779:1
		private static final int RATTUS_NORVEGICUS_CHROMLENGTH_8  = 129041809;		//chromosome:RGSC3.4:8:1:129,041,809:1
		private static final int RATTUS_NORVEGICUS_CHROMLENGTH_9  = 113440463;		//chromosome:RGSC3.4:9:1:113,440,463:1
		private static final int RATTUS_NORVEGICUS_CHROMLENGTH_10 = 110718848;		//chromosome:RGSC3.4:10:1:110,718,848:1
		private static final int RATTUS_NORVEGICUS_CHROMLENGTH_11 =  87759784;		//chromosome:RGSC3.4:11:1:87,759,784:1
		private static final int RATTUS_NORVEGICUS_CHROMLENGTH_12 =  46782294;		//chromosome:RGSC3.4:12:1:46,782,294:1
		private static final int RATTUS_NORVEGICUS_CHROMLENGTH_13 = 111154910;		//chromosome:RGSC3.4:13:1:111,154,910:1
		private static final int RATTUS_NORVEGICUS_CHROMLENGTH_14 = 112194335;		//chromosome:RGSC3.4:14:1:112,194,335:1
		private static final int RATTUS_NORVEGICUS_CHROMLENGTH_15 = 109758846;		//chromosome:RGSC3.4:15:1:109,758,846:1
		private static final int RATTUS_NORVEGICUS_CHROMLENGTH_16 =  90238779;		//chromosome:RGSC3.4:16:1:90,238,779:1
		private static final int RATTUS_NORVEGICUS_CHROMLENGTH_17 =  97296363;		//chromosome:RGSC3.4:17:1:97,296,363:1
		private static final int RATTUS_NORVEGICUS_CHROMLENGTH_18 =  87265094;		//chromosome:RGSC3.4:18:1:87,265,094:1
		private static final int RATTUS_NORVEGICUS_CHROMLENGTH_19 =  59218465;		//chromosome:RGSC3.4:19:1:59,218,465:1
		private static final int RATTUS_NORVEGICUS_CHROMLENGTH_20 =  55268282;		//chromosome:RGSC3.4:20:1:55,268,282:1
		private static final int RATTUS_NORVEGICUS_CHROMLENGTH_MT =     16300;		//chromosome:RGSC3.4:MT:1:16,300:1
		private static final int RATTUS_NORVEGICUS_CHROMLENGTH_X  = 160699376;		//chromosome:RGSC3.4:X:1:160,699,376:1
		private static final int RATTUS_NORVEGICUS_CHROMAMOUNT    =        22;	

		/* Homo Sapiens */
		private static final int HOMO_SAPIENS_CHROMLENGTH_1  = 247249719;
		private static final int HOMO_SAPIENS_CHROMLENGTH_2  = 242951149;
		private static final int HOMO_SAPIENS_CHROMLENGTH_3  = 199501827;
		private static final int HOMO_SAPIENS_CHROMLENGTH_4  = 191273063;
		private static final int HOMO_SAPIENS_CHROMLENGTH_5  = 180857866;
		private static final int HOMO_SAPIENS_CHROMLENGTH_6  = 170899992;
		private static final int HOMO_SAPIENS_CHROMLENGTH_7  = 158821424;
		private static final int HOMO_SAPIENS_CHROMLENGTH_8  = 146274826;
		private static final int HOMO_SAPIENS_CHROMLENGTH_9  = 140273252;
		private static final int HOMO_SAPIENS_CHROMLENGTH_10 = 135374737;
		private static final int HOMO_SAPIENS_CHROMLENGTH_11 = 134452384;
		private static final int HOMO_SAPIENS_CHROMLENGTH_12 = 132349534;
		private static final int HOMO_SAPIENS_CHROMLENGTH_13 = 114142980;
		private static final int HOMO_SAPIENS_CHROMLENGTH_14 = 106368585;
		private static final int HOMO_SAPIENS_CHROMLENGTH_15 = 100338915;
		private static final int HOMO_SAPIENS_CHROMLENGTH_16 =  88827254;
		private static final int HOMO_SAPIENS_CHROMLENGTH_17 =  78774742;
		private static final int HOMO_SAPIENS_CHROMLENGTH_18 =  76117153;
		private static final int HOMO_SAPIENS_CHROMLENGTH_19 =  63811651;
		private static final int HOMO_SAPIENS_CHROMLENGTH_20 =  62435964;
		private static final int HOMO_SAPIENS_CHROMLENGTH_21 =  46944323;
		private static final int HOMO_SAPIENS_CHROMLENGTH_22 =  49691432;
		private static final int HOMO_SAPIENS_CHROMLENGTH_MT =     16571;
		private static final int HOMO_SAPIENS_CHROMLENGTH_X  = 154913754;
	//	private static final int HOMO_SAPIENS_CHROMLENGTH_Y  =  57443437;
		private static final int HOMO_SAPIENS_CHROMAMOUNT    =        24;	

		
		private int genome = -1;
		private long genomeLength = 0;
		
		public static final int RATTUS_NORVEGICUS = 1001;
		public static final int HOMO_SAPIENS = 2001;
		
		
		
		public GenomeHandler(int genomeID) {
			genome = genomeID;

			genomeLength = 0;
			
			for (int c = 0; c < getChromosomeAmount(); c++) genomeLength += getChromosomeLength(c);
		}
		
		
		
		public int getChromosomeAmount() {
			
			switch (genome) {
			
				case RATTUS_NORVEGICUS: return RATTUS_NORVEGICUS_CHROMAMOUNT;
				case HOMO_SAPIENS: return HOMO_SAPIENS_CHROMAMOUNT;
			
			}
			
			return -1;
		}
		
		
		
		public int getChromosomeLength(int chr) {
			
			switch (genome) {
			
			case RATTUS_NORVEGICUS: {
				switch (chr) {
					case 0: return RATTUS_NORVEGICUS_CHROMLENGTH_MT;
					case 1: return RATTUS_NORVEGICUS_CHROMLENGTH_1;
					case 2: return RATTUS_NORVEGICUS_CHROMLENGTH_2;
					case 3: return RATTUS_NORVEGICUS_CHROMLENGTH_3;
					case 4: return RATTUS_NORVEGICUS_CHROMLENGTH_4;
					case 5: return RATTUS_NORVEGICUS_CHROMLENGTH_5;
					case 6: return RATTUS_NORVEGICUS_CHROMLENGTH_6;
					case 7: return RATTUS_NORVEGICUS_CHROMLENGTH_7;
					case 8: return RATTUS_NORVEGICUS_CHROMLENGTH_8;
					case 9: return RATTUS_NORVEGICUS_CHROMLENGTH_9;
					case 10: return RATTUS_NORVEGICUS_CHROMLENGTH_10;
					case 11: return RATTUS_NORVEGICUS_CHROMLENGTH_11;
					case 12: return RATTUS_NORVEGICUS_CHROMLENGTH_12;
					case 13: return RATTUS_NORVEGICUS_CHROMLENGTH_13;
					case 14: return RATTUS_NORVEGICUS_CHROMLENGTH_14;
					case 15: return RATTUS_NORVEGICUS_CHROMLENGTH_15;
					case 16: return RATTUS_NORVEGICUS_CHROMLENGTH_16;
					case 17: return RATTUS_NORVEGICUS_CHROMLENGTH_17;
					case 18: return RATTUS_NORVEGICUS_CHROMLENGTH_18;
					case 19: return RATTUS_NORVEGICUS_CHROMLENGTH_19;
					case 20: return RATTUS_NORVEGICUS_CHROMLENGTH_20;
					case 21: return RATTUS_NORVEGICUS_CHROMLENGTH_X;
					case 91: return RATTUS_NORVEGICUS_CHROMLENGTH_X;
				}
			}
			
			case HOMO_SAPIENS: {
				switch (chr) {
					case 0: return HOMO_SAPIENS_CHROMLENGTH_MT;
					case 1: return HOMO_SAPIENS_CHROMLENGTH_1;
					case 2: return HOMO_SAPIENS_CHROMLENGTH_2;
					case 3: return HOMO_SAPIENS_CHROMLENGTH_3;
					case 4: return HOMO_SAPIENS_CHROMLENGTH_4;
					case 5: return HOMO_SAPIENS_CHROMLENGTH_5;
					case 6: return HOMO_SAPIENS_CHROMLENGTH_6;
					case 7: return HOMO_SAPIENS_CHROMLENGTH_7;
					case 8: return HOMO_SAPIENS_CHROMLENGTH_8;
					case 9: return HOMO_SAPIENS_CHROMLENGTH_9;
					case 10: return HOMO_SAPIENS_CHROMLENGTH_10;
					case 11: return HOMO_SAPIENS_CHROMLENGTH_11;
					case 12: return HOMO_SAPIENS_CHROMLENGTH_12;
					case 13: return HOMO_SAPIENS_CHROMLENGTH_13;
					case 14: return HOMO_SAPIENS_CHROMLENGTH_14;
					case 15: return HOMO_SAPIENS_CHROMLENGTH_15;
					case 16: return HOMO_SAPIENS_CHROMLENGTH_16;
					case 17: return HOMO_SAPIENS_CHROMLENGTH_17;
					case 18: return HOMO_SAPIENS_CHROMLENGTH_18;
					case 19: return HOMO_SAPIENS_CHROMLENGTH_19;
					case 20: return HOMO_SAPIENS_CHROMLENGTH_20;
					case 21: return HOMO_SAPIENS_CHROMLENGTH_21;
					case 22: return HOMO_SAPIENS_CHROMLENGTH_22;
					case 23: return HOMO_SAPIENS_CHROMLENGTH_X;
		//			case 24: return HOMO_SAPIENS_CHROMLENGTH_Y;
		//			case 91: return HOMO_SAPIENS_CHROMLENGTH_X;
			//		case 92: return HOMO_SAPIENS_CHROMLENGTH_Y;
				}
			}
			
			}
			
			return -1;
		}
		
		
		
		public String getChromosomeName(int chr) {
			
			switch (genome) {
			
			case RATTUS_NORVEGICUS: {
				switch (chr) {
					case 0: return "MT";
					case 1: return "1";
					case 21: return "X";
					case 91: return "X";
					default: return ""+chr;
				}
			}
			
			case HOMO_SAPIENS: {
				switch (chr) {
					case 0: return "MT";
					case 23: return "X";
//					case 24: return "Y";
//					case 91: return "X";
//					case 92: return "Y";
					default: return ""+chr;
				}
			}
			
			}
			
			return null;
		}
		
		
		
		public int getSpecialChromosomeNo(char type) {
			
			switch (genome) {

				case RATTUS_NORVEGICUS: {
					switch(type) {
						case 'X': return 91;
						case 'M': return 0;
					}
					break;
				}
				
				case HOMO_SAPIENS: {
					switch(type) {
						case 'X': return 91;
//						case 'Y': return 92;
//						case 'M': return 0;
					}
					break;
				}
			}		
			
			return -1;
		}
		

		
		public int convertSpecialChromosomeNo(int type) {
			
			switch (genome) {

				case RATTUS_NORVEGICUS: {
					switch(type) {
						case 90: return 0;
						case 91: return 21;
					}
					break;
				}
				
				case HOMO_SAPIENS: {
					switch(type) {
						case 90: return 0;
						case 91: return 23;
//						case 92: return 24;
					}
					break;
				}
			}		
			
			return type;
		}
		
		
		
		public long getGenomeLength() {
			
			switch (genome) {
			
				case RATTUS_NORVEGICUS: return genomeLength;
				case HOMO_SAPIENS: return genomeLength;
			
			}
			
			return -1;
		}
		
	}
	
}