package circus;

import java.io.*;

public class DeFuse_1on1 {

	/**
	 * @param args
	 */
	public static void main(String[] args) {

		try {
			
			String samplenames = "/home/yhoang/workspace_Mariani/DeFuse_Perl/"+args[0] + ".txt";
			FileInputStream datei1 = new FileInputStream(samplenames);
			BufferedReader rdr1 = new BufferedReader(new InputStreamReader(
					datei1));
			String strLine;
			String[] samples = new String[1000];
			int n = 0;
			while ((strLine = rdr1.readLine()) != null) {
				samples[n] = strLine;
				// System.out.println(samples[n]);
				++n;
			}
			for (int m = 0; m < n; m++) {
				Circus circus = new Circus(Circus.GenomeHandler.HOMO_SAPIENS);
				circus.setSize(3000, 3000);		//do not change
				circus.setCenter(1500, 1500);	//do not change
				circus.setRadius(1200);
				circus.setTitle(samples[m], 15, 150, 200);
				
				String dateiname = "/home/yhoang/workspace_Mariani/DeFuse_Perl/compare_1on1/" + samples[m];
				FileInputStream datei = new FileInputStream(dateiname);
				BufferedReader rdr = new BufferedReader(new InputStreamReader(
						datei));

				int[][] chrom = new int[1000][5];
				int k = 0;
				while ((strLine = rdr.readLine()) != null) {
					String[] split = strLine.split("\t"); // dont save tabs
					for (int i = 0; i < split.length; i++) {
						if (split[i].equals("X")) {
							chrom[k][i] = 23;
						} else if (split[i].equals("MT")) {
							chrom[k][i] = 0;
						} else {
							chrom[k][i] = Integer.parseInt(split[i]);
						}
					}
					++k;
				}// 0:chr1, 1:chr2, 2:pos1, 3:pos2
				
				circus.init();

				float line_width = 10.0f;
				for (int l = 0; l < k; l++) {
					if (chrom[l][0] == chrom[l][1]) {
						// inside the cirlce, small arc
						// circus.initLines(0.0f, 4.0f, 1.0f, 10.0f);
						circus.initLines(0, 33, 249, 200); // blue
						circus.addLink(chrom[l][0], chrom[l][2], chrom[l][1],
								chrom[l][3], line_width, Circus.INNERLINE_TI);
					} else {
						// inside the circle, big arc
						circus.initLines(0.0f, 1.0f, 10f, 10.0f); // red
						circus.addLink(chrom[l][0], chrom[l][2], chrom[l][1],
								chrom[l][3], line_width, Circus.INNERLINE_TX);
					}
				}

				circus.saveImage("png_1on1/" + samples[m] + ".png");
				System.out.println("png_1on1/" + samples[m] + ".png created.");
			}
			System.out.println("Done.");
		} catch (Exception e) {
			e.printStackTrace();
		}

	}

}
/*
 * colors: circus.initLines(1.0f, 211.0f, 100.0f, 77.0f); orange
 * circus.initLines(211.0f, 100.0f, 77.0f, 0.8f); cream circus.initLines(0.0f,
 * 94.0f, 196.0f, 0.8f); brown circus.initLines(0.2f, 0.6f, 0.8f, 0.8f); light
 * green
 */

/*
 * // inside the colored arc, used to mark coverage of the // chromosome
 * circus.initLines(0.0f, 1.0f, 0.0f, 0.1f); circus.addLink(chrom[l][0],
 * chrom[l][3], chrom[l][1], chrom[l][3], Circus.OUTERLINE_INSIDE);
 * 
 * // outside the circle, used to mark fusions on the same // chromosome
 * circus.initLines(0.0f, 1.0f, 0.0f, 0.8f); circus.addLink(chrom[l][0],
 * chrom[l][2], chrom[l][1], chrom[l][3], Circus.OUTERLINE_OUTSIDE);
 */