package circus;

import java.io.*;

public class Highlight_denovo {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		String[] primarynames = new String[3];
		String[] metastasisnames = new String[6];
		int[][] primary = new int[200][4];
		int[][] metastasis = new int[200][4];
		int k = 0, n = 0, group = 4;
		boolean match = false;
		try {
			for (int q = 1; q <= group; q++) {
				// primary
				String primaryinput = "grouping/primary" + q + ".txt";
				FileInputStream datei1 = new FileInputStream(primaryinput);
				BufferedReader rdr1 = new BufferedReader(new InputStreamReader(
						datei1));
				String strLine;
				n = 0;
				k = 0;
				while ((strLine = rdr1.readLine()) != null) {
					primarynames[n] = strLine;
					++n;
				}
				for (int m = 0; m < n; m++) {
					String dateiname3 = "/home/yhoang/workspace_Mariani/DeFuse_Perl/compare_myself/"
							+ primarynames[m] + ".txt";
					FileInputStream datei3 = new FileInputStream(dateiname3);
					BufferedReader rdr3 = new BufferedReader(
							new InputStreamReader(datei3));
					while ((strLine = rdr3.readLine()) != null) {
						String[] split = strLine.split("\t"); // dont save tabs
						for (int i = 0; i <4; i++) {
								primary[k][i] = Integer.parseInt(split[i]);
						}
						++k;
					}// 0:chr1, 1:chr2, 2:pos1, 3:pos2
				}

				// metastasis
				String metastasisinput = "grouping/metastasis" + q + ".txt";
				FileInputStream datei2 = new FileInputStream(metastasisinput);
				BufferedReader rdr2 = new BufferedReader(new InputStreamReader(
						datei2));
				int v = 0;
				while ((strLine = rdr2.readLine()) != null) {
					metastasisnames[v] = strLine;
					++v;
				}

				for (int w = 0; w < v; w++) {
					Circus circus = new Circus(
							Circus.GenomeHandler.HOMO_SAPIENS);
					circus.setSize(3000, 3000); // do not change
					circus.setCenter(1500, 1500); // do not change
					circus.setRadius(1200);
					circus.setTitle(metastasisnames[w], 15, 150, 200);

					String dateiname4 = "/home/yhoang/workspace_Mariani/DeFuse_Perl/compare_myself/"
							+ metastasisnames[w] + ".txt";
					FileInputStream datei4 = new FileInputStream(dateiname4);
					BufferedReader rdr4 = new BufferedReader(
							new InputStreamReader(datei4));
					int a = 0;
					while ((strLine = rdr4.readLine()) != null) {
						String[] split = strLine.split("\t"); // dont save tabs
						for (int i = 0; i < 4; i++) {
								metastasis[a][i] = Integer.parseInt(split[i]);
						}
						++a;
					}// 0:chr1, 1:chr2, 2:pos1, 3:pos2

					circus.init();
					float line_width = 10.0f;
					for (int l = 0; l < a; l++) {
						for (int s = 0; s < k; s++) {
							if (metastasis[l][0] == primary[s][0]
									&& metastasis[l][1] == primary[s][1]
									&& metastasis[l][2] == primary[s][2]
									&& metastasis[l][3] == primary[s][3]) {
								match = true;
							}

							if (match == true) {
								if (metastasis[l][0] == metastasis[l][1]) {
									// inside the circle, small arc
									circus.initLines(0, 33, 249, 250); // blue
									circus.addLink(metastasis[l][0],
											metastasis[l][2], metastasis[l][1],
											metastasis[l][3], line_width,
											Circus.INNERLINE_TI);
								} else {// inside the circle, big arc
									circus.initLines(0.0f, 1.0f, 10f, 10.0f); // red
									circus.addLink(metastasis[l][0],
											metastasis[l][2], metastasis[l][1],
											metastasis[l][3], line_width,
											Circus.INNERLINE_TX);
								}
							} else {
								if (metastasis[l][0] == metastasis[l][1]) {
									// inside the circle, small arc
									circus.initLines(255, 255, 0, 200); // yellow
									circus.addLink(metastasis[l][0],
											metastasis[l][2], metastasis[l][1],
											metastasis[l][3], line_width,
											Circus.INNERLINE_TI);
								} else {// inside the circle, big arc
									circus.initLines(255, 255, 0, 200); // yellow
									circus.addLink(metastasis[l][0],
											metastasis[l][2], metastasis[l][1],
											metastasis[l][3], line_width,
											Circus.INNERLINE_TX);
								}
							}
						}
						match = false;
					}
					circus.saveImage("png_group" + q + "/" + metastasisnames[w]
							+ ".png");
					System.out.println("png_group" + q + "/"
							+ metastasisnames[w] + ".png created.");
				}
			}
			System.out.println("Done.");
		} catch (Exception e) {
			e.printStackTrace();
		}

	}

}
