

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

public class GISLog2Norm {

	public static void main(String[] args) throws IOException {
		String GIS_CHANEL = "134N_Ref";
		File quantIntFile = new File("tutorial/Acetyl_QuantMethod1_95.0.txt");
		BufferedReader BR = new BufferedReader(new FileReader(quantIntFile));
		BufferedWriter BWGIS = new BufferedWriter(new FileWriter(quantIntFile.getAbsolutePath().replace(".txt", ".rmNullRef.GISLog2Transformed.tsv")));
		BufferedWriter BW = new BufferedWriter(new FileWriter(quantIntFile.getAbsolutePath().replace(".txt", ".rmNullRef.tsv")));
		String line = null;
		
		String setinfo = BR.readLine();
		String[] setinfos = setinfo.split("\t");
		String chanel = BR.readLine();
		String[] chanels = chanel.split("\t");
		
		ArrayList<Integer> gisIndices = new ArrayList<Integer>();
		
		int totalSamples = 0;
		for(int i=0; i<chanels.length; i++) {
			
			// sample count
			if(!chanels[i].contains("Ref") && i!=0) {
				totalSamples++;
			}
			
			if(chanels[i].equalsIgnoreCase(GIS_CHANEL)) {
				gisIndices.add(i);
			} else {
				if(i!=0) {
					BWGIS.append("\t");
				}
				BWGIS.append(setinfos[i]);
			}
		}
		System.out.println(totalSamples);
		BWGIS.newLine();
		for(int i=0; i<chanels.length; i++) {
			if(!chanels[i].equalsIgnoreCase(GIS_CHANEL)) {
				if(i!=0) {
					BWGIS.append("\t");
				}
				BWGIS.append(chanels[i]);
			}
		}
		BWGIS.newLine();
		
		System.out.println("Size of GIS: "+gisIndices.size());
		
		BW.append(setinfo);
		BW.newLine();
		BW.append(chanel);
		BW.newLine();
		
		int[] setSites = new int[27];
		while((line = BR.readLine()) != null) {
			String[] fields = line.split("\t");
			String siteId = fields[0];
			int startIdx = -1;
			int endIdx = -1;
			
			ArrayList<String> quants = new ArrayList<String>();
			ArrayList<String> gisNorm = new ArrayList<String>();
			boolean isDropped = true;
			int setIdx = 0;
			for(int gisIdx = 0; gisIdx<gisIndices.size(); gisIdx++) {
				endIdx = gisIndices.get(gisIdx);
				startIdx = endIdx - 15;
				
				double sum = 0;
				double gis = 0;
				ArrayList<Double> quantValues = new ArrayList<Double>();
				for(int i=startIdx; i<=endIdx; i++) {
					double quantVal = Double.parseDouble(fields[i]);
					
					quantValues.add(quantVal);
					sum += quantVal;
					
					if(i==endIdx) {
						gis = quantVal;
					}
				}
				
				boolean isDroppedEntry = false;
				// drop not quantifiable
				if(sum == 0 || gis == 0) {
					isDroppedEntry = true;
				} else {
					isDropped = false;
					setSites[setIdx]++;
				}
				
				for(int i=startIdx; i<=endIdx; i++) {
					if(isDroppedEntry) {
						quants.add("NA");
					} else {
						double quantVal = Double.parseDouble(fields[i]);
						if(quantVal == 0) {
							quants.add("NA");
						} else {
							quants.add(quantVal+"");
						}
					}
				}
				
				for(int i=startIdx; i<endIdx; i++) {
					if(isDroppedEntry) {
						gisNorm.add("NA");
					} else {
						double quantVal = Double.parseDouble(fields[i]);
						if(quantVal == 0) {
							gisNorm.add("NA");
						} else {
							gisNorm.add((Math.log(quantVal/gis)/Math.log(2))+"");
						}
					}
				}
				
				setIdx++;
			}
			
			if(!isDropped) {
				BW.append(siteId);
				BWGIS.append(siteId);
				
				for(int i=0; i<quants.size(); i++) {
					BW.append("\t"+quants.get(i));
				}
				BW.newLine();
				for(int i=0; i<gisNorm.size(); i++) {
					BWGIS.append("\t"+gisNorm.get(i));
				}
				BWGIS.newLine();
			}
		}
		
		System.out.println();
		System.out.println("Batch\tNumOfSites");
		for(int i=0; i<setSites.length; i++) {
			int setNum = (i+1);
			String set = "Set";
			if(setNum < 10) {
				set += "0";
			}
			set += setNum;
			
			System.out.println(set+"\t"+setSites[i]);
		}
		
		BR.close();
		BW.close();
		BWGIS.close();
		
	}
	
}
