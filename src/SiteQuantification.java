

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Hashtable;

class STAT {
	ArrayList<String> sets = new ArrayList<String>();
	ArrayList<Integer> sites = new ArrayList<Integer>();
	ArrayList<Integer> psms = new ArrayList<Integer>();
}

public class SiteQuantification {

	public static String[] PSM_HEADER = null;
	
	public static boolean isUseShared = false;
	public static double ptmRSThreshold = 95;
	public static STAT stat = new STAT();
	
	
	public static String base = null;
	public static String base2 = null;
	
	public static int numOfBatches = 0;
	
	/**
	 * Phospho 
	 **/
	public static void loadPhosphoParam (int numOfBatches) {
		// input base file name
		base = "tutorial/BRCA_Phospho_";
		// output base file name
		base2 = "tutorial/Phospho_QuantMethod1_";
		
		SiteQuantification.numOfBatches = numOfBatches;
	}
	
	public static void loadAcetylParam (int numOfBatches) {
		// input base file name
		base = "tutorial/BRCA_Acetyl_";
		// output base file name
		base2 = "tutorial/Acetyl_QuantMethod1_";
		
		SiteQuantification.numOfBatches = numOfBatches;
	}
	
	
	public static void main(String[] args) throws IOException {
		
		// load param
		loadAcetylParam(1);
//		loadPhosphoParam(27);
		
		File fastaFile = new File("tutorial/RefSeq_GRCh38_2023_04_protein_with_381_contaminants.fasta");
		File channelFile = new File("tutorial/batch_table.tsv");
		Hashtable<String, String> proteinMap = loadProteinFasta(fastaFile);
		
		// set_channel => sampleName
		Hashtable<String, String> channelMap = sampleMapper(channelFile);
		
		System.out.println("ptmRS threhold: "+ptmRSThreshold);
		ArrayList<String> outputNames = new ArrayList<String>();
		for(int i=1; i<=numOfBatches; i++) {
			
			String set = i+"";
			if(set.length() < 2) {
				set = "0"+i;
			}
			set = "Set"+set;
			
			File quanFile = new File(base+set+"_QuanSpectra.txt");
			File psmFile = new File(base+set+"_PSMs.txt");
			File output = new File(psmFile.getAbsolutePath().replace("_PSMs.txt", "_QuantMethod1_ptm"+ptmRSThreshold+".txt"));
			
			
			Hashtable<String, String> quanMap = loadQuanSpectra(quanFile);
			
			Hashtable<String, ArrayList<String[]>> sites = loadPSM(psmFile, quanMap, proteinMap);
			
			// update stat
			stat.sets.add(set);
			writeMethod1(sites, output);
			outputNames.add(output.getAbsolutePath());
		}
		
		// merge as single file
		mergeSets(outputNames, channelMap, base2+ptmRSThreshold+".txt");
		
		System.out.println();
		System.out.println("Batch\tNumOfSites\tNumOfPSMs");
		for(int i=0; i<stat.sets.size(); i++) {
			System.out.println(stat.sets.get(i)+"\t"+stat.sites.get(i)+"\t"+stat.psms.get(i));
		}
	}
	
	// Set01_Chanel => SampleName
	public static Hashtable<String, String> sampleMapper (File file) throws IOException {
		Hashtable<String, String> map = new Hashtable<String, String>();
		
		BufferedReader BR = new BufferedReader(new FileReader(file));
		String line = null;
		
		String[] header = BR.readLine().split("\t");
		
		while((line = BR.readLine()) != null) {
			String[] fields = line.split("\t");
			
			String channel = fields[0];
			for(int i=1; i<fields.length; i++) {
				String set = header[i];
				String sampleName = fields[i];
				
				String key = set+"_"+channel;
				map.put(key, sampleName);
			}
		}
		
		
		BR.close();
		
		return map;
	}
	
	public static void mergeSets (ArrayList<String> listFiles, Hashtable<String, String> channelMap, String outputName) throws IOException {
		BufferedWriter BW = new BufferedWriter(new FileWriter(outputName));
		Hashtable<String, String[]> recordMap = new Hashtable<String, String[]>();
		
		BW.append("SetInfo");
		for(String set : stat.sets) {
			//16plex
			for(int i=0; i<16; i++) {
				BW.append("\t"+set);
			}
		}
		BW.newLine();
		
		
		BW.append("SiteId");
		for(int i=0; i<listFiles.size(); i++) {
			String set = stat.sets.get(i);
			String fileName = listFiles.get(i);
			
			BufferedReader BR = new BufferedReader(new FileReader(fileName));
			String line = null;
			String[] header = BR.readLine().split("\t");
			
			// idx 0 is "SiteId"
			// from 1 to 17 : channel id (but Abundance: 134N...)
			for(int idx=1; idx<header.length; idx++) {
				String key = set+"_"+header[idx].split("\\s")[1].trim();
				String sampleName = channelMap.get(key);
				
				if(sampleName == null) {
					System.out.println(key +" something is wrong! to map channel map");
				}
				BW.append("\t").append(sampleName);
			}
			
			while((line = BR.readLine()) != null) {
				String[] fields = line.split("\t");
				String siteId = fields[0];
				String[] quans = recordMap.get(siteId);
				if(quans == null) {
					quans = new String[listFiles.size()];
					recordMap.put(siteId, quans);
				}
				quans[i] = line;
			}
			
			BR.close();
		}
		BW.newLine();
		
		recordMap.forEach((siteId, quans) -> {
			try {
				BW.append(siteId);
				for(int i=0; i<quans.length; i++) {
					if(quans[i] == null) {
						for(int j=0; j<16; j++) {
							BW.append("\t0");
						}
					} else {
						String[] quanValues = quans[i].split("\t");
						for(int j=1; j<quanValues.length; j++) {
							BW.append("\t"+quanValues[j]);
						}
					}
				}
				BW.newLine();
			}catch(IOException ioe) {
				
			}
		});
		
		BW.close();
	}
	
	public static void writeMethod1 (Hashtable<String, ArrayList<String[]>> sites, File output) throws IOException  {
		BufferedWriter BW = new BufferedWriter(new FileWriter(output));

		int[] totalSites = {0};
		
		int[] quantStartIdx = new int[1];
		int[] quantSize = {16};
		for(int i=0; i<PSM_HEADER.length; i++) {
			if(PSM_HEADER[i].equalsIgnoreCase("Abundance: 126")) {
				quantStartIdx[0] = i;
			} 
		}
		
		BW.append("SiteId");
		for(int i=0; i<quantSize[0]; i++) {
			BW.append("\t").append(PSM_HEADER[i+quantStartIdx[0]]);
		}
		BW.newLine();
		
		
		sites.forEach((site, records) -> {
			double[] quants = new double[quantSize[0]];
			
			for(String[] fields : records) {
				for(int i=0; i<quantSize[0]; i++) {
					String quantValue = fields[i+quantStartIdx[0]];
					if(quantValue.length() == 0) {
						quantValue = "0";
					}
					quants[i] += Double.parseDouble(quantValue);
				}
			}
			
			double sum = 0;
			for(int i=0; i<quants.length; i++) {
				sum += quants[i];
			}
			
			if(sum != 0) {
				try {
					totalSites[0]++;
					BW.append(site);
					for(int i=0; i<quants.length; i++) {
						BW.append("\t"+quants[i]);
					}
					BW.newLine();
				}catch(IOException ioe) {
					System.out.println("Error when calculate quantification method1");
				}
			}
			
			
		});
		
		//record total sites
		stat.sites.add(totalSites[0]);
		BW.close();
	}
	
	// protein id => sequence
	public static Hashtable<String, String> loadProteinFasta (File file) throws IOException {
		BufferedReader BR = new BufferedReader(new FileReader(file));
		Hashtable<String, String> map = new Hashtable<String, String>();
		String line = null;
		
		StringBuilder sequence = new StringBuilder();
		String header = null;
		
		while((line = BR.readLine()) != null) {
			if(line.startsWith(">")) {
				if(header != null) {
					map.put(header, sequence.toString());
				}
				
				sequence.setLength(0);
				header = line.split("\\s")[0].substring(1);
				
			} else {
				sequence.append(line);
			}
		}
		
		map.put(header, sequence.toString());
		
		BR.close();
		
		
		return map;
	}
	
	
	public static Hashtable<String, ArrayList<String[]>> loadPSM(File file, 
			Hashtable<String, String> quanMap,
			Hashtable<String, String> proteinMap) throws IOException {
		BufferedReader BR = new BufferedReader(new FileReader(file));
		String line = null;
		
		Hashtable<String, ArrayList<String[]>> map = new Hashtable<String, ArrayList<String[]>>();
		String[] header = BR.readLine().split("\t");
		
		int annotatedSequenceIdx = -1;
		int fileIdx = -1;
		int scanIdx = -1;
		int quanInfoIdx = -1;
		int phosphoRSIdx = -1;
		int masterProteinAccessionsIdx = -1;
		
		for(int i=0; i<header.length; i++) {
			
			header[i] = header[i].replace("\"", "");
			
			if(header[i].equalsIgnoreCase("Quan Info")) {
				quanInfoIdx = i;
			} else if(header[i].equalsIgnoreCase("Spectrum File")) {
				fileIdx = i;
			} else if(header[i].equalsIgnoreCase("First Scan")) {
				scanIdx = i;
			} else if(header[i].equalsIgnoreCase("PhosphoRS: Best Site Probabilities") ||
					header[i].equalsIgnoreCase("ptmRS: Best Site Probabilities")) {
				phosphoRSIdx = i;
			} else if(header[i].equalsIgnoreCase("Annotated Sequence")) {
				annotatedSequenceIdx = i;
			} else if(header[i].equalsIgnoreCase("Master Protein Accessions")) {
				masterProteinAccessionsIdx = i;
			}
		}
		
		PSM_HEADER = header;
		
		int total = 0;
		int notQuantifiable = 0;
		int shared = 0;
		while((line = BR.readLine()) != null) {
			String[] fields = line.split("\t");
			
			
			for(int i=0; i<fields.length; i++) {
				fields[i] = fields[i].replace("\"", "");
			}
			
			String key = fields[fileIdx]+"_"+fields[scanIdx];
			total++;
			
			// reject "NOT QUANTIFIABLE"
			if(!fields[quanInfoIdx].equalsIgnoreCase("") || quanMap.get(key) == null) {
				notQuantifiable++;
				continue;
			}
			
			
			String[] phosphoSites = fields[phosphoRSIdx].split(";");
			trim(phosphoSites);
			
			// calculate site
			String[] proteins = fields[masterProteinAccessionsIdx].split("\\;");
			trim(proteins);
			
			if(proteins.length > 1) {
				shared++;
				if(!isUseShared) {
					continue;
				}
			}
			
			for(String protein : proteins) {
				String sequence = proteinMap.get(protein);
				if(sequence == null) {
					// Only contaminants will be removed in this step.
					continue;
				} else {
					String peptide = fields[annotatedSequenceIdx].split("\\.")[1].toUpperCase();
					int startIdx = sequence.indexOf(peptide);
					if(startIdx == -1) {
						System.out.println(" OMG! "+peptide+" & " +protein);
					} else {
						
						for(String phosphoSite : phosphoSites) {
							if(phosphoSite.length() == 0 || phosphoSite.equalsIgnoreCase("Too many isoforms")) {
								continue;
							}
							
							
							String site = phosphoSite.split("\\(")[0];
							String aa = site.charAt(0)+"";
							int idx = Integer.parseInt(site.substring(1));
							
							String siteKey = protein+"_"+aa+(idx+startIdx);
							
							// prob threshold
							double prob = Double.parseDouble(phosphoSite.split("\\:")[1].trim());
							if(prob < ptmRSThreshold) {
								continue;
							}
							
							ArrayList<String[]> records = map.get(siteKey);
							if(records == null) {
								records = new ArrayList<String[]>();
								map.put(siteKey, records);
							}
							records.add(fields);
						}
						
					}
				}
			}
			
		}
		
		BR.close();
		
		System.out.println(file.getName());
		System.out.println("Total: "+total);
		System.out.println("Quan rejected: "+notQuantifiable);
		System.out.println("Shared rejected: "+shared);
		System.out.println("Passed sites (this is not a final number): "+map.size());
		System.out.println("Note that passed sites do not consider zero quantifications! the zero quantification will be discarded at the end of process (writing).");
		stat.psms.add(total);
		return map;
	}
	
	// save quantifiable scan
	// key: "file_name"_"scan id"
	public static Hashtable<String, String> loadQuanSpectra(File file) throws IOException {
		BufferedReader BR = new BufferedReader(new FileReader(file));
		String line = null;
		
		Hashtable<String, String> map = new Hashtable<String, String>();
		String[] header = BR.readLine().split("\t");
		
		int fileIdx = -1;
		int scanIdx = -1;
		int quanInfoIdx = -1;
		
		for(int i=0; i<header.length; i++) {
			
			header[i] = header[i].replace("\"", "");
			
			if(header[i].equalsIgnoreCase("Quan Info")) {
				quanInfoIdx = i;
			} else if(header[i].equalsIgnoreCase("Spectrum File")) {
				fileIdx = i;
			} else if(header[i].equalsIgnoreCase("First Scan")) {
				scanIdx = i;
			}
		}
		
		int total = 0;
		int rejected = 0;
		while((line = BR.readLine()) != null) {
			String[] fields = line.split("\t");
			for(int i=0; i<fields.length; i++) {
				fields[i] = fields[i].replace("\"", "");
			}
			total++;
			if(fields[quanInfoIdx].equalsIgnoreCase("NoQuanValues") || fields[quanInfoIdx].equalsIgnoreCase("RejectedByMethod")) {
				rejected++;
				continue;
			}
			
			String key = fields[fileIdx]+"_"+fields[scanIdx];
			
			map.put(key, "");
		}
		
		
		BR.close();
		
		System.out.println(file.getName());
		System.out.println("Total: "+total);
		System.out.println("Rejected: "+rejected);
		System.out.println("Pass: "+map.size());
		double passRate = ((double) map.size())/((double) total);
		System.out.println("Pass rate: "+passRate);
		return map;
	}
	
	public static void trim (String[] data) {
		for(int i=0; i<data.length; i++) {
			data[i] = data[i].trim();
		}
	}
}
