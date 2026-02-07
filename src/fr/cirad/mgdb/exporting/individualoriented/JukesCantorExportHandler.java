/*******************************************************************************
 * MGDB Export - Mongo Genotype DataBase, export handlers
 * Copyright (C) 2016 - 2019, <CIRAD> <IRD>
 *
 * This program is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Affero General Public License, version 3 as published by
 * the Free Software Foundation.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more
 * details.
 *
 * See <http://www.gnu.org/licenses/agpl.html> for details about GNU General
 * Public License V3.
 *******************************************************************************/
package fr.cirad.mgdb.exporting.individualoriented;

import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.InputStream;
import java.io.OutputStream;
import java.text.DecimalFormatSymbols;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Scanner;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;

import org.apache.log4j.Logger;
import org.springframework.data.mongodb.core.MongoTemplate;
import org.springframework.data.mongodb.core.query.Criteria;
import org.springframework.data.mongodb.core.query.Query;

import com.traviswheeler.ninja.TreeBuilderBinHeap;
import com.traviswheeler.ninja.TreeNode;

import fr.cirad.mgdb.exporting.IExportHandler;
import fr.cirad.mgdb.exporting.tools.ExportManager.ExportOutputs;
import fr.cirad.mgdb.exporting.tools.dist.JukesCantorDistanceMatrixCalculator;
import fr.cirad.mgdb.model.mongo.maintypes.Assembly;
import fr.cirad.tools.ProgressIndicator;
import fr.cirad.tools.mgdb.VariantQueryWrapper;
import fr.cirad.tools.mongo.MongoTemplateManager;

/**
 * The Class JukesCantorNewickTreeExportHandler.
 */
public class JukesCantorExportHandler extends FastaPseudoAlignmentExportHandler {

	private int nMaxMissingDataPercentageForIndividuals = 50;
	private final static Pattern missingAllelePattern = Pattern.compile("[N-]");
	
    private static final Logger LOG = Logger.getLogger(JukesCantorExportHandler.class);

    public static String MATRIX_EXTENSION= "mtx";
    public static String NEWICK_EXTENSION = "nwk";
    
    public JukesCantorExportHandler() {
    }
    
    public JukesCantorExportHandler(int nMaxMissingDataPercentageForIndividuals) {
    	super();
    	this.nMaxMissingDataPercentageForIndividuals = nMaxMissingDataPercentageForIndividuals;
    }
    
    /* (non-Javadoc)
	 * @see fr.cirad.mgdb.exporting.IExportHandler#getExportFormatName()
     */
    @Override
    public String getExportFormatName() {
        return "JUKES-CANTOR";
    }

    /* (non-Javadoc)
	 * @see fr.cirad.mgdb.exporting.IExportHandler#getExportFormatDescription()
     */
    @Override
    public String getExportFormatDescription() {
    	return "Exports a zip archive featuring a <a target='_blank' href='https://treethinkers.org/jukes-cantor-model-of-dna-substitution/'>Jukes-Cantor</a> distance matrix and a <a target='_blank' href='https://en.wikipedia.org/wiki/Newick_format'>Newick</a> file containing a neighbour-joining tree, both based on a pseudo-alignment consisting in the concatenation of SNP alleles. Individuals with more than 50% missing data are automatically excluded.";
    }

    @Override
    public void exportData(OutputStream outputStream, String sModule, Integer nAssemblyId, ExportOutputs exportOutputs, boolean fDeleteSampleExportFilesOnExit, ProgressIndicator progress, String tmpVarCollName, VariantQueryWrapper varQueryWrapper, long markerCount, Map<String, String> markerSynonyms, Map<String, String> individualPopulations, Map<String, InputStream> readyToExportFiles) throws Exception {
		MongoTemplate mongoTemplate = MongoTemplateManager.get(sModule);
        int nQueryChunkSize = IExportHandler.computeQueryChunkSize(mongoTemplate, markerCount);

        // save existing warnings into a temp file so we can append to it
        File warningFile = File.createTempFile("export_warnings_", "");
        try {
	        FileOutputStream warningOS = new FileOutputStream(warningFile);
	        for (File f : exportOutputs.getWarningFiles()) {
		    	if (f != null && f.length() > 0) {
		            BufferedReader in = new BufferedReader(new FileReader(f));
		            String sLine;
		            while ((sLine = in.readLine()) != null)
		            	warningOS.write((sLine + "\n").getBytes());
		            in.close();
			    	f.delete();
		    	}
	        }
		        
	        List<String> exportedIndividuals = new ArrayList<>();
	        for (File indFile : exportOutputs.getGenotypeFiles())
	            try (Scanner scanner = new Scanner(indFile)) {
	                exportedIndividuals.add(scanner.nextLine());
	            }

	        // Stream sequences to temporary FASTA file
	        File tmpFastaFile = File.createTempFile("fasta_", ".tmp");
	        tmpFastaFile.deleteOnExit();
	        try (BufferedOutputStream fastaOS = new BufferedOutputStream(new FileOutputStream(tmpFastaFile))) {
	            fastaOS.write(getHeaderlines(exportOutputs.getGenotypeFiles().length, (int) markerCount).getBytes());
	            writeGenotypeFile(fastaOS, sModule, exportedIndividuals, nQueryChunkSize, markerSynonyms, exportOutputs.getGenotypeFiles(), warningOS, progress);
	            fastaOS.write(getFooterlines().getBytes());
	        }

	        progress.moveToNextStep();

	        // Process sequences line by line from temp file
	        List<String> sequenceNames = new ArrayList<>();
	        List<String> sequences = new ArrayList<>();
	        LinkedHashMap<String, String> seqMap = new LinkedHashMap<>();
	        StringBuilder missingDataWarnings = new StringBuilder();

	        try (BufferedReader reader = new BufferedReader(new FileReader(tmpFastaFile))) {
	            String line;
	            String currentName = null;
	            StringBuilder currentSeq = new StringBuilder();
	            while ((line = reader.readLine()) != null) {
	                line = line.trim();
	                if (line.startsWith(">")) {
	                    if (currentName != null) {
	                        processSequence(currentName, currentSeq.toString(), sequenceNames, sequences, seqMap, missingDataWarnings);
	                    }
	                    currentName = line.substring(1);
	                    currentSeq.setLength(0);
	                } else {
	                    currentSeq.append(line);
	                }
	            }
	            if (currentName != null) {
	                processSequence(currentName, currentSeq.toString(), sequenceNames, sequences, seqMap, missingDataWarnings);
	            }
	        }

	        // Compute distance matrix (upper-triangle safe)
	        double[][] distanceMatrix = new JukesCantorDistanceMatrixCalculator(sequences).calculate(progress);

			Assembly assembly = mongoTemplate.findOne(new Query(Criteria.where("_id").is(nAssemblyId)), Assembly.class);
	        String exportName = IExportHandler.buildExportName(sModule, assembly, markerCount, exportOutputs.getGenotypeFiles().length, exportOutputs.isWorkWithSamples());
	        ZipOutputStream zos = IExportHandler.createArchiveOutputStream(outputStream, readyToExportFiles, exportOutputs);
	        createMatrixEntry(sequenceNames, exportName, distanceMatrix, zos, progress);
	        createTreeEntry(sequenceNames, exportName, distanceMatrix, zos, progress);

	        // Add warnings if any
	        if (progress.getError() == null && !progress.isAborted()) {
	            if (warningFile.length() > 0 || missingDataWarnings.length() > 0) {
	                progress.addStep("Adding lines to warning file");
	                progress.moveToNextStep();
	                progress.setPercentageEnabled(false);
	                zos.putNextEntry(new ZipEntry(exportName + "-REMARKS.txt"));

	                if (missingDataWarnings.length() > 0)
	                    zos.write(missingDataWarnings.toString().getBytes());

	                if (warningFile.length() > 0) {
	                    try (BufferedReader in = new BufferedReader(new FileReader(warningFile))) {
	                        String sLine;
	                        int nWarningCount = 0;
	                        while ((sLine = in.readLine()) != null) {
	                            zos.write((sLine + "\n").getBytes());
	                            progress.setCurrentStepProgress(nWarningCount++);
	                        }
	                        LOG.info("Number of Warnings for export (" + exportName + "): " + nWarningCount);
	                    }
	                }

	                zos.closeEntry();
	            }
	        }

	        zos.finish();
	        zos.close();
	    }
        finally {
	        warningFile.delete();
	    }

	    progress.setPercentageEnabled(true);
	    progress.setCurrentStepProgress((short) 100);
	}
    
	private void processSequence(String name, String sequence, List<String> sequenceNames, List<String> sequences, Map<String, String> seqMap, StringBuilder missingDataWarnings) {
		Matcher matcher = missingAllelePattern.matcher(sequence);
		int missingAlleleCount = 0;
		while (matcher.find())
			missingAlleleCount++;

		int missingPercent = missingAlleleCount * 100 / sequence.length();
		if (missingPercent > nMaxMissingDataPercentageForIndividuals)
			missingDataWarnings.append("- Excluding individual ").append(name).append(" from NJ export, it has too much missing data: ").append(missingPercent).append("%\n");
		else {
			sequenceNames.add(name);
			sequences.add(sequence);
			seqMap.put(name, sequence);
		}
	}

	protected void createMatrixEntry(List<String> sequenceNames, String exportName, double[][] distanceMatrix, ZipOutputStream zos, ProgressIndicator progress) throws Exception {
	    zos.putNextEntry(new ZipEntry(exportName + "." + MATRIX_EXTENSION));
	    java.text.NumberFormat formatter = new java.text.DecimalFormat("#0.000000", new DecimalFormatSymbols(Locale.US)); 
	    
	    int n = sequenceNames.size();	    
	    zos.write(("\t" + n).getBytes());	// Write header with number of individuals
	    
	    // Write complete distance matrix (not just upper triangle)
	    for (int i = 0; i < n; i++) {
	        StringBuilder rowBuilder = new StringBuilder();
	        rowBuilder.append("\n").append(sequenceNames.get(i));
	        
	        for (int j = 0; j < n; j++) {
	            double distance;
	            if (i == j)
	                distance = 0.0;  // Diagonal: distance to self is 0
	            else if (i < j)
	            	distance = distanceMatrix[i][j - i - 1];	// Upper triangle: retrieve from distanceMatrix[i][j-i-1]
	            else
	                distance = distanceMatrix[j][i - j - 1];	// Lower triangle: symmetric to upper triangle
	            rowBuilder.append(" ").append(formatter.format(distance));
	        }
	        
	        zos.write(rowBuilder.toString().getBytes());
	        
	        if (n > 1000 && i % 500 == 0 && progress != null)	// Update progress for large datasets
	            progress.setCurrentStepProgress((int) ((i / (double) n) * 100));
	    }
	    
	    zos.closeEntry();
	}

    protected void createTreeEntry(List<String> sequenceNames, String exportName, double[][] distanceMatrix, ZipOutputStream zos, ProgressIndicator progress) throws Exception {
        progress.moveToNextStep();

        final int n = distanceMatrix.length;

        // Upper-triangle int distance matrix
        int[][] intDistanceMatrix = new int[n][];
        for (int i = 0; i < n; i++) {
            intDistanceMatrix[i] = new int[n - 1 - i];
            for (int j = i + 1; j < n; j++)
                intDistanceMatrix[i][j - i - 1] = (int) (distanceMatrix[i][j - i - 1] * 100000000);
        }

        StringBuffer sb = new StringBuffer();

        TreeBuilderBinHeap tb = new TreeBuilderBinHeap(sequenceNames.toArray(new String[sequenceNames.size()]), intDistanceMatrix);

        TreeNode[] nodes = tb.build();

//          nodes[nodes.length-1].rootTreeAt("CR1062");
        nodes[nodes.length - 1].buildTreeString(sb);

        zos.putNextEntry(new ZipEntry(exportName + "." + NEWICK_EXTENSION));
        String treeString = sb.toString() + ";";
        zos.write(treeString.getBytes());

//          NewickTreeRerooter rerooter = new NewickTreeRerooter();
//          TreeNode root = rerooter.parseNewick(treeString);
//          TreeNode newRoot = rerooter.reroot(root, "CX280");
//          String newNewick = rerooter.toNewick(newRoot);
//          zos.write(newNewick.getBytes());

//          zos.write(NewickRooter.rootTree(treeString, /* \"IRIS_313-10729\" */ \"CX280\").getBytes());
//          StringBuffer rootedTreeSB = new StringBuffer();
//          rootTree(treeString, /* \"IRIS_313-10729\" */ \"CX280\").buildTreeString(rootedTreeSB);
//          zos.write(rootedTreeSB.toString().getBytes());

        zos.closeEntry();
    }

	/* (non-Javadoc)
	 * @see fr.cirad.mgdb.exporting.IExportHandler#getStepList()
     */
    @Override
    public List<String> getStepList() {
    	return Arrays.asList(new String[]{"Converting data to FASTA format", "Calculating Jukes-Cantor distance matrix", "Building neighbor-joining Newick tree using NINJA algorithm"});
    }

	@Override
	public String[] getExportDataFileExtensions() {
		return new String[] {MATRIX_EXTENSION, NEWICK_EXTENSION};
	}
}