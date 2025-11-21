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

import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
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

import fr.cirad.mgdb.exporting.IExportHandler;
import fr.cirad.mgdb.exporting.tools.ExportManager.ExportOutputs;
import fr.cirad.mgdb.exporting.tools.nj.JukesCantorDistanceMatrixCalculator;
import fr.cirad.mgdb.model.mongo.maintypes.Assembly;
import fr.cirad.tools.ProgressIndicator;
import fr.cirad.tools.mgdb.VariantQueryWrapper;
import fr.cirad.tools.mongo.MongoTemplateManager;
import htsjdk.variant.variantcontext.VariantContext.Type;

/**
 * The Class JukesCantorDistanceMatrixExportHandler.
 */
public class JukesCantorDistanceMatrixExportHandler extends FastaPseudoAlignmentExportHandler {

	private int nMaxMissingDataPercentageForIndividuals = 50;
	
	private final static Pattern missingAllelePattern = Pattern.compile("[N-]");
	
    /**
     * The Constant LOG.
     */
    private static final Logger LOG = Logger.getLogger(JukesCantorDistanceMatrixExportHandler.class);

    /**
     * The supported variant types.
     */
    private static List<String> supportedVariantTypes;

    static {
        supportedVariantTypes = new ArrayList<String>();
        supportedVariantTypes.add(Type.SNP.toString());
    }

    public JukesCantorDistanceMatrixExportHandler() {
    }
    
    public JukesCantorDistanceMatrixExportHandler(int nMaxMissingDataPercentageForIndividuals) {
    	super();
    	this.nMaxMissingDataPercentageForIndividuals = nMaxMissingDataPercentageForIndividuals;
    }
    
    /* (non-Javadoc)
	 * @see fr.cirad.mgdb.exporting.IExportHandler#getExportFormatName()
     */
    @Override
    public String getExportFormatName() {
        return "JK-DIST";
    }

    /* (non-Javadoc)
	 * @see fr.cirad.mgdb.exporting.IExportHandler#getExportFormatDescription()
     */
    @Override
    public String getExportFormatDescription() {
    	return "Exports a zipped distance matrix file calculated according to the <a target='_blank' href='https://www.sglp.uzh.ch/apps/static/MLS/stemmatology/Jukes-Cantor-model_229150204.html'>Jukes-Cantor</a> model, based on a pseudo-alignment consisting in the concatenation of SNP alleles. Individuals with more than 50% missing data are automatically excluded.";
    }

	/* (non-Javadoc)
	 * @see fr.cirad.mgdb.exporting.IExportHandler#getSupportedVariantTypes()
	 */
    @Override
    public List<String> getSupportedVariantTypes() {
        return supportedVariantTypes;
    }
    
	@Override
	public String getExportArchiveExtension() {
		return "zip";
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
	
	        ZipOutputStream zos = IExportHandler.createArchiveOutputStream(outputStream, readyToExportFiles, exportOutputs);
			Assembly assembly = mongoTemplate.findOne(new Query(Criteria.where("_id").is(nAssemblyId)), Assembly.class);
	        
	        ArrayList<String> exportedIndividuals = new ArrayList<>();
	        for (File indFile : exportOutputs.getGenotypeFiles())
	        	try (Scanner scanner = new Scanner(indFile)) {
	        		exportedIndividuals.add(scanner.nextLine());
	        	}
	        
	    	OutputStream fastaOS = new OutputStream() {
	    	    private StringBuilder string = new StringBuilder();
	
	    	    @Override
	    	    public void write(int b) throws IOException {
	    	        this.string.append((char) b );
	    	    }
	
	    	    public String toString() {
	    	        return this.string.toString();
	    	    }
	    	};
	    	fastaOS.write(getHeaderlines(exportOutputs.getGenotypeFiles().length, (int) markerCount).getBytes());
	        writeGenotypeFile(fastaOS, sModule, exportedIndividuals, nQueryChunkSize, markerSynonyms, exportOutputs.getGenotypeFiles(), warningOS, progress);
	        fastaOS.write(getFooterlines().getBytes());
	
	        progress.moveToNextStep();
	        try (Scanner scanner = new Scanner(fastaOS.toString())) {
	            List<String> sequenceNames = new ArrayList<>(), sequences = new ArrayList<>();
	            LinkedHashMap<String, String> seqMap = new LinkedHashMap<>();
	            
	            String currentSequenceName = null;
	            StringBuilder currentSequence = new StringBuilder(), missingDataWarnings = new StringBuilder();
	
	            while (scanner.hasNextLine()) {
	                String line = scanner.nextLine().trim();
	
	                if (line.startsWith(">")) {
	                    // Sequence header line
	                    if (currentSequenceName != null) {
	                    	String sequence = currentSequence.toString();
	                    	Matcher matcher = missingAllelePattern.matcher(sequence);
	                        int missingAlleleCount = 0;
	                        while (matcher.find())
	                            missingAlleleCount++;
	                        if (missingAlleleCount * 100 / sequence.length() > nMaxMissingDataPercentageForIndividuals)
	                        	missingDataWarnings.append("- Excluding individual " + currentSequenceName + " from NJ export, it has too much missing data: " + (missingAlleleCount * 100 / sequence.length()) + "%\n");
	                        else {
		                        sequenceNames.add(currentSequenceName);
		                        sequences.add(sequence);
		                        seqMap.put(currentSequenceName, sequence);
	                        }
	                    }
	                    currentSequenceName = line.substring(1);
	                    currentSequence = new StringBuilder();
	                } else {
	                    // Sequence data line
	                   	currentSequence.append(line);
	                }
	            }
	
	            // Add the last sequence
	            if (currentSequenceName != null) {
	            	String sequence = currentSequence.toString();
	            	Matcher matcher = missingAllelePattern.matcher(sequence);
	                int missingAlleleCount = 0;
	                while (matcher.find())
	                    missingAlleleCount++;
	                if (missingAlleleCount * 100 / sequence.length() > nMaxMissingDataPercentageForIndividuals)
	                	missingDataWarnings.append("- Excluding individual " + currentSequenceName + " from NJ export, it has too much missing data: " + (missingAlleleCount * 100 / sequence.length()) + "%\n");
	                else {
	                    sequenceNames.add(currentSequenceName);
	                    sequences.add(sequence);
	                    seqMap.put(currentSequenceName, sequence);
	                }
	            }
	
	            double[][] distanceMatrix = JukesCantorDistanceMatrixCalculator.calculateDistanceMatrix(sequences, progress);
	            String exportName = IExportHandler.buildExportName(sModule, assembly, markerCount, exportOutputs.getGenotypeFiles().length, exportOutputs.isWorkWithSamples());
	            
	            finalizeExportUsingDistanceMatrix(sequenceNames, exportedIndividuals, exportName, distanceMatrix, zos, progress);

	            if (progress.getError() != null || progress.isAborted())
	        		return;
	
	        	warningOS.close();
		        if (warningFile.length() > 0 || missingDataWarnings.length() > 0) {
		            progress.addStep("Adding lines to warning file");
		            progress.moveToNextStep();
		            progress.setPercentageEnabled(false);
		            zos.putNextEntry(new ZipEntry(exportName + "-REMARKS.txt"));
		            int nWarningCount = 0;
		            if (missingDataWarnings.length() > 0)
		            	zos.write(missingDataWarnings.toString().getBytes());
		            if (warningFile.length() > 0) {
			            BufferedReader in = new BufferedReader(new FileReader(warningFile));
			            String sLine;
			            while ((sLine = in.readLine()) != null) {
			                zos.write((sLine + "\n").getBytes());
			                progress.setCurrentStepProgress(nWarningCount++);
			            }
			            in.close();
		            }
		            LOG.info("Number of Warnings for export (" + exportName + "): " + nWarningCount);
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

	protected void finalizeExportUsingDistanceMatrix(List<String> sequenceNames, ArrayList<String> exportedIndividuals, String exportName, double[][] distanceMatrix, ZipOutputStream zos, ProgressIndicator progress) throws Exception {
        zos.putNextEntry(new ZipEntry(exportName + "." + getExportDataFileExtensions()[0]));
        java.text.NumberFormat formatter = new java.text.DecimalFormat("#0.000000", new DecimalFormatSymbols(Locale.US)); 
        zos.write(("\t" + exportedIndividuals.size()).getBytes());
        int k = 0;
        for (double[] row : distanceMatrix) {
        	zos.write(("\n" + exportedIndividuals.get(k++) + " ").getBytes());
        	for (double cell : row)
        		zos.write((formatter.format(cell) + " ").getBytes());
        }
        zos.closeEntry();
	}

	/* (non-Javadoc)
	 * @see fr.cirad.mgdb.exporting.IExportHandler#getStepList()
     */
    @Override
    public List<String> getStepList() {
        return Arrays.asList(new String[]{"Converting data to FASTA format", "Calculating distance matrix"});
    }

	@Override
	public String[] getExportDataFileExtensions() {
		return new String[] {"mtx"};
	}

    @Override
    public int[] getSupportedPloidyLevels() {
        return new int[] {2};
    }
}