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
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.TreeSet;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;

import org.apache.log4j.Logger;
import org.bson.Document;
import org.springframework.data.mongodb.core.MongoTemplate;
import org.springframework.data.mongodb.core.query.Criteria;
import org.springframework.data.mongodb.core.query.Query;

import com.mongodb.BasicDBList;

import com.mongodb.client.MongoCursor;
import com.traviswheeler.ninja.TreeBuilderBinHeap;
import com.traviswheeler.ninja.TreeNode;

import fr.cirad.mgdb.exporting.IExportHandler;
import fr.cirad.mgdb.exporting.markeroriented.EigenstratExportHandler;
import fr.cirad.mgdb.exporting.tools.ExportManager.ExportOutputs;
import fr.cirad.mgdb.exporting.tools.dist.AlleleSharingDistanceMatrixCalculator;
import fr.cirad.mgdb.model.mongo.maintypes.Assembly;
import fr.cirad.mgdb.model.mongo.maintypes.VariantData;
import fr.cirad.mgdb.model.mongo.subtypes.Callset;
import fr.cirad.mgdb.model.mongo.subtypes.ReferencePosition;
import fr.cirad.tools.AlphaNumericComparator;
import fr.cirad.tools.Helper;
import fr.cirad.tools.ProgressIndicator;
import fr.cirad.tools.mgdb.VariantQueryWrapper;
import fr.cirad.tools.mongo.MongoTemplateManager;

/**
 * The Class ASDNewickTreeExportHandler.
 */
public class AlleleSharingDistanceExportHandler extends EigenstratExportHandler {

    private int nMaxMissingDataPercentageForIndividuals = 50;

    private static final Logger LOG = Logger.getLogger(AlleleSharingDistanceExportHandler.class);
    
    public AlleleSharingDistanceExportHandler() {
    }
    
    public AlleleSharingDistanceExportHandler(int nMaxMissingDataPercentageForIndividuals) {
        super();
        this.nMaxMissingDataPercentageForIndividuals = nMaxMissingDataPercentageForIndividuals;
    }
    
    /* (non-Javadoc)
	 * @see fr.cirad.mgdb.exporting.IExportHandler#getExportFormatName()
     */
    @Override
    public String getExportFormatName() {
        return "ASD";
    }

    /* (non-Javadoc)
	 * @see fr.cirad.mgdb.exporting.IExportHandler#getExportFormatDescription()
     */
    @Override
    public String getExportFormatDescription() {
    	return "Exports a zip archive featuring an Allele Sharing Distance (ASD) matrix and a <a target='_blank' href='https://en.wikipedia.org/wiki/Newick_format'>Newick</a> file containing a neighbour-joining tree, both based on an Eigenstrat data extraction. Individuals with more than 50% missing data are automatically excluded.";
    }

	protected void createMatrixEntry(List<String> sequenceNames, String exportName, double[][] distanceMatrix, ZipOutputStream zos, ProgressIndicator progress) throws Exception {
        zos.putNextEntry(new ZipEntry(exportName + "." + JukesCantorExportHandler.MATRIX_EXTENSION));
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

    protected void createTreeEntry(List<String> sequenceNames, String exportName, double[][] distanceMatrix, ZipOutputStream zos) throws Exception {
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

        zos.putNextEntry(new ZipEntry(exportName + "." + JukesCantorExportHandler.NEWICK_EXTENSION));
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

    public static byte[][] readGenotypeDataAndFilterByMissingness(
            File genoFile,
            List<String> allIndividuals,
            int numMarkers,
            int maxMissingPct,
            List<String> keptIndividuals,
            StringBuilder missingDataWarnings,
            ProgressIndicator progress) throws IOException {

        final int nInd = allIndividuals.size();

        byte[][] genotypeAll = new byte[nInd][numMarkers];
        int[] missingCounts = new int[nInd];

        try (BufferedReader reader = new BufferedReader(new FileReader(genoFile))) {
            String line;
            int markerIndex = 0;

            while ((line = reader.readLine()) != null && markerIndex < numMarkers) {
                line = line.trim();
                if (line.length() < nInd)
                    throw new IllegalStateException("Line " + markerIndex + " too short: expected " + nInd + " chars");

                for (int i = 0; i < nInd; i++) {
                    char c = line.charAt(i);
                    if (c >= '0' && c <= '2')
                        genotypeAll[i][markerIndex] = (byte) (c - '0');
                    else {
                        genotypeAll[i][markerIndex] = 9;
                        missingCounts[i]++;
                    }
                }

                markerIndex++;

                if (progress != null && markerIndex % 100 == 0) {
                    int pct = (int) ((markerIndex / (double) numMarkers) * 100);
                    progress.setCurrentStepProgress(pct);

                    if (progress.isAborted() || progress.getError() != null)
                        return null;
                }
            }

            if (markerIndex != numMarkers)
                throw new IllegalStateException("Expected " + numMarkers + " markers, got " + markerIndex);
        }

        // Decide which individuals to keep
        List<Integer> keptIndices = new ArrayList<>();
        for (int i = 0; i < nInd; i++) {
            String ind = allIndividuals.get(i);
            int missingPct = (int) ((missingCounts[i] * 100.0) / numMarkers);

            if (missingPct > maxMissingPct)
                missingDataWarnings.append("- Excluding individual ").append(ind).append(" from ASD export, it has too much missing data: ").append(missingPct).append("%\n");
            else {
                keptIndices.add(i);
                keptIndividuals.add(ind);
            }
        }

        if (keptIndices.size() < 2)
            return null;

        // Compact genotype matrix
        byte[][] genotypeFiltered = new byte[keptIndices.size()][numMarkers];
        for (int fi = 0; fi < keptIndices.size(); fi++) {
        	Integer index = keptIndices.get(fi);
            genotypeFiltered[fi] = genotypeAll[index];
            genotypeAll[index] = null;	// save memory
        }

        return genotypeFiltered;
    }

    @Override
    public void exportData(OutputStream outputStream, String sModule, Integer nAssemblyId, String sExportingUser, ProgressIndicator progress, String tmpVarCollName, VariantQueryWrapper varQueryWrapper, long markerCount, Map<String, String> markerSynonyms, Map<String, Collection<String>> individualsByPop, boolean workWithSamples, Map<String, HashMap<String, Float>> annotationFieldThresholds, Collection<Callset> callSetsToExport, Collection<String> individualMetadataFieldsToExport, Map<String, InputStream> readyToExportFiles) throws Exception {
        MongoTemplate mongoTemplate = MongoTemplateManager.get(sModule);
        
        // Build sorted individual list
        TreeSet<String> indSet = new TreeSet<>(new AlphaNumericComparator<String>());
        for (Callset cs : callSetsToExport)
            indSet.add(workWithSamples ? cs.getSampleId() : cs.getIndividual());
        List<String> sortedIndividuals = new ArrayList<>(indSet);

        if (progress.isAborted())
            return;

        // Build variant query
        Collection<BasicDBList> variantDataQueries = varQueryWrapper.getVariantDataQueries();
        Document variantQueryForTargetCollection = variantDataQueries.isEmpty() ? new Document() : 
            (tmpVarCollName == null ? new Document("$and", variantDataQueries.iterator().next()) : 
                (varQueryWrapper.getBareQueries().iterator().hasNext() ? 
                    new Document("$and", varQueryWrapper.getBareQueries().iterator().next()) : new Document()));

        // Save warnings to temp file
        File warningFile = File.createTempFile("export_warnings_", "");
        File tempGenoFile = File.createTempFile("eigenstrat_for_asd__", "");
        
        try {
            // Write Eigenstrat genotype file
            OutputStream genoOS = new FileOutputStream(tempGenoFile);
            ExportOutputs exportOutputs = writeGenotypeFile(genoOS, sModule, nAssemblyId, workWithSamples, annotationFieldThresholds, progress, tmpVarCollName, variantQueryForTargetCollection, markerCount, markerSynonyms, callSetsToExport, individualsByPop);
            genoOS.close();

            // Collect warnings from export
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
            warningOS.close();

            progress.moveToNextStep();

            // Process genotype data from temp file
            StringBuilder missingDataWarnings = new StringBuilder();
            List<String> filteredIndividuals = new ArrayList<>();
            byte[][] genotypeBytes = readGenotypeDataAndFilterByMissingness(tempGenoFile, sortedIndividuals, (int) markerCount, nMaxMissingDataPercentageForIndividuals, filteredIndividuals, missingDataWarnings, progress);

            if (genotypeBytes == null || filteredIndividuals.size() < 2)
                throw new Exception("Not enough individuals with sufficient data for ASD calculation.");

            progress.moveToNextStep();

            // Compute ASD distance matrix
            double[][] distanceMatrix = new AlleleSharingDistanceMatrixCalculator(genotypeBytes/*, filteredIndividuals.size()*/, (int) markerCount).calculate(progress);
            if (distanceMatrix == null || progress.isAborted() || progress.getError() != null)
                return;
            
            // Create export name
            Assembly assembly = mongoTemplate.findOne(new Query(Criteria.where("_id").is(nAssemblyId)), Assembly.class);
            String exportName = sModule + (assembly != null && assembly.getName() != null ? "__" + assembly.getName() : "") + "__" + markerCount + "variants__" + filteredIndividuals.size() + "individuals";
            ZipOutputStream zos = IExportHandler.createArchiveOutputStream(outputStream, readyToExportFiles, null);
            
            // Add metadata if requested
            if (individualMetadataFieldsToExport == null || !individualMetadataFieldsToExport.isEmpty())
                IExportHandler.addMetadataEntryIfAny(sModule + "__" + filteredIndividuals.size() + (workWithSamples ? "sample" : "individual") + "s_metadata.tsv", sModule, sExportingUser, new TreeSet<>(filteredIndividuals), individualMetadataFieldsToExport, zos, (workWithSamples ? "sample" : "individual"), workWithSamples);

            progress.moveToNextStep();
            createMatrixEntry(filteredIndividuals, exportName, distanceMatrix, zos, progress);
            progress.moveToNextStep();
            createTreeEntry(filteredIndividuals, exportName, distanceMatrix, zos);

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
            
	        zos.putNextEntry(new ZipEntry(exportName + ".map"));
	        String refPosPath = Assembly.getVariantRefPosPath(nAssemblyId);
	        int nMarkerIndex = 0;
	        ArrayList<Comparable> unassignedMarkers = new ArrayList<>();
	    	String refPosPathWithTrailingDot = Assembly.getThreadBoundVariantRefPosPath() + ".";
	    	Document projectionAndSortDoc = new Document(refPosPathWithTrailingDot + ReferencePosition.FIELDNAME_SEQUENCE, 1).append(refPosPathWithTrailingDot + ReferencePosition.FIELDNAME_START_SITE, 1);
	    	try (MongoCursor<Document> markerCursor = IExportHandler.getMarkerCursorWithCorrectCollation(mongoTemplate.getCollection(tmpVarCollName != null ? tmpVarCollName : mongoTemplate.getCollectionName(VariantData.class)), variantQueryForTargetCollection, projectionAndSortDoc, 10000)) {
	            progress.addStep("Writing map file");
	            progress.moveToNextStep();
		        while (markerCursor.hasNext()) {
		            Document exportVariant = markerCursor.next();
		            Document refPos = (Document) Helper.readPossiblyNestedField(exportVariant, refPosPath, ";", null);
		            Long pos = refPos == null ? null : ((Number) refPos.get(ReferencePosition.FIELDNAME_START_SITE)).longValue();
		            String chrom = refPos == null ? null : (String) refPos.get(ReferencePosition.FIELDNAME_SEQUENCE);
	                String markerId = (String) exportVariant.get("_id");
		            if (chrom == null)
		            	unassignedMarkers.add(markerId);
		            String exportedId = markerSynonyms == null ? markerId : markerSynonyms.get(markerId);
		            zos.write(((chrom == null ? "0" : chrom) + " " + exportedId + " " + 0 + " " + (pos == null ? 0 : pos) + LINE_SEPARATOR).getBytes());
	
	                progress.setCurrentStepProgress(nMarkerIndex++ * 100 / markerCount);
		        }
			}
	        zos.closeEntry();
	        
	        IExportHandler.writeZipEntryFromChunkFiles(zos, exportOutputs.getAnnotationFiles(), exportName + ".ann", IExportHandler.VEP_LIKE_HEADER_LINE);

            zos.finish();
            zos.close();
        }
        catch (OutOfMemoryError oome) {
            LOG.error("Not enough RAM for ASD calculation", oome);
            progress.setError("Unable to compute ASD (dataset is too large): " + oome.getMessage());
        }
        finally {
            warningFile.delete();
            tempGenoFile.delete();
        }

        progress.setPercentageEnabled(true);
        progress.setCurrentStepProgress((short) 100);
    }
    
    /**
     * Reads genotype data from a file into byte arrays, producing a matrix of shape [nIndividuals][nMarkers].
     * Supports filtering individuals by name.
     *
     * File format:
     * - Each row = a SNP / marker
     * - Each column = an individual
     * - Characters = '0', '1', '2', '9' (missing)
     *
     * @param genoFile          input file
     * @param allIndividuals    list of all individual names in file order
     * @param filteredIndividuals subset to keep (may be same as allIndividuals)
     * @param numMarkers        number of markers / rows in the file
     * @param progress          progress indicator
     * @return [nFilteredIndividuals][numMarkers] numeric byte matrix
     */
    public static byte[][] readFilteredGenotypeData(
            File genoFile,
            List<String> allIndividuals,
            List<String> filteredIndividuals,
            int numMarkers,
            ProgressIndicator progress) throws IOException {

        int nAll = allIndividuals.size();
        int nFiltered = filteredIndividuals.size();

        // Map individual name → column index
        Map<String, Integer> individualIndexMap = new HashMap<>();
        for (int i = 0; i < nAll; i++) {
            individualIndexMap.put(allIndividuals.get(i), i);
        }

        // Prepare filtered indices (columns)
        int[] filteredIndices = new int[nFiltered];
        for (int i = 0; i < nFiltered; i++) {
            filteredIndices[i] = individualIndexMap.get(filteredIndividuals.get(i));
        }

        // Allocate output: [nFilteredIndividuals][numMarkers]
        byte[][] genotypeBytes = new byte[nFiltered][numMarkers];

        try (BufferedReader reader = new BufferedReader(new FileReader(genoFile))) {
            String line;
            int markerIndex = 0;

            while ((line = reader.readLine()) != null && markerIndex < numMarkers) {
                line = line.trim();
                if (line.length() < nAll) {
                    throw new IllegalStateException("Line " + markerIndex + " too short: expected " + nAll + " chars");
                }

                // Copy only filtered individuals
                for (int fi = 0; fi < nFiltered; fi++) {
                    int col = filteredIndices[fi];
                    char c = line.charAt(col);
                    if (c >= '0' && c <= '2') {
                        genotypeBytes[fi][markerIndex] = (byte) (c - '0');
                    } else {
                        genotypeBytes[fi][markerIndex] = 9; // missing
                    }
                }

                markerIndex++;

                // Progress update
                if (progress != null && markerIndex % 100 == 0) {
                    int pct = 50 + (int) ((markerIndex / (double) numMarkers) * 25);
                    progress.setCurrentStepProgress(pct);

                    if (progress.isAborted() || progress.getError() != null) {
                        System.out.println("ASD export aborted during data reading");
                        return null;
                    }
                }
            }

            if (markerIndex != numMarkers) {
                throw new IllegalStateException("Expected " + numMarkers + " markers, but got " + markerIndex);
            }
        }

        return genotypeBytes;
    }

	/* (non-Javadoc)
	 * @see fr.cirad.mgdb.exporting.IExportHandler#getStepList()
     */
    @Override
    public List<String> getStepList() {
    	return Arrays.asList(new String[]{"Exporting to Eigenstrat format",  "Filtering individuals with missing data", "Calculating Allele Sharing Distance matrix", "Writing matrix to archive", "Building neighbor-joining Newick tree using NINJA algorithm"});
    }

	@Override
	public String[] getExportDataFileExtensions() {
		return new String[] {JukesCantorExportHandler.MATRIX_EXTENSION, JukesCantorExportHandler.NEWICK_EXTENSION};
	}
}