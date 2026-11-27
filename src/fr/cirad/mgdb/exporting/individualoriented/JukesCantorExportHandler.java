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
import fr.cirad.mgdb.exporting.tools.dist.JukesCantorDistanceMatrixCalculator;
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
 * Exports a Jukes-Cantor distance matrix and a neighbour-joining Newick tree, computed from an Eigenstrat data extraction
 * (see JukesCantorDistanceMatrixCalculator for why this is equivalent to working on a pseudo-alignment of SNP alleles).
 */
public class JukesCantorExportHandler extends EigenstratExportHandler {

	private int nMaxMissingDataPercentageForIndividuals = 50;

    private static final Logger LOG = Logger.getLogger(JukesCantorExportHandler.class);

    public static String MATRIX_EXTENSION = "mtx";
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
    public void exportData(OutputStream outputStream, String sModule, Integer nAssemblyId, String sExportingUser, ProgressIndicator progress, String tmpVarCollName, VariantQueryWrapper varQueryWrapper, long markerCount, Map<String, String> markerSynonyms, Map<String, Collection<String>> individualsByPop, boolean workWithSamples, Map<String, HashMap<String, Float>> annotationFieldThresholds, Collection<Callset> callSetsToExport, Collection<String> individualMetadataFieldsToExport, Map<String, InputStream> readyToExportFiles) throws Exception {
        MongoTemplate mongoTemplate = MongoTemplateManager.get(sModule);

        // Build sorted individual list (same order as the columns of the Eigenstrat genotype file)
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

        File warningFile = File.createTempFile("export_warnings_", "");
        File tempGenoFile = File.createTempFile("eigenstrat_for_jc__", "");

        try {
            // Write Eigenstrat genotype file
            ExportOutputs exportOutputs;
            try (OutputStream genoOS = new FileOutputStream(tempGenoFile)) {
            	exportOutputs = writeGenotypeFile(genoOS, sModule, nAssemblyId, workWithSamples, annotationFieldThresholds, progress, tmpVarCollName, variantQueryForTargetCollection, markerCount, markerSynonyms, callSetsToExport, individualsByPop);
            }

            // Collect warnings from export
            try (FileOutputStream warningOS = new FileOutputStream(warningFile)) {
	            for (File f : exportOutputs.getWarningFiles()) {
	                if (f != null && f.length() > 0) {
	                    try (BufferedReader in = new BufferedReader(new FileReader(f))) {
		                    String sLine;
		                    while ((sLine = in.readLine()) != null)
		                        warningOS.write((sLine + "\n").getBytes());
	                    }
	                    f.delete();
	                }
	            }
            }

            if (progress.isAborted() || progress.getError() != null)
            	return;

            progress.moveToNextStep();

            // Read genotypes, leaving out individuals with too much missing data
            StringBuilder missingDataWarnings = new StringBuilder();
            List<String> filteredIndividuals = new ArrayList<>();
            byte[][] genotypeBytes = readGenotypeDataAndFilterByMissingness(tempGenoFile, sortedIndividuals, (int) markerCount, nMaxMissingDataPercentageForIndividuals, filteredIndividuals, missingDataWarnings, progress, "Jukes-Cantor");
            if (progress.isAborted() || progress.getError() != null)
            	return;

            if (genotypeBytes == null || filteredIndividuals.size() < 2)
                throw new Exception("Not enough individuals with sufficient data for Jukes-Cantor calculation.");

            progress.moveToNextStep();

            // Compute distance matrix (upper triangle)
            double[][] distanceMatrix = new JukesCantorDistanceMatrixCalculator(genotypeBytes, (int) markerCount).calculate(progress);
            if (distanceMatrix == null || progress.isAborted() || progress.getError() != null)
                return;

			Assembly assembly = mongoTemplate.findOne(new Query(Criteria.where("_id").is(nAssemblyId)), Assembly.class);
	        String exportName = IExportHandler.buildExportName(sModule, assembly, markerCount, sortedIndividuals.size(), workWithSamples);
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
	    	String refPosPathWithTrailingDot = Assembly.getThreadBoundVariantRefPosPath() + ".";
	    	Document projectionAndSortDoc = new Document(refPosPathWithTrailingDot + ReferencePosition.FIELDNAME_SEQUENCE, 1).append(refPosPathWithTrailingDot + ReferencePosition.FIELDNAME_START_SITE, 1);
	    	try (MongoCursor<Document> markerCursor = IExportHandler.getMarkerCursorWithCorrectCollation(mongoTemplate.getCollection(tmpVarCollName != null ? tmpVarCollName : mongoTemplate.getCollectionName(VariantData.class)), variantQueryForTargetCollection, projectionAndSortDoc, IExportHandler.computeQueryChunkSize(mongoTemplate, markerCount))) {
	            progress.addStep("Writing map file");
	            progress.moveToNextStep();
		        while (markerCursor.hasNext()) {
		            Document exportVariant = markerCursor.next();
		            Document refPos = (Document) Helper.readPossiblyNestedField(exportVariant, refPosPath, ";", null);
		            Long pos = refPos == null || refPos.isEmpty() ? null : ((Number) refPos.get(ReferencePosition.FIELDNAME_START_SITE)).longValue();
		            String chrom = refPos == null || refPos.isEmpty() ? null : (String) refPos.get(ReferencePosition.FIELDNAME_SEQUENCE);
	                String markerId = (String) exportVariant.get("_id");
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
            LOG.error("Not enough RAM for Jukes-Cantor calculation", oome);
            progress.setError("Unable to compute Jukes-Cantor distances (dataset is too large): " + oome.getMessage());
        }
        finally {
	        warningFile.delete();
	        tempGenoFile.delete();
	    }

	    progress.setPercentageEnabled(true);
	    progress.setCurrentStepProgress((short) 100);
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
        nodes[nodes.length - 1].buildTreeString(sb);

        zos.putNextEntry(new ZipEntry(exportName + "." + NEWICK_EXTENSION));
        zos.write((sb.toString() + ";").getBytes());
        zos.closeEntry();
    }

	/* (non-Javadoc)
	 * @see fr.cirad.mgdb.exporting.IExportHandler#getStepList()
     */
    @Override
    public List<String> getStepList() {
    	return Arrays.asList(new String[]{"Exporting to Eigenstrat format", "Filtering individuals with missing data", "Calculating Jukes-Cantor distance matrix", "Writing matrix to archive", "Building neighbor-joining Newick tree using NINJA algorithm"});
    }

	@Override
	public String[] getExportDataFileExtensions() {
		return new String[] {MATRIX_EXTENSION, NEWICK_EXTENSION};
	}
}
