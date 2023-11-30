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
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Scanner;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;

import org.apache.commons.lang.StringUtils;
import org.apache.log4j.Logger;
import org.bson.Document;
import org.springframework.data.mongodb.core.MongoTemplate;
import org.springframework.data.mongodb.core.query.Criteria;
import org.springframework.data.mongodb.core.query.Query;

import com.mongodb.BasicDBList;
import com.mongodb.BasicDBObject;
import com.mongodb.client.MongoCursor;

import fr.cirad.mgdb.exporting.IExportHandler;
import fr.cirad.mgdb.model.mongo.maintypes.Assembly;
import fr.cirad.mgdb.model.mongo.maintypes.VariantData;
import fr.cirad.mgdb.model.mongo.subtypes.ReferencePosition;
import fr.cirad.mgdb.model.mongodao.MgdbDao;
import fr.cirad.tools.Helper;
import fr.cirad.tools.ProgressIndicator;
import fr.cirad.tools.mgdb.VariantQueryWrapper;
import fr.cirad.tools.mongo.MongoTemplateManager;
import htsjdk.variant.variantcontext.VariantContext.Type;

/**
 * The Class PLinkExportHandler.
 */
public class PLinkExportHandler extends AbstractIndividualOrientedExportHandler {

    /**
     * The Constant LOG.
     */
    private static final Logger LOG = Logger.getLogger(PLinkExportHandler.class);

    /**
     * The supported variant types.
     */
    private static List<String> supportedVariantTypes;

    static {
        supportedVariantTypes = new ArrayList<String>();
        supportedVariantTypes.add(Type.SNP.toString());
        supportedVariantTypes.add(Type.MNP.toString());
        supportedVariantTypes.add(Type.INDEL.toString());
        supportedVariantTypes.add(Type.MIXED.toString());
        supportedVariantTypes.add(Type.NO_VARIATION.toString());
    }

    /* (non-Javadoc)
	 * @see fr.cirad.mgdb.exporting.IExportHandler#getExportFormatName()
     */
    @Override
    public String getExportFormatName() {
        return "PLINK";
    }

    /* (non-Javadoc)
	 * @see fr.cirad.mgdb.exporting.IExportHandler#getExportFormatDescription()
     */
    @Override
    public String getExportFormatDescription() {
    	return "Exports zipped PED and MAP files. See <a target='_blank' href='http://zzz.bwh.harvard.edu/plink/data.shtml#ped'>http://zzz.bwh.harvard.edu/plink/data.shtml#ped</a> for more details";
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
    public void exportData(OutputStream outputStream, String sModule, Integer nAssemblyId, String sExportingUser, File[] individualExportFiles, boolean fDeleteSampleExportFilesOnExit, ProgressIndicator progress, String tmpVarCollName, VariantQueryWrapper varQueryWrapper, long markerCount, Map<String, String> markerSynonyms, Collection<String> individualMetadataFieldsToExport, Map<String, String> individualPopulations, Map<String, InputStream> readyToExportFiles) throws Exception {
		MongoTemplate mongoTemplate = MongoTemplateManager.get(sModule);
        int nQueryChunkSize = IExportHandler.computeQueryChunkSize(mongoTemplate, markerCount);
        File warningFile = File.createTempFile("export_warnings_", "");
        FileWriter warningFileWriter = new FileWriter(warningFile);
        ZipOutputStream zos = IExportHandler.createArchiveOutputStream(outputStream, readyToExportFiles);
		Assembly assembly = mongoTemplate.findOne(new Query(Criteria.where("_id").is(nAssemblyId)), Assembly.class);
        String exportName = sModule + (assembly != null && assembly.getName() != null ? "__" + assembly.getName() : "") + "__" + markerCount + "variants__" + individualExportFiles.length + "individuals";
        
        ArrayList<String> exportedIndividuals = new ArrayList<>();
        for (File indFile : individualExportFiles)
        	try (Scanner scanner = new Scanner(indFile)) {
        		exportedIndividuals.add(scanner.nextLine());
        	}

        if (individualMetadataFieldsToExport != null && !individualMetadataFieldsToExport.isEmpty()) {
        	zos.putNextEntry(new ZipEntry(sModule + "__" + individualExportFiles.length + "individuals_metadata.tsv"));
        	zos.write("individual".getBytes());
	        IExportHandler.writeMetadataFile(sModule, sExportingUser, exportedIndividuals, individualMetadataFieldsToExport, zos);
	    	zos.closeEntry();
        }
        
        Collection<BasicDBList> variantDataQueries = varQueryWrapper.getVariantDataQueries();
        BasicDBObject varQuery = !variantDataQueries.isEmpty() ? new BasicDBObject("$and", variantDataQueries.iterator().next()) : new BasicDBObject();

        zos.putNextEntry(new ZipEntry(exportName + ".ped"));
        writeGenotypeFile(zos, sModule, exportedIndividuals, individualPopulations, nQueryChunkSize, markerSynonyms, individualExportFiles, warningFileWriter, progress);
    	zos.closeEntry();

        zos.putNextEntry(new ZipEntry(exportName + ".map"));
        String refPosPath = Assembly.getVariantRefPosPath(nAssemblyId);
        int nMarkerIndex = 0;
        ArrayList<Comparable> unassignedMarkers = new ArrayList<>();
    	String refPosPathWithTrailingDot = Assembly.getThreadBoundVariantRefPosPath() + ".";
    	Document projectionAndSortDoc = new Document(refPosPathWithTrailingDot + ReferencePosition.FIELDNAME_SEQUENCE, 1).append(refPosPathWithTrailingDot + ReferencePosition.FIELDNAME_START_SITE, 1);
		try (MongoCursor<Document> markerCursor = IExportHandler.getMarkerCursorWithCorrectCollation(mongoTemplate.getCollection(tmpVarCollName != null ? tmpVarCollName : mongoTemplate.getCollectionName(VariantData.class)), tmpVarCollName != null ? new Document() : new Document(varQuery), projectionAndSortDoc, nQueryChunkSize)) {
            progress.addStep("Writing map file");
            progress.moveToNextStep();
	        while (markerCursor.hasNext()) {
	            Document exportVariant = markerCursor.next();
	            Document refPos = (Document) Helper.readPossiblyNestedField(exportVariant, refPosPath, ";", null);
	            Long pos = refPos == null ? null : ((Number) refPos.get(ReferencePosition.FIELDNAME_START_SITE)).longValue();
	            String chrom = refPos == null ? null : (String) refPos.get(ReferencePosition.FIELDNAME_SEQUENCE);
                String markerId = ((String) exportVariant.get("_id"));
	            if (chrom == null)
	            	unassignedMarkers.add(markerId);
	            String exportedId = markerSynonyms == null ? markerId : markerSynonyms.get(markerId);
	            zos.write(((chrom == null ? "0" : chrom) + " " + exportedId + " " + 0 + " " + (pos == null ? 0 : pos) + LINE_SEPARATOR).getBytes());

                progress.setCurrentStepProgress(nMarkerIndex++ * 100 / markerCount);
	        }
		}
        zos.closeEntry();
        
        if (unassignedMarkers.size() > 0)
        	LOG.info("No chromosomal position found for " + unassignedMarkers.size() + " markers " + StringUtils.join(unassignedMarkers, ", "));

        if (warningFile.length() > 0) {
            progress.addStep("Adding lines to warning file");
            progress.moveToNextStep();
            progress.setPercentageEnabled(false);
            zos.putNextEntry(new ZipEntry(exportName + "-REMARKS.txt"));
            int nWarningCount = 0;
            BufferedReader in = new BufferedReader(new FileReader(warningFile));
            String sLine;
            while ((sLine = in.readLine()) != null) {
                zos.write((sLine + "\n").getBytes());
                progress.setCurrentStepProgress(nWarningCount++);
            }
            LOG.info("Number of Warnings for export (" + exportName + "): " + nWarningCount);
            in.close();
            zos.closeEntry();
        }
        warningFile.delete();

        zos.finish();
        zos.close();
        progress.setPercentageEnabled(true);
        progress.setCurrentStepProgress((short) 100);
    }
    
    public void writeGenotypeFile(OutputStream os, String sModule, Collection<String> individualsToExport, Map<String, String> individualPops, int nQueryChunkSize, Map<String, String> markerSynonyms, File[] individualExportFiles, FileWriter warningFileWriter, ProgressIndicator progress) throws IOException, InterruptedException {
        short nProgress = 0, nPreviousProgress = 0;
        
        int i = 0, nNConcurrentThreads = Math.max(1, Runtime.getRuntime().availableProcessors());	// use multiple threads so we can prepare several lines at once
        HashMap<Integer, StringBuilder> individualLines = new HashMap<>(nNConcurrentThreads);
        final ArrayList<Thread> threadsToWaitFor = new ArrayList<>(nNConcurrentThreads);
        final AtomicInteger initialStringBuilderCapacity = new AtomicInteger();        

        try
        {
            int nWrittenIndividualCount = 0;
	        for (final File f : individualExportFiles) {
	            if (progress.isAborted() || progress.getError() != null)
	            	return;

	        	final int nThreadIndex = i % nNConcurrentThreads;
	            Thread thread = new Thread() {
	                @Override
	                public void run() {
                        StringBuilder indLine = individualLines.get(nThreadIndex);
                        if (indLine == null) {
                            indLine = new StringBuilder((int) f.length() / 3 /* rough estimation */);
                            individualLines.put(nThreadIndex, indLine);
                        }

	                	BufferedReader in = null;
	    	            try {
		    	            in = new BufferedReader(new FileReader(f));
			                String individualId, line = in.readLine();	// read sample id
			                if (line != null) {
			                    individualId = line;
			                    String population = individualPops.get(line);
			                    String individualInfo = (population == null ? "." : population) + " " + individualId;
			                    indLine.append(individualInfo).append(" 0 0 0 ").append(getIndividualGenderCode(sModule, individualId));
			                } else {
			                    throw new Exception("Unable to read first line of temp export file " + f.getName());
			                }
			
			                int nMarkerIndex = 0;
			                while ((line = in.readLine()) != null) {
			                	String mostFrequentGenotype = findOutMostFrequentGenotype(line, warningFileWriter, nMarkerIndex, individualId);
			                    String[] alleles = mostFrequentGenotype == null ? new String[0] : mostFrequentGenotype.split(" ");
			                    if (alleles.length > 2 && warningFileWriter != null)
		                    		warningFileWriter.write("- More than 2 alleles found for variant n. " + nMarkerIndex + ", individual " + individualId + ". Exporting only the first 2 alleles.\n");
			
			                    String all1 = alleles.length == 0 ? "0" : alleles[0];
			                    String all2 = alleles.length == 0 ? "0" : alleles[alleles.length == 1 ? 0 : 1];
		                        indLine.append(" ").append(all1).append(" ").append(all2);
			
			                    nMarkerIndex++;
			                }
	    	            } catch (Exception e) {
	    	                LOG.error("Error exporting data", e);
	    	                progress.setError("Error exporting data: " + e.getClass().getSimpleName() + (e.getMessage() != null ? " - " + e.getMessage() : ""));
	    	            } finally {
	    	            	if (in != null)
								try {
									in.close();
								} catch (IOException e) {
									LOG.warn(e);
								}               
	    	            }

	    	            indLine.append(LINE_SEPARATOR);
	                }
	            };
	            
	            thread.start();
	            threadsToWaitFor.add(thread);

	            if (++i % nNConcurrentThreads == 0 || i == individualExportFiles.length) {
	                for (Thread t : threadsToWaitFor) // wait for all previously launched async threads
	               		t.join();
	                
                    for (int j=0; j<nNConcurrentThreads && nWrittenIndividualCount++ < individualExportFiles.length; j++) {
                        StringBuilder indLine = individualLines.get(j);
                        if (indLine == null || indLine.length() == 0)
                            LOG.warn("No line to export for individual " + j);
                        else {
                            os.write(indLine.toString().getBytes());
                            individualLines.put(j, new StringBuilder(initialStringBuilderCapacity.get()));
                        }
                    }

    	            nProgress = (short) (i * 100 / individualExportFiles.length);
    	            if (nProgress > nPreviousProgress) {
    	                progress.setCurrentStepProgress(nProgress);
    	                nPreviousProgress = nProgress;
    	            }
	                threadsToWaitFor.clear();
	            }
	        }
        }
        finally
        {
        	for (File f : individualExportFiles)
                if (!f.delete()) {
                    f.deleteOnExit();
                    LOG.info("Unable to delete tmp export file " + f.getAbsolutePath());
                }
        }
    	if (warningFileWriter != null)
    		warningFileWriter.close();
    }

    /* (non-Javadoc)
	 * @see fr.cirad.mgdb.exporting.IExportHandler#getStepList()
     */
    @Override
    public List<String> getStepList() {
        return Arrays.asList(new String[]{"Exporting data to PLINK format"});
    }

    /**
     * Gets the individual gender code.
     *
     * @param sModule the module
     * @param individual the individual
     * @return the individual gender code
     */
    protected String getIndividualGenderCode(String sModule, String individual) {
        return "U";
    }
    
	@Override
	public String[] getExportDataFileExtensions() {
		return new String[] {"ped", "map"};
	}

    @Override
    public int[] getSupportedPloidyLevels() {
        return new int[] {2};
    }
}