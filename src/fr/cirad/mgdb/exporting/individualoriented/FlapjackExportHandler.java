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
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;
import java.util.Map;
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
import com.mongodb.client.MongoCollection;
import com.mongodb.client.MongoCursor;

import fr.cirad.mgdb.exporting.IExportHandler;
import fr.cirad.mgdb.exporting.tools.ExportManager.ExportOutputs;
import fr.cirad.mgdb.model.mongo.maintypes.Assembly;
import fr.cirad.mgdb.model.mongo.maintypes.Individual;
import fr.cirad.mgdb.model.mongo.maintypes.VariantData;
import fr.cirad.mgdb.model.mongo.maintypes.VariantRunData;
import fr.cirad.mgdb.model.mongo.subtypes.ReferencePosition;
import fr.cirad.mgdb.model.mongo.subtypes.VariantRunDataId;
import fr.cirad.tools.Helper;
import fr.cirad.tools.ProgressIndicator;
import fr.cirad.tools.mgdb.VariantQueryWrapper;
import fr.cirad.tools.mongo.MongoTemplateManager;

/**
 * The Class FlapjackExportHandler.
 */
public class FlapjackExportHandler extends AbstractIndividualOrientedExportHandler {

    /**
     * The Constant LOG.
     */
    private static final Logger LOG = Logger.getLogger(FlapjackExportHandler.class);

    /* (non-Javadoc)
     * @see fr.cirad.mgdb.exporting.IExportHandler#getExportFormatName()
     */
    @Override
    public String getExportFormatName() {
        return "FLAPJACK";
    }

    /* (non-Javadoc)
     * @see fr.cirad.mgdb.exporting.IExportHandler#getExportFormatDescription()
     */
    @Override
    public String getExportFormatDescription() {
        return "Exports zipped GENOTYPE and MAP files (plus PHENOTYPE file if metadata selected). See <a target='_blank' href='https://ics.hutton.ac.uk/wiki/index.php/Flapjack_Help_-_Projects_and_Data_Formats'>https://ics.hutton.ac.uk/wiki/index.php/Flapjack_Help_-_Projects_and_Data_Formats</a> for more details";
    }
    
    @Override
    public String getExportArchiveExtension() {
        return "fjzip";
    }
    
    @Override
    public String getMetadataContentsPrefix() {
        return "# fjFile = PHENOTYPE" + LINE_SEPARATOR;
    }

    @Override
    public String getMetadataFileExtension() {
        return "phenotype";
    }

    @Override
    public void exportData(OutputStream outputStream, String sModule, Integer nAssemblyId, ExportOutputs exportOutputs, boolean fDeleteSampleExportFilesOnExit, ProgressIndicator progress, String tmpVarCollName, VariantQueryWrapper varQueryWrapper, long markerCount, Map<String, String> markerSynonyms, Map<String, String> individualPopulations, Map<String, InputStream> readyToExportFiles) throws Exception {
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
	        MongoTemplate mongoTemplate = MongoTemplateManager.get(sModule);
	        
	        Assembly assembly = mongoTemplate.findOne(new Query(Criteria.where("_id").is(nAssemblyId)), Assembly.class);
	        String exportName = IExportHandler.buildExportName(sModule, assembly, markerCount, exportOutputs.getGenotypeFiles().length, exportOutputs.isWorkWithSamples());
	        int nQueryChunkSize = IExportHandler.computeQueryChunkSize(mongoTemplate, markerCount);
	
	        
	        Collection<BasicDBList> variantDataQueries = varQueryWrapper.getVariantDataQueries();
	        Document variantQueryForTargetCollection = variantDataQueries.isEmpty() ? new Document() : (tmpVarCollName == null ? new Document("$and", variantDataQueries.iterator().next()) : (varQueryWrapper.getBareQueries().iterator().hasNext() ? new Document("$and", varQueryWrapper.getBareQueries().iterator().next()) : new Document()));

	        zos.putNextEntry(new ZipEntry(exportName + ".genotype"));
	        MongoCollection<Document> varColl = mongoTemplate.getCollection(tmpVarCollName != null ? tmpVarCollName : mongoTemplate.getCollectionName(VariantData.class));
	        writeGenotypeFile(zos, nQueryChunkSize, varColl, variantQueryForTargetCollection, markerSynonyms, exportOutputs.getGenotypeFiles(), warningOS, progress);
	    	zos.closeEntry();
	
	
	        zos.putNextEntry(new ZipEntry(exportName + ".map"));
	        zos.write(("# fjFile = MAP" + LINE_SEPARATOR).getBytes());
	        int nMarkerIndex = 0;
			byte[] separatorBytes = "\t".getBytes();
			String refPosPath = Assembly.getVariantRefPosPath(nAssemblyId);
	        ArrayList<String> unassignedMarkers = new ArrayList<>();
	
	    	String refPosPathWithTrailingDot = Assembly.getThreadBoundVariantRefPosPath() + ".";
	    	Document projectionAndSortDoc = new Document(refPosPathWithTrailingDot + ReferencePosition.FIELDNAME_SEQUENCE, 1).append(refPosPathWithTrailingDot + ReferencePosition.FIELDNAME_START_SITE, 1);
	    	try (MongoCursor<Document> markerCursor = IExportHandler.getMarkerCursorWithCorrectCollation(mongoTemplate.getCollection(tmpVarCollName != null ? tmpVarCollName : mongoTemplate.getCollectionName(VariantData.class)), variantQueryForTargetCollection, projectionAndSortDoc, nQueryChunkSize)) {
	            progress.addStep("Writing map file");
	            progress.moveToNextStep();
	            while (markerCursor.hasNext()) {
	                Document exportVariant = markerCursor.next();
	                Document refPos = (Document) Helper.readPossiblyNestedField(exportVariant, refPosPath, ";", null);
	                String markerId = (String) exportVariant.get("_id");
	                String chrom = refPos == null ? null : (String) refPos.get(ReferencePosition.FIELDNAME_SEQUENCE);
	                Long pos = chrom == null ? null : ((Number) refPos.get(ReferencePosition.FIELDNAME_START_SITE)).longValue();
	                if (chrom == null) 
	                    unassignedMarkers.add(markerId);
	    
	                String exportedId = markerSynonyms == null ? markerId : markerSynonyms.get(markerId);
	                zos.write(exportedId.getBytes());
	                zos.write(separatorBytes);
	                zos.write((chrom == null ? "0" : chrom).getBytes());
	                zos.write(separatorBytes);
	                zos.write((pos == null ? "0" : pos.toString()).getBytes());
	                zos.write(LINE_SEPARATOR.getBytes());
	    
	                progress.setCurrentStepProgress(nMarkerIndex++ * 100 / markerCount);
	            }
	        }
	        zos.closeEntry();
	        
	        if (unassignedMarkers.size() > 0)
	            LOG.info("No chromosomal position found for " + unassignedMarkers.size() + " markers " + StringUtils.join(unassignedMarkers, ", "));
	
	        warningOS.close();
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
	        
	        IExportHandler.writeZipEntryFromChunkFiles(zos, exportOutputs.getAnnotationFiles(), exportName + ".ann", IExportHandler.VEP_LIKE_HEADER_LINE);

	        zos.finish();
	        zos.close();
	    }
        finally {
            warningFile.delete();
        }
        progress.setPercentageEnabled(true);
        progress.setCurrentStepProgress((short) 100);
    }

    public void writeGenotypeFile(OutputStream os, int nQueryChunkSize, MongoCollection<Document> varColl, Document varQuery, Map<String, String> markerSynonyms, File[] individualExportFiles, OutputStream warningOS, ProgressIndicator progress) throws IOException, InterruptedException {
   		os.write(("# fjFile = GENOTYPE" + LINE_SEPARATOR).getBytes());
        
   		boolean fWorkingOnRuns = varColl.getNamespace().getCollectionName().equals(MongoTemplateManager.getMongoCollectionName(VariantRunData.class));
    	String refPosPathWithTrailingDot = Assembly.getThreadBoundVariantRefPosPath() + ".";
    	Document projectionAndSortDoc = new Document(refPosPathWithTrailingDot + ReferencePosition.FIELDNAME_SEQUENCE, 1).append(refPosPathWithTrailingDot + ReferencePosition.FIELDNAME_START_SITE, 1);
   		try (MongoCursor<Document> markerCursor = IExportHandler.getMarkerCursorWithCorrectCollation(varColl, varQuery, projectionAndSortDoc, nQueryChunkSize)) {
	        while (markerCursor.hasNext()) {
	            Document exportVariant = markerCursor.next();
	            String markerId = (String) (fWorkingOnRuns ? Helper.readPossiblyNestedField(exportVariant, "_id." + VariantRunDataId.FIELDNAME_VARIANT_ID, ";", null) : exportVariant.get("_id"));
	            String exportedId = markerSynonyms == null ? markerId : markerSynonyms.get(markerId);
	            os.write(("\t" + exportedId).getBytes());
	        }
		}
		catch (IOException ioe) {
			LOG.error("Error occurred while writing Flapjack genotype file. Deleting " + individualExportFiles.length + " temporary individual-oriented files.", ioe);
			for (File f : individualExportFiles)
				f.delete();
		}

        os.write(LINE_SEPARATOR.getBytes());

        short nProgress = 0, nPreviousProgress = 0;
        int i = 0, nNConcurrentThreads = Math.max(1, Runtime.getRuntime().availableProcessors());   // use multiple threads so we can prepare several lines at once
        StringBuilder[] individualLines = new StringBuilder[nNConcurrentThreads];

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
                        StringBuilder indLine = individualLines[nThreadIndex];
                        if (indLine == null) {
                            indLine = new StringBuilder((int) f.length() / 3 /* rough estimation */);
                            individualLines[nThreadIndex] = indLine;
                        }

                        BufferedReader in = null;
                        try {
                            in = new BufferedReader(new FileReader(f));
                            String individualId, line = in.readLine();  // read sample id
                            if (line != null) {
                                individualId = line;
                                indLine.append(individualId);
                            } else {
                                throw new Exception("Unable to read first line of temp export file " + f.getName());
                            }

                            int nMarkerIndex = 0;
                            while ((line = in.readLine()) != null) {
                            	String mostFrequentGenotype = findOutMostFrequentGenotype(line, warningOS, nMarkerIndex, individualId);
                                List<String> alleles = mostFrequentGenotype == null ? new ArrayList<>() : Helper.split(mostFrequentGenotype, " ");
                                if (alleles.size() == 0 || (alleles.size() == 1 && alleles.get(0).length() == 0))
                                    indLine.append("\t-");
                                else
                                {
                                    boolean fHomozygous = alleles.stream().distinct().count() == 1;
                                    for (int i=0; i<(fHomozygous ? 1 : alleles.size()); i++)
                                        indLine.append(i == 0 ? "\t" : "/").append(alleles.get(i));
                                }

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
                        if (initialStringBuilderCapacity.get() == 0)
                            initialStringBuilderCapacity.set(indLine.length());
                    }
                };
                
                thread.start();
                threadsToWaitFor.add(thread);

                if (++i % nNConcurrentThreads == 0 || i == individualExportFiles.length) {
                    for (Thread t : threadsToWaitFor) // wait for all previously launched async threads
                        t.join();
                    
                    for (int j=0; j<nNConcurrentThreads && nWrittenIndividualCount++ < individualExportFiles.length; j++) {
                        StringBuilder indLine = individualLines[j];
                        if (indLine == null || indLine.length() == 0)
                            LOG.warn("No line to export for individual " + j);
                        else {
                            os.write(indLine.toString().getBytes());
                            individualLines[j] = new StringBuilder(initialStringBuilderCapacity.get());
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
    }

    /* (non-Javadoc)
     * @see fr.cirad.mgdb.exporting.IExportHandler#getStepList()
     */
    @Override
    public List<String> getStepList() {
        return Arrays.asList(new String[]{"Exporting data to Flapjack format"});
    }

    /**
     * Gets the individual population.
     *
     * @param sModule the module
     * @param individual the individual
     * @return the individual population
     */
    protected String getIndividualPopulation(final String sModule, final String individual) {
        MongoTemplate mongoTemplate = MongoTemplateManager.get(sModule);
        return mongoTemplate.findById(individual, Individual.class).getPopulation();
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
        return new String[] {"genotype", "map", "phenotype"};
    }

	@Override
	public String getExportContentType() {
		return "application/x-fjzip";
	}

    @Override
    public int[] getSupportedPloidyLevels() {
        return null;
    }
}