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
import java.util.Scanner;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;

import org.apache.log4j.Logger;
import org.bson.Document;
import org.springframework.data.mongodb.core.MongoTemplate;
import org.springframework.data.mongodb.core.query.Criteria;
import org.springframework.data.mongodb.core.query.Query;

import com.mongodb.BasicDBList;
import com.mongodb.client.MongoCursor;

import fr.cirad.mgdb.exporting.IExportHandler;
import fr.cirad.mgdb.exporting.tools.ExportManager.ExportOutputs;
import fr.cirad.mgdb.model.mongo.maintypes.Assembly;
import fr.cirad.mgdb.model.mongo.maintypes.GenotypingProject;
import fr.cirad.mgdb.model.mongo.maintypes.VariantData;
import fr.cirad.mgdb.model.mongo.subtypes.ReferencePosition;
import fr.cirad.tools.ProgressIndicator;
import fr.cirad.tools.mgdb.VariantQueryWrapper;
import fr.cirad.tools.mongo.MongoTemplateManager;

/**
 * The Class DARwinExportHandler.
 */
public class DARwinExportHandler extends AbstractIndividualOrientedExportHandler {

    /**
     * The Constant LOG.
     */
    private static final Logger LOG = Logger.getLogger(DARwinExportHandler.class);

    /* (non-Javadoc)
     * @see fr.cirad.mgdb.exporting.IExportHandler#getExportFormatName()
     */
    @Override
    public String getExportFormatName() {
        return "DARwin";
    }

    /* (non-Javadoc)
     * @see fr.cirad.mgdb.exporting.IExportHandler#getExportFormatDescription()
     */
    @Override
    public String getExportFormatDescription() {
        return "Exports data in DARwin Format. See <a target='_blank' href='http://darwin.cirad.fr/'>http://darwin.cirad.fr/</a> for more details";
    }

    /* (non-Javadoc)
     * @see fr.cirad.mgdb.exporting.IExportHandler#getStepList()
     */
    @Override
    public List<String> getStepList() {
        return Arrays.asList(new String[]{"Exporting data to DARWIN format"});
    }
    
    @Override
    public String getExportArchiveExtension() {
        return "zip";
    }

	@Override
	public void exportData(OutputStream outputStream, String sModule, Integer nAssemblyId, ExportOutputs exportOutputs, boolean fDeleteSampleExportFilesOnExit, ProgressIndicator progress, String tmpVarCollName, VariantQueryWrapper varQueryWrapper, long markerCount, Map<String, String> markerSynonyms, Map<String, String> individualPopulations, Map<String, InputStream> readyToExportFiles) throws Exception {
        MongoTemplate mongoTemplate = MongoTemplateManager.get(sModule);
        GenotypingProject aProject = mongoTemplate.findOne(new Query(Criteria.where(GenotypingProject.FIELDNAME_PLOIDY_LEVEL).exists(true)), GenotypingProject.class);
        if (aProject == null) {
            LOG.warn("Unable to find a project containing ploidy level information! Assuming ploidy level is 2.");
        }

        int ploidy = aProject == null ? 2 : aProject.getPloidyLevel();
        
        // save existing warnings into a temp file so we can append to it
        File warningFile = File.createTempFile("export_warnings_", "");
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

        ZipOutputStream os = IExportHandler.createArchiveOutputStream(outputStream, readyToExportFiles, null);
		Assembly assembly = mongoTemplate.findOne(new Query(Criteria.where("_id").is(nAssemblyId)), Assembly.class);
		String exportName = IExportHandler.buildExportName(sModule, assembly, markerCount, exportOutputs.getGenotypeFiles().length, exportOutputs.isWorkWithSamples());
        
        String missingGenotype = "";
        for (int j = 0; j < ploidy; j++)
            missingGenotype += "\tN";
        String finalMissingGenotype = missingGenotype;

        os.putNextEntry(new ZipEntry(exportName + ".var"));
        os.write(("@DARwin 5.0 - ALLELIC - " + ploidy + LINE_SEPARATOR + exportOutputs.getGenotypeFiles().length + "\t" + markerCount * ploidy + LINE_SEPARATOR + "N°").getBytes());

        short nProgress = 0, nPreviousProgress = 0;

        Collection<BasicDBList> variantDataQueries = varQueryWrapper.getVariantDataQueries();
        Document variantQueryForTargetCollection = variantDataQueries.isEmpty() ? new Document() : (tmpVarCollName == null ? new Document("$and", variantDataQueries.iterator().next()) : (varQueryWrapper.getBareQueries().iterator().hasNext() ? new Document("$and", varQueryWrapper.getBareQueries().iterator().next()) : new Document()));

        String refPosPathWithTrailingDot = Assembly.getThreadBoundVariantRefPosPath() + ".";
    	Document projectionAndSortDoc = new Document(refPosPathWithTrailingDot + ReferencePosition.FIELDNAME_SEQUENCE, 1).append(refPosPathWithTrailingDot + ReferencePosition.FIELDNAME_START_SITE, 1);
        int nQueryChunkSize = IExportHandler.computeQueryChunkSize(mongoTemplate, markerCount);
        try (MongoCursor<Document> markerCursor = IExportHandler.getMarkerCursorWithCorrectCollation(mongoTemplate.getCollection(tmpVarCollName != null ? tmpVarCollName : mongoTemplate.getCollectionName(VariantData.class)), variantQueryForTargetCollection, projectionAndSortDoc, nQueryChunkSize)) {
            while (markerCursor.hasNext()) {
                Document exportVariant = markerCursor.next();
                String markerId = (String) exportVariant.get("_id");
    
                if (markerSynonyms != null) {
                    String syn = markerSynonyms.get(markerId);
                    if (syn != null) {
                        markerId = syn;
                    }
                }
                for (int j = 0; j < ploidy; j++) {
                    os.write(("\t" + markerId).getBytes());
                }
            }
        }
        
        ArrayList<String> exportedIndividuals = new ArrayList<>();
        for (File indFile : exportOutputs.getGenotypeFiles())
            try (Scanner scanner = new Scanner(indFile)) {
                exportedIndividuals.add(scanner.nextLine());
            }
//        LinkedHashMap<String, Individual> indMap = MgdbDao.getInstance().loadIndividualsWithAllMetadata(sModule, sExportingUser, null, exportedIndividuals, null);
        ArrayList<String> distinctAlleles = new ArrayList<String>();    // the index of each allele will be used as its code
//        String[] donFileContents = new String[indMap.size()];
        
//        final Collection<String> finalMdFieldsToExport = individualMetadataFieldsToExport != null ? individualMetadataFieldsToExport : indMap.values().stream().flatMap(individual -> individual.getAdditionalInfo().keySet().stream()).collect(Collectors.toSet());

        int i = 0, nNConcurrentThreads = Math.max(1, Runtime.getRuntime().availableProcessors());    // use multiple threads so we can prepare several lines at once
        StringBuilder[] individualLines = new StringBuilder[nNConcurrentThreads];

        final ArrayList<Thread> threadsToWaitFor = new ArrayList<>(nNConcurrentThreads);
        final AtomicInteger initialStringBuilderCapacity = new AtomicInteger();
        
        try {
            int nWrittenIndividualCount = 0;
            for (final File f : exportOutputs.getGenotypeFiles()) {
                if (progress.isAborted() || progress.getError() != null)
                    return;

                final int nThreadIndex = i % nNConcurrentThreads;
                final int count = i + 1;
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
                            String individualId, line = in.readLine();    // read sample id

                            if (line != null)
                                individualId = line;
                            else
                                throw new Exception("Unable to read first line of temp export file " + f.getName());
                            
//                            StringBuilder donFileLineSB = new StringBuilder();
//                            donFileLineSB.append(count).append("\t").append(individualId);
//                            for (String header : finalMdFieldsToExport)
//                                donFileLineSB.append("\t").append(Helper.nullToEmptyString(indMap.get(individualId).getAdditionalInfo().get(header)));
//                            donFileLineSB.append(LINE_SEPARATOR);
//                            donFileContents[count - 1] = donFileLineSB.toString();

                            indLine.append(LINE_SEPARATOR).append(count);
                            int nMarkerIndex = 0;

                            while ((line = in.readLine()) != null) {
                            	String mostFrequentGenotype = findOutMostFrequentGenotype(line, warningOS, nMarkerIndex, individualId);
                                String codedGenotype = "";
                                if (mostFrequentGenotype != null && mostFrequentGenotype.length() > 0) {
                                    for (String allele : mostFrequentGenotype.split(" ")) {
                                        if (!distinctAlleles.contains(allele)) {
                                            distinctAlleles.add(allele);
                                        }
                                        codedGenotype += "\t" + distinctAlleles.indexOf(allele);
                                    }
                                } else
                                    codedGenotype = finalMissingGenotype.replaceAll("N", "-1");    // missing data is coded as -1        
                                indLine.append(codedGenotype);

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
                    }
                };
                
                thread.start();
                threadsToWaitFor.add(thread);

                if (++i % nNConcurrentThreads == 0 || i == exportOutputs.getGenotypeFiles().length) {
                    for (Thread t : threadsToWaitFor) // wait for all previously launched async threads
                           t.join();
                    
                    for (int j=0; j<nNConcurrentThreads && nWrittenIndividualCount++ < exportOutputs.getGenotypeFiles().length; j++) {
                        StringBuilder indLine = individualLines[j];
                        if (indLine == null || indLine.length() == 0)
                            LOG.warn("No line to export for individual " + j);
                        else {
                            os.write(indLine.toString().getBytes());
                            individualLines[j] = new StringBuilder(initialStringBuilderCapacity.get());
                        }
                    }

                    nProgress = (short) (i * 100 / exportOutputs.getGenotypeFiles().length);
                    if (nProgress > nPreviousProgress) {
                        progress.setCurrentStepProgress(nProgress);
                        nPreviousProgress = nProgress;
                    }
                    threadsToWaitFor.clear();
                }
            }
        }
        finally {
            for (File f : exportOutputs.getGenotypeFiles())
                if (!f.delete()) {
                    f.deleteOnExit();
                    LOG.info("Unable to delete tmp export file " + f.getAbsolutePath());
                }
            warningFile.delete();
        }
        os.closeEntry();

        os.putNextEntry(new ZipEntry(exportName + ".don"));
        String mdFileContents = exportOutputs == null || exportOutputs.getMetadataFileContents() == null ? "" : exportOutputs.getMetadataFileContents();
        if (mdFileContents.isEmpty())
        	mdFileContents = "individual" + LINE_SEPARATOR + String.join(LINE_SEPARATOR, exportedIndividuals);
    	String[] mdLines = mdFileContents.split(LINE_SEPARATOR);
    	for (int j=0; j<mdLines.length; j++) {
    		if (j == 0)
    			os.write(("@DARwin 5.0 - DON -" + LINE_SEPARATOR + exportOutputs.getGenotypeFiles().length + "\t" + (mdLines[0].split("\t", -1).length) + LINE_SEPARATOR + "N°" + "\t").getBytes());
    		else
    			os.write((j + "\t").getBytes());
        	os.write((mdLines[j] + LINE_SEPARATOR).getBytes());
    	}
        os.closeEntry();

        warningOS.close();
        if (warningFile.length() > 0) {
            progress.addStep("Adding lines to warning file");
            progress.moveToNextStep();
            progress.setPercentageEnabled(false);
            os.putNextEntry(new ZipEntry(exportName + "-REMARKS.txt"));
            int nWarningCount = 0;
            BufferedReader in = new BufferedReader(new FileReader(warningFile));
            String sLine;
            while ((sLine = in.readLine()) != null) {
                os.write((sLine + "\n").getBytes());
                progress.setCurrentStepProgress(nWarningCount++);
            }
            LOG.info("Number of Warnings for export (" + exportName + "): " + nWarningCount);
            in.close();
            os.closeEntry();
        }
        warningFile.delete();

        os.finish();
        os.close();
        progress.setPercentageEnabled(true);
        progress.setCurrentStepProgress((short) 100);
    }

    @Override
    public String[] getExportDataFileExtensions() {
        return new String[] {"don", "var"};
    }
    
    @Override
    public int[] getSupportedPloidyLevels() {
        return null;
    }
}