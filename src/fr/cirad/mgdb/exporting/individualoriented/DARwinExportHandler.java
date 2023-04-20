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
import java.util.Comparator;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Scanner;
import java.util.TreeMap;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.Collectors;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;

import org.apache.log4j.Logger;
import org.bson.Document;
import org.springframework.data.mongodb.core.MongoTemplate;
import org.springframework.data.mongodb.core.query.Criteria;
import org.springframework.data.mongodb.core.query.Query;

import com.mongodb.client.MongoCollection;
import com.mongodb.client.MongoCursor;

import fr.cirad.mgdb.exporting.IExportHandler;
import fr.cirad.mgdb.model.mongo.maintypes.Assembly;
import fr.cirad.mgdb.model.mongo.maintypes.GenotypingProject;
import fr.cirad.mgdb.model.mongo.maintypes.Individual;
import fr.cirad.mgdb.model.mongo.maintypes.VariantData;
import fr.cirad.mgdb.model.mongo.subtypes.ReferencePosition;
import fr.cirad.mgdb.model.mongodao.MgdbDao;
import fr.cirad.tools.Helper;
import fr.cirad.tools.ProgressIndicator;
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
	public void exportData(OutputStream outputStream, String sModule, Integer nAssemblyId, String sExportingUser, File[] individualExportFiles, boolean fDeleteSampleExportFilesOnExit, ProgressIndicator progress, String tmpVarCollName, Document varQuery, long markerCount, Map<String, String> markerSynonyms, Collection<String> individualMetadataFieldsToExport, Map<String, InputStream> readyToExportFiles) throws Exception {
        MongoTemplate mongoTemplate = MongoTemplateManager.get(sModule);
        GenotypingProject aProject = mongoTemplate.findOne(new Query(Criteria.where(GenotypingProject.FIELDNAME_PLOIDY_LEVEL).exists(true)), GenotypingProject.class);
        if (aProject == null) {
            LOG.warn("Unable to find a project containing ploidy level information! Assuming ploidy level is 2.");
        }

        int ploidy = aProject == null ? 2 : aProject.getPloidyLevel();

        File warningFile = File.createTempFile("export_warnings_", "");
        FileWriter warningFileWriter = new FileWriter(warningFile);

        ZipOutputStream os = IExportHandler.createArchiveOutputStream(outputStream, readyToExportFiles);
		Assembly assembly = mongoTemplate.findOne(new Query(Criteria.where("_id").is(nAssemblyId)), Assembly.class);
		String exportName = sModule + (assembly != null ? "__" + assembly.getName() : "") + "__" + markerCount + "variants__" + individualExportFiles.length + "individuals";
        
        String missingGenotype = "";
        for (int j = 0; j < ploidy; j++)
            missingGenotype += "\tN";
        String finalMissingGenotype = missingGenotype;

        os.putNextEntry(new ZipEntry(exportName + ".var"));
        os.write(("@DARwin 5.0 - ALLELIC - " + ploidy + LINE_SEPARATOR + individualExportFiles.length + "\t" + markerCount * ploidy + LINE_SEPARATOR + "N°").getBytes());

        short nProgress = 0, nPreviousProgress = 0;
        MongoCollection<Document> varColl = mongoTemplate.getCollection(tmpVarCollName != null ? tmpVarCollName : mongoTemplate.getCollectionName(VariantData.class));

    	String refPosPathWithTrailingDot = Assembly.getThreadBoundVariantRefPosPath() + ".";
    	Document projectionAndSortDoc = new Document(refPosPathWithTrailingDot + ReferencePosition.FIELDNAME_SEQUENCE, 1).append(refPosPathWithTrailingDot + ReferencePosition.FIELDNAME_START_SITE, 1);
        int nQueryChunkSize = IExportHandler.computeQueryChunkSize(mongoTemplate, markerCount);
		try (MongoCursor<Document> markerCursor = IExportHandler.getMarkerCursorWithCorrectCollation(varColl, varQuery, projectionAndSortDoc, nQueryChunkSize)) {
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
        for (File indFile : individualExportFiles)
            try (Scanner scanner = new Scanner(indFile)) {
                exportedIndividuals.add(scanner.nextLine());
            }
        LinkedHashMap<String, Individual> indMap = MgdbDao.getInstance().loadIndividualsWithAllMetadata(sModule, sExportingUser, null, exportedIndividuals);

        TreeMap<Integer, Comparable> problematicMarkerIndexToNameMap = new TreeMap<Integer, Comparable>();
        ArrayList<String> distinctAlleles = new ArrayList<String>();    // the index of each allele will be used as its code
        String[] donFileContents = new String[indMap.size()];

        int i = 0, nNConcurrentThreads = Math.max(1, Runtime.getRuntime().availableProcessors());    // use multiple threads so we can prepare several lines at once
        HashMap<Integer, StringBuilder> individualLines = new HashMap<>(nNConcurrentThreads);

        final ArrayList<Thread> threadsToWaitFor = new ArrayList<>(nNConcurrentThreads);
        final AtomicInteger initialStringBuilderCapacity = new AtomicInteger();
        
        try {
            int nWrittenIndividualCount = 0;
            for (final File f : individualExportFiles) {
                if (progress.isAborted() || progress.getError() != null)
                    return;

                final int nThreadIndex = i % nNConcurrentThreads;
                final int count = i + 1;
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
                            String individualId, line = in.readLine();    // read sample id

                            if (line != null)
                                individualId = line;
                            else
                                throw new Exception("Unable to read first line of temp export file " + f.getName());
                            
                            StringBuilder donFileLineSB = new StringBuilder();
                            donFileLineSB.append(count).append("\t").append(individualId);
                            for (String header : individualMetadataFieldsToExport)
                                donFileLineSB.append("\t").append(Helper.nullToEmptyString(indMap.get(individualId).getAdditionalInfo().get(header)));
                            donFileLineSB.append(LINE_SEPARATOR);
                            donFileContents[count - 1] = donFileLineSB.toString();

                            indLine.append(LINE_SEPARATOR).append(count);
                            int nMarkerIndex = 0;

                            while ((line = in.readLine()) != null) {
                                String mostFrequentGenotype = null;
                                if (!line.isEmpty()) {
                                    List<String> genotypes = Helper.split(line, "|");
                                    if (genotypes.size() == 1)
                                        mostFrequentGenotype = genotypes.get(0);
                                    else {
                                        HashMap<Object, Integer> genotypeCounts = new HashMap<Object, Integer>();   // will help us to keep track of missing genotypes
                                        int highestGenotypeCount = 0;
        
                                        for (String genotype : genotypes) {
                                            if (genotype == null) {
                                                continue;   /* skip missing genotypes */
                                            }
                    
                                            int gtCount = 1 + Helper.getCountForKey(genotypeCounts, genotype);
                                            if (gtCount > highestGenotypeCount) {
                                                highestGenotypeCount = gtCount;
                                                mostFrequentGenotype = genotype;
                                            }
                                            genotypeCounts.put(genotype, gtCount);
                                        }
                    
                                        if (genotypeCounts.size() > 1) {
                                            if (warningFileWriter != null){
                                                List<Integer> reverseSortedGtCounts = genotypeCounts.values().stream().sorted(Comparator.reverseOrder()).collect(Collectors.toList());
                                                if (reverseSortedGtCounts.get(0) == reverseSortedGtCounts.get(1))
                                                    mostFrequentGenotype = null;
                                                warningFileWriter.write("- Dissimilar genotypes found for variant n. " + nMarkerIndex + ", individual " + individualId + ". " + (mostFrequentGenotype == null ? "Exporting as missing data" : "Exporting most frequent: " + mostFrequentGenotype) + "\n");
                                            }
                                            problematicMarkerIndexToNameMap.put(nMarkerIndex, "");
                                        }
                                    }
                                }

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
                if (!f.delete())
                {
                    f.deleteOnExit();
                    LOG.info("Unable to delete tmp export file " + f.getAbsolutePath());
                }
        }
        os.closeEntry();

        os.putNextEntry(new ZipEntry(exportName + ".don"));
        String donFileHeader = "@DARwin 5.0 - DON -" + LINE_SEPARATOR + individualExportFiles.length + "\t" + (1 + individualMetadataFieldsToExport.size()) + LINE_SEPARATOR + "N°" + "\t" + "individual";
        for (String header : individualMetadataFieldsToExport)
            donFileHeader += "\t" + header;
        os.write((donFileHeader + LINE_SEPARATOR).getBytes());
        for (String donIndLine : donFileContents)
            os.write(donIndLine.getBytes());

//        // now read variant names for those that induced warnings
//        int nMarkerIndex = 0;
//        try (MongoCursor<Document> markerCursor = IExportHandler.getMarkerCursorWithCorrectCollation(varColl, varQuery, nQueryChunkSize)) {
//            while (markerCursor.hasNext()) {
//                Document exportVariant = markerCursor.next();
//                if (problematicMarkerIndexToNameMap.containsKey(nMarkerIndex)) {
//                    Comparable markerId = (Comparable) exportVariant.get("_id");
//    
//                    if (markerSynonyms != null) {
//                        Comparable syn = markerSynonyms.get(markerId);
//                        if (syn != null)
//                            markerId = syn;
//                    }
//                    for (int j = 0; j < ploidy; j++)
//                        zos.write(("\t" + markerId).getBytes());
//    
//                    problematicMarkerIndexToNameMap.put(nMarkerIndex, markerId);
//                }
//            }
//        }
        os.closeEntry();

        warningFileWriter.close();
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