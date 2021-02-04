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
import java.io.InputStream;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;

import org.apache.commons.collections.CollectionUtils;
import org.apache.log4j.Logger;
import org.bson.Document;
import org.springframework.data.mongodb.core.MongoTemplate;
import org.springframework.data.mongodb.core.query.Criteria;
import org.springframework.data.mongodb.core.query.Query;

import com.mongodb.client.MongoCollection;
import com.mongodb.client.MongoCursor;

import fr.cirad.mgdb.model.mongo.maintypes.GenotypingProject;
import fr.cirad.mgdb.model.mongo.maintypes.Individual;
import fr.cirad.mgdb.model.mongo.maintypes.VariantRunData;
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

    /* (non-Javadoc)
	 * @see fr.cirad.mgdb.exporting.individualoriented.AbstractIndividualOrientedExportHandler#exportData(java.io.OutputStream, java.lang.String, java.util.Collection, boolean, fr.cirad.tools.ProgressIndicator, com.mongodb.DBCursor, java.util.Map, java.util.Map)
     */
    @Override
    public void exportData(OutputStream outputStream, String sModule, Integer nAssemblyId,Collection<File> individualExportFiles, boolean fDeleteSampleExportFilesOnExit, ProgressIndicator progress, MongoCollection<Document> varColl, Document varQuery, Map<String, String> markerSynonyms, Map<String, InputStream> readyToExportFiles) throws Exception {
        MongoTemplate mongoTemplate = MongoTemplateManager.get(sModule);
        GenotypingProject aProject = mongoTemplate.findOne(new Query(Criteria.where(GenotypingProject.FIELDNAME_PLOIDY_LEVEL).exists(true)), GenotypingProject.class);
        if (aProject == null) {
            LOG.warn("Unable to find a project containing ploidy level information! Assuming ploidy level is 2.");
        }

        int ploidy = aProject == null ? 2 : aProject.getPloidyLevel();

        File warningFile = File.createTempFile("export_warnings_", "");
        FileWriter warningFileWriter = new FileWriter(warningFile);

        ZipOutputStream zos = new ZipOutputStream(outputStream);

        if (readyToExportFiles != null) {
            for (String readyToExportFile : readyToExportFiles.keySet()) {
                zos.putNextEntry(new ZipEntry(readyToExportFile));
                InputStream inputStream = readyToExportFiles.get(readyToExportFile);
                byte[] dataBlock = new byte[1024];
                int count = inputStream.read(dataBlock, 0, 1024);
                while (count != -1) {
                    zos.write(dataBlock, 0, count);
                    count = inputStream.read(dataBlock, 0, 1024);
                }
                zos.closeEntry();
            }
        }

		long markerCount = varColl.countDocuments(varQuery);
        String exportName = sModule + "__" + markerCount + "variants__" + individualExportFiles.size() + "individuals";

        StringBuffer donFileContents = new StringBuffer();

        int count = 0;
        String missingGenotype = "";
        for (int j = 0; j < ploidy; j++) {
            missingGenotype += "\tN";
        }

        zos.putNextEntry(new ZipEntry(exportName + ".var"));
        zos.write(("@DARwin 5.0 - ALLELIC - " + ploidy + LINE_SEPARATOR + individualExportFiles.size() + "\t" + markerCount * ploidy + LINE_SEPARATOR + "N°").getBytes());

        short nProgress = 0, nPreviousProgress = 0;
        Number avgObjSize = (Number) mongoTemplate.getDb().runCommand(new Document("collStats", mongoTemplate.getCollectionName(VariantRunData.class))).get("avgObjSize");
    	int nQueryChunkSize = (int) (nMaxChunkSizeInMb * 1024 * 1024 / avgObjSize.doubleValue());

        int nMarkerIndex = 0;
        try (MongoCursor<Document> markerCursor = varColl.find(varQuery).projection(projectionDoc(nAssemblyId)).sort(sortDoc(nAssemblyId)).noCursorTimeout(true).collation(collationObj).batchSize(nQueryChunkSize).iterator()) {
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
	                zos.write(("\t" + markerId).getBytes());
	            }
			}
        }

        List<String> indInfoHeaders = new ArrayList<>();
        TreeMap<Integer, Comparable> problematicMarkerIndexToNameMap = new TreeMap<Integer, Comparable>();
        ArrayList<String> distinctAlleles = new ArrayList<String>();	// the index of each allele will be used as its code
        int i = 0;
        try {
            for (File f : individualExportFiles) {
                BufferedReader in = new BufferedReader(new FileReader(f));
                try {
                    String individualId, line = in.readLine();	// read sample id

                    if (line != null)
                        individualId = line;
                    else
                        throw new Exception("Unable to read first line of temp export file " + f.getName());
                    
                    HashMap<String, Comparable> indInfo = mongoTemplate.findById(individualId, Individual.class).getAdditionalInfo();
                    if (indInfoHeaders.size() == 0)
                    	indInfoHeaders.addAll(indInfo.keySet());
                    else
                    {	// deal with fields we might have not yet encountered
	                    Collection<String> newHeaders = CollectionUtils.subtract(indInfo.keySet(), indInfoHeaders);
	                    if (newHeaders.size() > 0)
	                    	indInfoHeaders.addAll(newHeaders);
                    }
                    	
                    donFileContents.append(++count + "\t" + individualId);
                    for (String header : indInfoHeaders)
                    	donFileContents.append("\t" + Helper.nullToEmptyString(indInfo.get(header)));
                    donFileContents.append(LINE_SEPARATOR);

                    zos.write((LINE_SEPARATOR + count).getBytes());
                    nMarkerIndex = 0;

                    while ((line = in.readLine()) != null) {
                        List<String> genotypes = Helper.split(line, "|");
                        HashMap<Object, Integer> genotypeCounts = new HashMap<Object, Integer>();	// will help us to keep track of missing genotypes
                        int highestGenotypeCount = 0;
                        String mostFrequentGenotype = null;
                        for (String genotype : genotypes) {
                            if (genotype == null)
                                continue;	/* skip missing genotypes */

                            int gtCount = 1 + Helper.getCountForKey(genotypeCounts, genotype);
                            if (gtCount > highestGenotypeCount) {
                                highestGenotypeCount = gtCount;
                                mostFrequentGenotype = genotype;
                            }
                            genotypeCounts.put(genotype, gtCount);
                        }

                        if (genotypeCounts.size() > 1) {
                            warningFileWriter.write("- Dissimilar genotypes found for variant __" + nMarkerIndex + "__, individual " + individualId + ". Exporting most frequent: " + mostFrequentGenotype + "\n");
                            problematicMarkerIndexToNameMap.put(nMarkerIndex, "");
                        }

                        String codedGenotype = "";
                        if (mostFrequentGenotype != null && mostFrequentGenotype.length() > 0) {
                            for (String allele : mostFrequentGenotype.split(" ")) {
                                if (!distinctAlleles.contains(allele)) {
                                    distinctAlleles.add(allele);
                                }
                                codedGenotype += "\t" + distinctAlleles.indexOf(allele);
                            }
                        } else {
                            codedGenotype = missingGenotype.replaceAll("N", "-1");	// missing data is coded as -1		
                        }
                        zos.write(codedGenotype.getBytes());

                        nMarkerIndex++;
                    }
                } catch (Exception e) {
                    LOG.error("Error exporting data", e);
                    progress.setError("Error exporting data: " + e.getClass().getSimpleName() + (e.getMessage() != null ? " - " + e.getMessage() : ""));
                    return;
                } finally {
                    in.close();
                    warningFileWriter.close();
                }

                if (progress.isAborted()) {
                    return;
                }

                nProgress = (short) (++i * 100 / individualExportFiles.size());
                if (nProgress > nPreviousProgress) {
                    //				LOG.debug("============= doDARwinExport (" + i + "): " + nProgress + "% =============");
                    progress.setCurrentStepProgress(nProgress);
                    nPreviousProgress = nProgress;
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
        zos.closeEntry();

        zos.putNextEntry(new ZipEntry(exportName + ".don"));
        
        String donFileHeader = "@DARwin 5.0 - DON -" + LINE_SEPARATOR + individualExportFiles.size() + "\t" + (1 + indInfoHeaders.size()) + LINE_SEPARATOR + "N°" + "\t" + "individual";
        for (String header : indInfoHeaders)
        	donFileHeader += "\t" + header;
        zos.write((donFileHeader + LINE_SEPARATOR).getBytes());
        zos.write(donFileContents.toString().getBytes());

        // now read variant names for those that induced warnings
        nMarkerIndex = 0;
        try (MongoCursor<Document> markerCursor = varColl.find(varQuery).projection(projectionDoc(nAssemblyId)).sort(sortDoc(nAssemblyId)).noCursorTimeout(true).collation(collationObj).batchSize(nQueryChunkSize).iterator()) {
	        while (markerCursor.hasNext()) {
	            Document exportVariant = markerCursor.next();
	            if (problematicMarkerIndexToNameMap.containsKey(nMarkerIndex)) {
	                Comparable markerId = (Comparable) exportVariant.get("_id");
	
	                if (markerSynonyms != null) {
	                    Comparable syn = markerSynonyms.get(markerId);
	                    if (syn != null)
	                        markerId = syn;
	                }
	                for (int j = 0; j < ploidy; j++)
	                    zos.write(("\t" + markerId).getBytes());
	
	                problematicMarkerIndexToNameMap.put(nMarkerIndex, markerId);
	            }
	        }
		}

        warningFileWriter.close();
        if (warningFile.length() > 0) {
            zos.putNextEntry(new ZipEntry(exportName + "-REMARKS.txt"));
            int nWarningCount = 0;
            BufferedReader in = new BufferedReader(new FileReader(warningFile));
            String sLine;
            while ((sLine = in.readLine()) != null) {
                for (Integer aMarkerIndex : problematicMarkerIndexToNameMap.keySet()) {
                    sLine = sLine.replaceAll("__" + aMarkerIndex + "__", problematicMarkerIndexToNameMap.get(aMarkerIndex).toString());
                }
                zos.write((sLine + "\n").getBytes());
                nWarningCount++;
            }
            LOG.info("Number of Warnings for export (" + exportName + "): " + nWarningCount);
            in.close();
            zos.closeEntry();
        }
        warningFile.delete();

        zos.finish();
        zos.close();
        progress.setCurrentStepProgress((short) 100);
    }

	@Override
	public String[] getExportDataFileExtensions() {
		return new String[] {"don", "var"};
	}
}
