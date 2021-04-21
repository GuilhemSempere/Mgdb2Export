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

import fr.cirad.mgdb.exporting.IExportHandler;
import fr.cirad.mgdb.model.mongo.maintypes.Individual;
import fr.cirad.mgdb.model.mongo.maintypes.VariantData;
import fr.cirad.mgdb.model.mongo.maintypes.VariantDataV2;
import fr.cirad.mgdb.model.mongo.maintypes.VariantRunData;
import fr.cirad.mgdb.model.mongo.subtypes.ReferencePosition;
import fr.cirad.tools.Helper;
import fr.cirad.tools.ProgressIndicator;
import fr.cirad.tools.mongo.MongoTemplateManager;

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

import org.apache.commons.lang.StringUtils;
import org.apache.log4j.Logger;
import org.bson.Document;
import org.springframework.data.mongodb.core.MongoTemplate;

import com.mongodb.client.FindIterable;
import com.mongodb.client.MongoCollection;
import com.mongodb.client.MongoCursor;

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
        return "Exports zipped GENOTYPE and MAP files. See <a target='_blank' href='https://ics.hutton.ac.uk/wiki/index.php/Flapjack_Help_-_Projects_and_Data_Formats'>https://ics.hutton.ac.uk/wiki/index.php/Flapjack_Help_-_Projects_and_Data_Formats</a> for more details";
    }
    
	@Override
	public String getExportArchiveExtension() {
		return "fjzip";
	}

	/* (non-Javadoc)
	 * @see fr.cirad.mgdb.exporting.individualoriented.AbstractIndividualOrientedExportHandler#exportData(java.io.OutputStream, java.lang.String, java.util.Collection, boolean, fr.cirad.tools.ProgressIndicator, com.mongodb.DBCursor, java.util.Map, java.util.Map)
     */
    @Override
    public void exportData(OutputStream outputStream, String sModule, Integer nAssemblyId, Collection<File> individualExportFiles, boolean fDeleteSampleExportFilesOnExit, ProgressIndicator progress, MongoCollection<Document> varColl, Document varQuery, long variantCount, Map<String, String> markerSynonyms, Map<String, InputStream> readyToExportFiles) throws Exception {
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

        MongoTemplate mongoTemplate = MongoTemplateManager.get(sModule);
        Number avgObjSize = (Number) mongoTemplate.getDb().runCommand(new Document("collStats", mongoTemplate.getCollectionName(VariantRunData.class))).get("avgObjSize");
        int nQueryChunkSize = (int) (nMaxChunkSizeInMb * 1024 * 1024 / avgObjSize.doubleValue());
        
		long markerCount = varColl.countDocuments(varQuery);
        String exportName = sModule + "__" + markerCount + "variants__" + individualExportFiles.size() + "individuals";
        zos.putNextEntry(new ZipEntry(exportName + ".genotype"));
        zos.write(("# fjFile = GENOTYPE" + LINE_SEPARATOR).getBytes());
        
        try (MongoCursor<Document> markerCursor = IExportHandler.getMarkerCursorWithCorrectCollation(varColl, Document.class, varQuery, null, nAssemblyId, nQueryChunkSize)) {
	        while (markerCursor.hasNext()) {
	            Document exportVariant = markerCursor.next();
	            Comparable markerId = (Comparable) exportVariant.get("_id");
	            Comparable exportedId = markerSynonyms == null ? markerId : markerSynonyms.get(markerId);
	            zos.write(("\t" + exportedId).getBytes());
	        }
		}
        zos.write(LINE_SEPARATOR.getBytes());

        TreeMap<Integer, Comparable> problematicMarkerIndexToNameMap = new TreeMap<Integer, Comparable>();
        short nProgress = 0, nPreviousProgress = 0;
        int i = 0;
        try
        {
	        for (File f : individualExportFiles) {
	            BufferedReader in = new BufferedReader(new FileReader(f));
	            try {
	                String individualId, line = in.readLine();	// read sample id
	                if (line != null) {
	                    individualId = line;
	                    zos.write(individualId.getBytes());
	                } else {
	                    throw new Exception("Unable to read first line of temp export file " + f.getName());
	                }
	
	                int nMarkerIndex = 0;
	                while ((line = in.readLine()) != null) {
	                    List<String> genotypes = Helper.split(line, "|");
	                    HashMap<Object, Integer> genotypeCounts = new HashMap<Object, Integer>();	// will help us to keep track of missing genotypes
	                    int highestGenotypeCount = 0;
	                    String mostFrequentGenotype = null;
	                    for (String genotype : genotypes) {
	                        if (genotype == null) {
	                            continue;	/* skip missing genotypes */
	                        }
	
	                        int gtCount = 1 + Helper.getCountForKey(genotypeCounts, genotype);
	                        if (gtCount > highestGenotypeCount) {
	                            highestGenotypeCount = gtCount;
	                            mostFrequentGenotype = genotype;
	                        }
	                        genotypeCounts.put(genotype, gtCount);
	                    }
	
	                    if (genotypeCounts.size() > 1) {
	                        warningFileWriter.write("- Dissimilar genotypes found for variant " + nMarkerIndex + ", individual " + individualId + ". Exporting most frequent: " + mostFrequentGenotype + "\n");
	                        problematicMarkerIndexToNameMap.put(nMarkerIndex, "");
	                    }
	
	                    String[] alleles = mostFrequentGenotype == null ? new String[0] : mostFrequentGenotype.split(" ");
	                    if (alleles.length > 2) {
	                        warningFileWriter.write("- More than 2 alleles found for variant " + nMarkerIndex + ", individual " + individualId + ". Exporting only the first 2 alleles.\n");
	                        problematicMarkerIndexToNameMap.put(nMarkerIndex, "");
	                    }

	                    if (alleles.length == 0 || (alleles.length == 1 && alleles[0].length() == 0))
	                    	zos.write(("\t-").getBytes());
	                    else
	                    {
		                    String all1 = alleles[0];
		                    String all2 = alleles[alleles.length == 1 ? 0 : 1];
	                        zos.write(("\t" + all1 + (!all2.equals(all1) ? "/" + all2 : "")).getBytes());
	                    }

	                    nMarkerIndex++;
	                }
	            } catch (Exception e) {
	                LOG.error("Error exporting data", e);
	                progress.setError("Error exporting data: " + e.getClass().getSimpleName() + (e.getMessage() != null ? " - " + e.getMessage() : ""));
	                return;
	            } finally {
	                in.close();               
	            }
	
	            if (progress.isAborted()) {
	                return;
	            }
	
	            nProgress = (short) (++i * 100 / individualExportFiles.size());
	            if (nProgress > nPreviousProgress) {
	                progress.setCurrentStepProgress(nProgress);
	                nPreviousProgress = nProgress;
	            }
	            zos.write('\n');
	        }
        }
        finally
        {
        	for (File f : individualExportFiles)
                if (!f.delete()) {
                    f.deleteOnExit();
                    LOG.info("Unable to delete tmp export file " + f.getAbsolutePath());
                }
        	zos.closeEntry();
        }
        warningFileWriter.close();

        zos.putNextEntry(new ZipEntry(exportName + ".map"));
        zos.write(("# fjFile = MAP" + LINE_SEPARATOR).getBytes());

        int nMarkerIndex = 0;
        ArrayList<String> unassignedMarkers = new ArrayList<>();
        try (MongoCursor<Document> markerCursor = IExportHandler.getMarkerCursorWithCorrectCollation(varColl, Document.class, varQuery, null, nAssemblyId, nQueryChunkSize)) {
	        while (markerCursor.hasNext()) {
	            Document exportVariant = markerCursor.next();
	            Document refPos = (Document) Helper.readPossiblyNestedField(exportVariant, (nAssemblyId != null ? VariantData.FIELDNAME_REFERENCE_POSITION + "." + nAssemblyId : VariantDataV2.FIELDNAME_REFERENCE_POSITION), "; ", null);
	            String markerId = (String) exportVariant.get("_id");
	            String chrom = refPos == null ? null : (String) refPos.get(ReferencePosition.FIELDNAME_SEQUENCE);
	            Long pos = refPos == null ? null : ((Number) refPos.get(ReferencePosition.FIELDNAME_START_SITE)).longValue();
	            if (chrom == null) 
	            	unassignedMarkers.add(markerId);
	
	            String exportedId = markerSynonyms == null ? markerId : markerSynonyms.get(markerId);
	            zos.write((exportedId + "\t" + (chrom == null ? "0" : chrom) + "\t" + (pos == null ? 0 : pos) + LINE_SEPARATOR).getBytes());
	
	            if (problematicMarkerIndexToNameMap.containsKey(nMarkerIndex)) {	// we are going to need this marker's name for the warning file
	            	String variantName = markerId;
	                if (markerSynonyms != null) {
	                	String syn = markerSynonyms.get(markerId);
	                    if (syn != null)
	                        variantName = syn;
	                }
	                problematicMarkerIndexToNameMap.put(nMarkerIndex, variantName);
	            }
	            nMarkerIndex++;
	        }
	    }
        zos.closeEntry();
        
        if (unassignedMarkers.size() > 0)
        	LOG.info("No chromosomal position found for " + unassignedMarkers.size() + " markers " + StringUtils.join(unassignedMarkers, ", "));

        if (warningFile.length() > 0) {
            zos.putNextEntry(new ZipEntry(exportName + "-REMARKS.txt"));
            int nWarningCount = 0;
            BufferedReader in = new BufferedReader(new FileReader(warningFile));
            String sLine;
            while ((sLine = in.readLine()) != null) {
                for (Integer aMarkerIndex : problematicMarkerIndexToNameMap.keySet())
                    sLine = sLine.replaceAll("__" + aMarkerIndex + "__", problematicMarkerIndexToNameMap.get(aMarkerIndex).toString());
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
		return new String[] {"genotype", "map"};
	}

	@Override
	public String getExportContentType() {
		return "application/x-fjzip";
	}
}
