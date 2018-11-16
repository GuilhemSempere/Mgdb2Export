/*******************************************************************************
 * MGDB Export - Mongo Genotype DataBase, export handlers
 * Copyright (C) 2016, 2018, <CIRAD> <IRD>
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
package fr.cirad.mgdb.exporting.markeroriented;

import fr.cirad.mgdb.exporting.tools.AsyncExportTool;
import fr.cirad.mgdb.exporting.tools.AsyncExportTool.AbstractDataOutputHandler;
import fr.cirad.mgdb.model.mongo.maintypes.GenotypingSample;
import fr.cirad.mgdb.model.mongo.maintypes.Individual;
import fr.cirad.mgdb.model.mongo.maintypes.VariantData;
import fr.cirad.mgdb.model.mongo.maintypes.VariantRunData;
import fr.cirad.mgdb.model.mongo.subtypes.SampleGenotype;
import fr.cirad.mgdb.model.mongodao.MgdbDao;
import fr.cirad.tools.AlphaNumericComparator;
import fr.cirad.tools.Helper;
import fr.cirad.tools.ProgressIndicator;
import fr.cirad.tools.mongo.MongoTemplateManager;
import htsjdk.variant.variantcontext.VariantContext.Type;

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
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.stream.Collectors;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;

import org.apache.commons.lang.StringUtils;
import org.apache.log4j.Logger;
import org.springframework.data.mongodb.core.MongoTemplate;

import com.mongodb.DBCursor;
import com.mongodb.DBObject;

/**
 * The Class EigenstratExportHandler.
 */
public class EigenstratExportHandler extends AbstractMarkerOrientedExportHandler {

    /**
     * The Constant LOG.
     */
    private static final Logger LOG = Logger.getLogger(EigenstratExportHandler.class);

    /**
     * The Constant EIGENSTRAT_FORMAT.
     */
    public static final String EIGENSTRAT_FORMAT = "EIGENSTRAT";

    /**
     * The supported variant types.
     */
    private static List<String> supportedVariantTypes;

    static {
        supportedVariantTypes = new ArrayList<String>();
        supportedVariantTypes.add(Type.SNP.toString());
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

    /**
     * Gets the populations from samples.
     *
     * @param sModule the module
     * @param samples the samples
     * @return the populations from samples
     */
    @SuppressWarnings("unchecked")
    protected List<String> getPopulationsFromSamples(final String sModule, final List<GenotypingSample> samples) {
        ArrayList<String> result = new ArrayList<String>();
        for (Individual individual : MgdbDao.getIndividualsFromSamples(sModule, samples)) {
            result.add(individual.getPopulation());
        }
        return result;
    }

    /* (non-Javadoc)
	 * @see fr.cirad.mgdb.exporting.IExportHandler#getExportFormatName()
     */
    @Override
    public String getExportFormatName() {
        return EIGENSTRAT_FORMAT;
    }

    /* (non-Javadoc)
	 * @see fr.cirad.mgdb.exporting.IExportHandler#getExportFormatDescription()
     */
    @Override
    public String getExportFormatDescription() {
        return "Exports zipped ind, snp and eigenstratgeno files, along with an optional remark-file. See <a target='_blank' href='https://github.com/argriffing/eigensoft/blob/master/CONVERTF/README'>https://github.com/argriffing/eigensoft/blob/master/CONVERTF/README</a> for more details";
    }

    /* (non-Javadoc)
	 * @see fr.cirad.mgdb.exporting.markeroriented.AbstractMarkerOrientedExportHandler#getSupportedVariantTypes()
     */
    @Override
    public List<String> getSupportedVariantTypes() {
        return supportedVariantTypes;
    }

	@Override
	public String getExportArchiveExtension() {
		return "zip";
	}

    /* (non-Javadoc)
	 * @see fr.cirad.mgdb.exporting.markeroriented.AbstractMarkerOrientedExportHandler#exportData(java.io.OutputStream, java.lang.String, java.util.List, fr.cirad.tools.ProgressIndicator, com.mongodb.DBCursor, java.util.Map, int, int, java.util.Map)
     */
    @Override
    public void exportData(OutputStream outputStream, String sModule, Collection<GenotypingSample> samples1, Collection<GenotypingSample> samples2, ProgressIndicator progress, DBCursor markerCursor, Map<String, String> markerSynonyms, HashMap<String, Integer> annotationFieldThresholds, HashMap<String, Integer> annotationFieldThresholds2, List<GenotypingSample> samplesToExport, Map<String, InputStream> readyToExportFiles) throws Exception {
        File warningFile = File.createTempFile("export_warnings_", "");
        FileWriter warningFileWriter = new FileWriter(warningFile);
        File snpFile = null;

        try {
            snpFile = File.createTempFile("snpFile", "");
            FileWriter snpFileWriter = new FileWriter(snpFile);

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
            int markerCount = markerCursor.count();

    		List<String> individuals1 = MgdbDao.getIndividualsFromSamples(sModule, samples1).stream().map(ind -> ind.getId()).collect(Collectors.toList());	
    		List<String> individuals2 = MgdbDao.getIndividualsFromSamples(sModule, samples2).stream().map(ind -> ind.getId()).collect(Collectors.toList());

    		List<String> sortedIndividuals = samplesToExport.stream().map(gs -> gs.getIndividual()).distinct().sorted(new AlphaNumericComparator<String>()).collect(Collectors.toList());
 
            String exportName = sModule + "__" + markerCount + "variants__" + sortedIndividuals.size() + "individuals";
            
            zos.putNextEntry(new ZipEntry(exportName + ".eigenstratgeno"));

			final Map<Integer, String> sampleIdToIndividualMap = new HashMap<>();
			for (GenotypingSample gs : samplesToExport)
				sampleIdToIndividualMap.put(gs.getId(), gs.getIndividual());

            ArrayList<Comparable> unassignedMarkers = new ArrayList<>();
            
    		AbstractDataOutputHandler<Integer, LinkedHashMap<VariantData, Collection<VariantRunData>>> dataOutputHandler = new AbstractDataOutputHandler<Integer, LinkedHashMap<VariantData, Collection<VariantRunData>>>() {				
    			@Override
    			public Void call() {
    				for (VariantData variant : variantDataChunkMap.keySet())
    					try
    					{
    						String variantId = variant.getId();

    	                    if (variant.getReferencePosition() == null)
    	                    	unassignedMarkers.add(variantId);
    	                    // LOG.debug(marker + "\t" + (chromAndPos.length == 0 ? "0" : chromAndPos[0]) + "\t" + 0 + "\t" + (chromAndPos.length == 0 ? 0l : Long.parseLong(chromAndPos[1])) + LINE_SEPARATOR);
    	                    if (markerSynonyms != null) {
    	                    	String syn = markerSynonyms.get(variantId);
    	                        if (syn != null) {
    	                            variantId = syn;
    	                        }
    	                    }
    	                    snpFileWriter.write(variantId + "\t" + (variant.getReferencePosition() == null ? 0 : variant.getReferencePosition().getSequence()) + "\t" + 0 + "\t" + (variant.getReferencePosition() == null ? 0 : variant.getReferencePosition().getStartSite()) + LINE_SEPARATOR);

    	                    Map<String, List<String>> individualGenotypes = new TreeMap<String, List<String>>(new AlphaNumericComparator<String>());
    	                    Collection<VariantRunData> runs = variantDataChunkMap.get(variant);
    	                    if (runs != null) {
    	                        for (VariantRunData run : runs) {
    	    						for (Integer sampleId : run.getSampleGenotypes().keySet()) {
    	    							SampleGenotype sampleGenotype = run.getSampleGenotypes().get(sampleId);
    	                                String individualId = sampleIdToIndividualMap.get(sampleId);
    	                                List<String> storedIndividualGenotypes = individualGenotypes.get(individualId);
    	                                if (storedIndividualGenotypes == null) {
    	                                    storedIndividualGenotypes = new ArrayList<String>();
    	                                    individualGenotypes.put(individualId, storedIndividualGenotypes);
    	                                }
    	                                
    	    							if (!VariantData.gtPassesVcfAnnotationFilters(individualId, sampleGenotype, individuals1, annotationFieldThresholds, individuals2, annotationFieldThresholds2))
    	    								continue;	// skip genotype

    	                                storedIndividualGenotypes.add(sampleGenotype.getCode());
    	                            }
    	                        }
    	                    }

    	                    boolean fFirstLoopExecution = true;
    	                    if (individualGenotypes.isEmpty())
    	                    {
    	                    	byte[] missingDataBytes = "9".getBytes();
    	                    	for (int i=0; i<sortedIndividuals.size(); i++)
    	                    		zos.write(missingDataBytes);
    	                    	warningFileWriter.write("- No genotypes found for variant " + (variantId == null ? variant.getId() : variantId) + "\n");
    	                    }
    	                    else
    		                    for (String individualId : sortedIndividuals /* we use this object because it has the proper ordering*/) {
    		                        List<String> genotypes = individualGenotypes.get(individualId);
    		                        HashMap<Object, Integer> genotypeCounts = new HashMap<Object, Integer>(); // will help us to keep track of missing genotypes
    		                        int highestGenotypeCount = 0;
    		                        String mostFrequentGenotype = null;
    		                        if (genotypes != null) {
    		                            for (String genotype : genotypes) {
    		                                if (genotype == null) {
    		                                    continue; /* skip missing genotypes */
    		                                }
    		
    		                                int gtCount = 1 + Helper.getCountForKey(genotypeCounts, genotype);
    		                                if (gtCount > highestGenotypeCount) {
    		                                    highestGenotypeCount = gtCount;
    		                                    mostFrequentGenotype = genotype;
    		                                }
    		                                genotypeCounts.put(genotype, gtCount);
    		                            }
    		                        }
    		
    		                        List<String> alleles = mostFrequentGenotype == null ? new ArrayList<String>() : variant.getAllelesFromGenotypeCode(mostFrequentGenotype);
    		
    		                        int nOutputCode = 0;
    		                        if (mostFrequentGenotype == null)
    		                            nOutputCode = 9;
    		                        else
    		                            for (String all : Helper.split(mostFrequentGenotype, "/"))
    		                                if ("0".equals(all))
    		                                    nOutputCode++;
    		
    		                        if (fFirstLoopExecution && variant.getKnownAlleleList().size() > 2) {
    		                            warningFileWriter.write("- Variant " + variant.getId() + " is multi-allelic. Make sure Eigenstrat genotype encoding specifications are suitable for you.\n");
    		                        }
    		                        zos.write(("" + nOutputCode).getBytes());
    		
    		                        if (genotypeCounts.size() > 1 || alleles.size() > 2) {
    		                            if (genotypeCounts.size() > 1) {
    		                                warningFileWriter.write("- Dissimilar genotypes found for variant " + (variantId == null ? variant.getId() : variantId) + ", individual " + individualId + ". Exporting most frequent: " + nOutputCode + "\n");
    		                            }
    		                            if (alleles.size() > 2) {
    		                                warningFileWriter.write("- More than 2 alleles found for variant " + (variantId == null ? variant.getId() : variantId) + ", individual " + individualId + ". Exporting only the first 2 alleles.\n");
    		                            }
    		                        }
    		                        fFirstLoopExecution = false;
    		                    }
    		                    zos.write((LINE_SEPARATOR).getBytes());
    		                }
    					catch (Exception e)
    					{
    						if (progress.getError() == null)	// only log this once
    							LOG.debug("Unable to export " + variant.getId(), e);
    						progress.setError("Unable to export " + variant.getId() + ": " + e.getMessage());
    						
     	                    try
     	                    {
     	   	                    warningFileWriter.close();
								snpFileWriter.close();
							}
     	                    catch (IOException ignored) {}
    					}
    				return null;
    			}
    		};
    		
            int avgObjSize = (Integer) mongoTemplate.getCollection(mongoTemplate.getCollectionName(VariantRunData.class)).getStats().get("avgObjSize");
    		int nQueryChunkSize = Math.max(1, (nMaxChunkSizeInMb*1024*1024 / avgObjSize) / AsyncExportTool.INITIAL_NUMBER_OF_SIMULTANEOUS_QUERY_THREADS);
    		AsyncExportTool syncExportTool = new AsyncExportTool(markerCursor, markerCount, nQueryChunkSize, mongoTemplate, samplesToExport, dataOutputHandler, progress);
    		syncExportTool.launch();

    		while (progress.getCurrentStepProgress() < 100 && !progress.isAborted())
    			Thread.sleep(500);

            zos.closeEntry();
            snpFileWriter.close();
            
            if (unassignedMarkers.size() > 0)
            	LOG.info("No chromosomal position found for " + unassignedMarkers.size() + " markers " + StringUtils.join(unassignedMarkers, ", "));
            
            StringBuffer indFileContents = new StringBuffer();
            for (String individual : sortedIndividuals)
            {
            	String pop = MgdbDao.getIndividualPopulation(sModule, individual);
                indFileContents.append(individual + "\t" + getIndividualGenderCode(sModule, individual) + "\t" + (pop == null ? "." : pop) + LINE_SEPARATOR);
            }

            zos.putNextEntry(new ZipEntry(exportName + ".ind"));
            zos.write(indFileContents.toString().getBytes());
            zos.closeEntry();
            
            zos.putNextEntry(new ZipEntry(exportName + ".snp"));
            BufferedReader in = new BufferedReader(new FileReader(snpFile));
            String sLine;
            while ((sLine = in.readLine()) != null) {
                zos.write((sLine + "\n").getBytes());
            }
            in.close();
            zos.closeEntry();

            warningFileWriter.close();
            if (warningFile.length() > 0) {
                zos.putNextEntry(new ZipEntry(exportName + "-REMARKS.txt"));
                int nWarningCount = 0;
                in = new BufferedReader(new FileReader(warningFile));
                while ((sLine = in.readLine()) != null) {
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
        } finally {
            if (snpFile != null && snpFile.exists()) {
                snpFile.delete();
            }
        }
    }

    /* (non-Javadoc)
	 * @see fr.cirad.mgdb.exporting.IExportHandler#getStepList()
     */
    @Override
    public List<String> getStepList() {
        return Arrays.asList(new String[]{"Exporting data to EIGENSTRAT format"});
    }
    
	@Override
	public String[] getExportDataFileExtensions() {
		return new String[] {"snp", "ind", "eigenstratgeno"};
	}
}
