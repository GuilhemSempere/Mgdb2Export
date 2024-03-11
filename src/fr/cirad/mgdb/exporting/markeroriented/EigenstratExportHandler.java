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
package fr.cirad.mgdb.exporting.markeroriented;

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
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.Collectors;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;

import org.apache.commons.lang.StringUtils;
import org.apache.log4j.Logger;
import org.springframework.data.mongodb.core.MongoTemplate;
import org.springframework.data.mongodb.core.query.Criteria;
import org.springframework.data.mongodb.core.query.Query;

import com.mongodb.BasicDBList;
import com.mongodb.client.MongoCollection;

import fr.cirad.mgdb.exporting.AbstractExportWritingThread;
import fr.cirad.mgdb.exporting.IExportHandler;
import fr.cirad.mgdb.exporting.tools.ExportManager;
import fr.cirad.mgdb.model.mongo.maintypes.Assembly;
import fr.cirad.mgdb.model.mongo.maintypes.GenotypingSample;
import fr.cirad.mgdb.model.mongo.maintypes.Individual;
import fr.cirad.mgdb.model.mongo.maintypes.VariantData;
import fr.cirad.mgdb.model.mongo.maintypes.VariantRunData;
import fr.cirad.mgdb.model.mongo.subtypes.AbstractVariantData;
import fr.cirad.mgdb.model.mongo.subtypes.ReferencePosition;
import fr.cirad.mgdb.model.mongo.subtypes.SampleGenotype;
import fr.cirad.mgdb.model.mongodao.MgdbDao;
import fr.cirad.tools.AlphaNumericComparator;
import fr.cirad.tools.Helper;
import fr.cirad.tools.ProgressIndicator;
import fr.cirad.tools.mgdb.VariantQueryWrapper;
import fr.cirad.tools.mongo.MongoTemplateManager;
import htsjdk.variant.variantcontext.VariantContext.Type;

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
    
	public static final byte missingData = 9;

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

    /* (non-Javadoc)
	 * @see fr.cirad.mgdb.exporting.markeroriented.AbstractMarkerOrientedExportHandler#exportData(java.io.OutputStream, java.lang.String, java.util.List, fr.cirad.tools.ProgressIndicator, com.mongodb.DBCursor, java.util.Map, int, int, java.util.Map)
     */
    @Override
	public void exportData(OutputStream outputStream, String sModule, Integer nAssemblyId, String sExportingUser, ProgressIndicator progress, String tmpVarCollName, VariantQueryWrapper varQueryWrapper, long markerCount, Map<String, String> markerSynonyms, Map<String, Collection<String>> individuals, Map<String, HashMap<String, Float>> annotationFieldThresholds, List<GenotypingSample> samplesToExport, Collection<String> individualMetadataFieldsToExport, Map<String, InputStream> readyToExportFiles) throws Exception {
        File warningFile = File.createTempFile("export_warnings_", "");
        FileWriter warningFileWriter = new FileWriter(warningFile);
        File snpFile = null;

        try {
            snpFile = File.createTempFile("snpFile", "");
            FileWriter snpFileWriter = new FileWriter(snpFile);

            MongoTemplate mongoTemplate = MongoTemplateManager.get(sModule);
            ZipOutputStream zos = IExportHandler.createArchiveOutputStream(outputStream, readyToExportFiles);
    		MongoCollection collWithPojoCodec = mongoTemplate.getDb().withCodecRegistry(ExportManager.pojoCodecRegistry).getCollection(tmpVarCollName != null ? tmpVarCollName : mongoTemplate.getCollectionName(VariantRunData.class));

    		Map<String, Integer> individualPositions = new LinkedHashMap<>();
    		for (String ind : samplesToExport.stream().map(gs -> gs.getIndividual()).distinct().sorted(new AlphaNumericComparator<String>()).collect(Collectors.toList()))
    			individualPositions.put(ind, individualPositions.size());
 
    		Assembly assembly = mongoTemplate.findOne(new Query(Criteria.where("_id").is(nAssemblyId)), Assembly.class);
            String exportName = sModule + (assembly != null && assembly.getName() != null ? "__" + assembly.getName() : "") + "__" + markerCount + "variants__" + individualPositions.size() + "individuals";
            
            if (individualMetadataFieldsToExport == null || !individualMetadataFieldsToExport.isEmpty())
            	IExportHandler.addMetadataEntryIfAny(sModule + "__" + individualPositions.size() + "individuals_metadata.tsv", sModule, sExportingUser, individualPositions.keySet(), individualMetadataFieldsToExport, zos, "individual");
            
            zos.putNextEntry(new ZipEntry(exportName + ".eigenstratgeno"));
            final Map<Integer, String> sampleIdToIndividualMap = samplesToExport.stream().collect(Collectors.toMap(GenotypingSample::getId, sp -> sp.getIndividual()));
            ArrayList<Comparable> unassignedMarkers = new ArrayList<>();

    		int nQueryChunkSize = IExportHandler.computeQueryChunkSize(mongoTemplate, markerCount);
    	    final AtomicInteger initialStringBuilderCapacity = new AtomicInteger();
    		
    		AbstractExportWritingThread writingThread = new AbstractExportWritingThread() {
    			public void run() {
    				final Iterator<String> exportedVariantIterator = orderedMarkerIDs.iterator();
//    				final Iterator<String> exportedVariantIterator = orderedMarkerIDs != null ? orderedMarkerIDs.iterator() : markerRunsToWrite.;
    				markerRunsToWrite.forEach(runsToWrite -> {
                    	String idOfVarToWrite = exportedVariantIterator.next();
    					if (progress.isAborted() || progress.getError() != null)
    						return;
    					
    					AbstractVariantData variant = runsToWrite == null || runsToWrite.isEmpty() ? mongoTemplate.findById(idOfVarToWrite, VariantData.class) : runsToWrite.iterator().next();

    					StringBuilder sb = new StringBuilder(initialStringBuilderCapacity.get() == 0 ? 3 * individualPositions.size() /* rough estimation */ : initialStringBuilderCapacity.get());
    					try
    					{
    		                if (markerSynonyms != null) {
    		                	String syn = markerSynonyms.get(idOfVarToWrite);
    		                    if (syn != null)
    		                    	idOfVarToWrite = syn;
    		                }

    		                ReferencePosition rp = variant.getReferencePosition(nAssemblyId);
    	                    snpFileWriter.write(idOfVarToWrite + "\t" + (rp == null ? 0 : rp.getSequence()) + "\t" + 0 + "\t" + (rp == null ? 0 : rp.getStartSite()) + LINE_SEPARATOR);

    		                List<String>[] individualGenotypes = new ArrayList[individualPositions.size()];
    		                if (runsToWrite != null)
	    		                runsToWrite.forEach(run -> {
	    	                    	for (Integer sampleId : run.getSampleGenotypes().keySet()) {
	                                    String individualId = sampleIdToIndividualMap.get(sampleId);
	                                    Integer individualIndex = individualPositions.get(individualId);
	                                    if (individualIndex == null)
	                                        continue;   // unwanted sample
	
	    								SampleGenotype sampleGenotype = run.getSampleGenotypes().get(sampleId);
	    	                            String gtCode = sampleGenotype.getCode();
	    	                            
	    								if (gtCode == null || !VariantData.gtPassesVcfAnnotationFilters(individualId, sampleGenotype, individuals, annotationFieldThresholds))
	    									continue;	// skip genotype
	
	    								if (individualGenotypes[individualIndex] == null)
	    									individualGenotypes[individualIndex] = new ArrayList<String>();
	    								individualGenotypes[individualIndex].add(gtCode);
	    	                        }
	    	                    });

    	                    boolean fFirstLoopExecution = true;
    		                for (String individual : individualPositions.keySet() /* we use this list because it has the proper ordering */) {
    		                	String mostFrequentGenotype = null;
    		                    LinkedHashMap<Object, Integer> genotypeCounts = sortGenotypesFromMostFound(individualGenotypes[individualPositions.get(individual)]);
                                if (genotypeCounts.size() == 1 || genotypeCounts.values().stream().limit(2).distinct().count() == 2)
                                	mostFrequentGenotype = genotypeCounts.keySet().iterator().next().toString();

                                long nAlleleCount = 0;
                                byte nOutputCode = 0;
                                if (mostFrequentGenotype == null)
                                    nOutputCode = missingData;
                                else {
                                    for (String all : Helper.split(mostFrequentGenotype, "/")) {
                                        nAlleleCount++;
                                        if ("0".equals(all))
                                            nOutputCode++;
                                    }
                                }

                                if (genotypeCounts.size() > 1)
                                	warningFileWriter.write("- Dissimilar genotypes found for variant " + idOfVarToWrite + ", individual " + individual + ". " + (mostFrequentGenotype == null ? "Exporting as missing data" : "Exporting most frequent: " + nOutputCode) + "\n");
                                if (nAlleleCount > 2)
                                    warningFileWriter.write("- More than 2 alleles found for variant " + idOfVarToWrite + ", individual " + individual + ". Exporting only the first 2 alleles.\n");

                                if (fFirstLoopExecution && variant.getKnownAlleles().size() > 2)
                                    warningFileWriter.write("- Variant " + variant.getVariantId() + " is multi-allelic. Make sure Eigenstrat genotype encoding specifications are suitable for you.\n");

                                sb.append(nOutputCode);        
                                fFirstLoopExecution = false;
		                    }
		                    sb.append(LINE_SEPARATOR);
	                        if (initialStringBuilderCapacity.get() == 0)
	                            initialStringBuilderCapacity.set(sb.length());
	        				zos.write(sb.toString().getBytes());
    	                }
    					catch (Exception e)
    					{
    						if (progress.getError() == null)	// only log this once
    							LOG.error("Unable to export " + idOfVarToWrite, e);
    						progress.setError("Unable to export " + idOfVarToWrite + ": " + e.getMessage());
    					}
    				});
    			}
    		};

    		Collection<BasicDBList> variantRunDataQueries = varQueryWrapper.getVariantRunDataQueries();
    		ExportManager exportManager = new ExportManager(mongoTemplate, nAssemblyId, collWithPojoCodec, VariantRunData.class, !variantRunDataQueries.isEmpty() ? variantRunDataQueries.iterator().next() : new BasicDBList(), samplesToExport, true, nQueryChunkSize, writingThread, markerCount, warningFileWriter, progress);
    		exportManager.readAndWrite();
            zos.closeEntry();            
            
            if (unassignedMarkers.size() > 0)
            	LOG.info("No chromosomal position found for " + unassignedMarkers.size() + " markers " + StringUtils.join(unassignedMarkers, ", "));
            
            progress.addStep("Generating .ind file");
            progress.moveToNextStep();
            StringBuilder indFileContents = new StringBuilder(individualPositions.size() * 10);
            Map<String, String> individualPops = IExportHandler.getIndividualPopulations(individuals, true);
            for (String individual : individualPositions.keySet()) {
            	String pop = individualPops.get(individual);
                indFileContents.append(individual).append("\t").append(getIndividualGenderCode(sModule, individual)).append("\t").append((pop == null ? "." : pop)).append(LINE_SEPARATOR);
            }

            zos.putNextEntry(new ZipEntry(exportName + ".ind"));
            zos.write(indFileContents.toString().getBytes());
            zos.closeEntry();
            
            snpFileWriter.close();
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
	
    @Override
    public int[] getSupportedPloidyLevels() {
        return new int[] {2};
    }
}
