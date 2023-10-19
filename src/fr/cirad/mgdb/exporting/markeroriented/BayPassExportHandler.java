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
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Collectors;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;

import org.apache.commons.lang.StringUtils;
import org.apache.log4j.Logger;
import org.bson.Document;
import org.springframework.data.annotation.Id;
import org.springframework.data.mongodb.core.MongoTemplate;
import org.springframework.data.mongodb.core.query.Criteria;
import org.springframework.data.mongodb.core.query.Query;

import com.mongodb.client.MongoCollection;

import fr.cirad.mgdb.exporting.AbstractExportWritingThread;
import fr.cirad.mgdb.exporting.IExportHandler;
import fr.cirad.mgdb.exporting.tools.ExportManager;
import fr.cirad.mgdb.model.mongo.maintypes.Assembly;
import fr.cirad.mgdb.model.mongo.maintypes.GenotypingProject;
import fr.cirad.mgdb.model.mongo.maintypes.GenotypingSample;
import fr.cirad.mgdb.model.mongo.maintypes.Individual;
import fr.cirad.mgdb.model.mongo.maintypes.VariantData;
import fr.cirad.mgdb.model.mongo.maintypes.VariantRunData;
import fr.cirad.mgdb.model.mongo.subtypes.AbstractVariantData;
import fr.cirad.mgdb.model.mongo.subtypes.ReferencePosition;
import fr.cirad.mgdb.model.mongo.subtypes.SampleGenotype;
import fr.cirad.mgdb.model.mongodao.MgdbDao;
import fr.cirad.tools.AlphaNumericComparator;
import fr.cirad.tools.ProgressIndicator;
import fr.cirad.tools.mongo.MongoTemplateManager;
import htsjdk.variant.variantcontext.VariantContext.Type;

/**
 * The Class HapMapExportHandler.
 */
public class BayPassExportHandler extends AbstractMarkerOrientedExportHandler {

    /**
     * The Constant LOG.
     */
    private static final Logger LOG = Logger.getLogger(BayPassExportHandler.class);

    public static final String missingGenotype = "NN";
    
    /**
     * The supported variant types.
     */
    private static List<String> supportedVariantTypes;

    static {
        supportedVariantTypes = new ArrayList<String>();
        supportedVariantTypes.add(Type.SNP.toString());
//        supportedVariantTypes.add(Type.MNP.toString());
//        supportedVariantTypes.add(Type.INDEL.toString());
//        supportedVariantTypes.add(Type.MIXED.toString());
        supportedVariantTypes.add(Type.NO_VARIATION.toString());
    }

    /* (non-Javadoc)
	 * @see fr.cirad.mgdb.exporting.IExportHandler#getExportFormatName()
     */
    @Override
    public String getExportFormatName() {
        return "BAYPASS";
    }

    /* (non-Javadoc)
	 * @see fr.cirad.mgdb.exporting.IExportHandler#getExportFormatDescription()
     */
    @Override
    public String getExportFormatDescription() {
        return "Exports data in BAYPASS Format. See <a target='_blank' href='https://bitbucket.org/tasseladmin/tassel-5-source/wiki/UserManual/Load/Load'>https://bitbucket.org/tasseladmin/tassel-5-source/wiki/UserManual/Load/Load</a> for more details";
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
    public void exportData(OutputStream outputStream, String sModule, Integer nAssemblyId, String sExportingUser, Collection<String> individuals1, Collection<String> individuals2, ProgressIndicator progress, String tmpVarCollName, Document varQuery, long markerCount, Map<String, String> markerSynonyms, HashMap<String, Float> annotationFieldThresholds, HashMap<String, Float> annotationFieldThresholds2, List<GenotypingSample> samplesToExport, Collection<String> individualMetadataFieldsToExport, Map<String, InputStream> readyToExportFiles) throws Exception {
		Map<String, Integer> individualPositions = new LinkedHashMap<>();
        Map<String, String> allIndividualPops = new LinkedHashMap<>();//<individual,population>
		MongoTemplate mongoTemplate = MongoTemplateManager.get(sModule);
		Query query = new Query();
//        List<Individual> individualDataList = mongoTemplate.find(query, Individual.class);
        for (Individual individualData : mongoTemplate.find(query, Individual.class)) {
        	allIndividualPops.put(individualData.getId(), individualData.getPopulation());
        }
		for (String ind : samplesToExport.stream().map(gs -> gs.getIndividual()).distinct().sorted(new AlphaNumericComparator<String>()).collect(Collectors.toList()))
			individualPositions.put(ind, individualPositions.size());
		
        File warningFile = File.createTempFile("export_warnings_", "");
        FileWriter warningFileWriter = new FileWriter(warningFile);
        ZipOutputStream zos = IExportHandler.createArchiveOutputStream(outputStream, readyToExportFiles);
		MongoCollection collWithPojoCodec = mongoTemplate.getDb().withCodecRegistry(ExportManager.pojoCodecRegistry).getCollection(tmpVarCollName != null ? tmpVarCollName : mongoTemplate.getCollectionName(VariantRunData.class));
        Assembly assembly = mongoTemplate.findOne(new Query(Criteria.where("_id").is(nAssemblyId)), Assembly.class);
        String exportName = sModule + (assembly != null && assembly.getName() != null ? "__" + assembly.getName() : "") + "__" + markerCount + "variants__" + individualPositions.size() + "individuals";

        if (individualMetadataFieldsToExport != null && !individualMetadataFieldsToExport.isEmpty()) {
        	zos.putNextEntry(new ZipEntry(sModule + (assembly != null ? "__" + assembly.getName() : "") + "__" + individualPositions.size()+ "individuals_metadata.tsv"));
        	zos.write("individual".getBytes());
	        IExportHandler.writeMetadataFile(sModule, sExportingUser, individualPositions.keySet(), individualMetadataFieldsToExport, zos);
	    	zos.closeEntry();
        }
        final Map<Integer, String> sampleIdToIndividualMap = samplesToExport.stream().collect(Collectors.toMap(GenotypingSample::getId, sp -> sp.getIndividual()));
        List<String> samplePops = new ArrayList<>();
        zos.putNextEntry(new ZipEntry(exportName + ".popnames"));
       
        for (String individual : sampleIdToIndividualMap.values()) {
			String pop = allIndividualPops.get(individual);
			if (pop != null) {
                if (!samplePops.contains(pop.toString())) {
                	samplePops.add(pop.toString());
                }
            }
        }
        for( String pop : samplePops.stream().sorted().collect(Collectors.toList())) {
        	zos.write((pop + "\n").getBytes());

        }
//        for (GenotypingSample sample : samplesToExport) {
//			String pop = individualPops.get(sample.getIndividual());
//			if (pop != null) {
//                if (!samplePops.contains(pop.toString())) {
//                	samplePops.add(pop.toString());
//                	zos.write((pop + "\n").getBytes());
//                }
//            }
//        }
        zos.closeEntry();

        
        zos.putNextEntry(new ZipEntry(exportName + ".baypass"));

		
		final AtomicInteger initialStringBuilderCapacity = new AtomicInteger();

		int nQueryChunkSize = IExportHandler.computeQueryChunkSize(mongoTemplate, markerCount);
		AbstractExportWritingThread writingThread = new AbstractExportWritingThread() {
			public void run() {		
				final Iterator<String> exportedVariantIterator = orderedMarkerIDs.iterator();
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

//		                boolean fIsSNP = variant.getType().equals(Type.SNP.toString());
		                ReferencePosition rp = variant.getReferencePosition(nAssemblyId);
//						sb.append(idOfVarToWrite).append("\t").append(StringUtils.join(variant.getKnownAlleles(), "/") + "\t" + (rp == null ? 0 : rp.getSequence()) + "\t" + (rp == null ? 0 : rp.getStartSite()) + "\t" + "+\t" + (assembly == null ? "NA" : assembly.getName()) + "\tNA\tNA\tNA\tNA\tNA");

		                List<String>[] individualGenotypes = new ArrayList[individualPositions.size()];

		                if (runsToWrite != null)
		                	runsToWrite.forEach( run -> {
		                    	for (Integer sampleId : run.getSampleGenotypes().keySet()) {
	                                String individualId = sampleIdToIndividualMap.get(sampleId);
	                                Integer individualIndex = individualPositions.get(individualId);
	                                if (individualIndex == null)
	                                    continue;   // unwanted sample
	
									SampleGenotype sampleGenotype = run.getSampleGenotypes().get(sampleId);
		                            String gtCode = sampleGenotype.getCode();
		                            
									if (gtCode == null || !VariantData.gtPassesVcfAnnotationFilters(individualId, sampleGenotype, individuals1, annotationFieldThresholds, individuals2, annotationFieldThresholds2))
										continue;	// skip genotype
									
									if (individualGenotypes[individualIndex] == null)
										individualGenotypes[individualIndex] = new ArrayList<String>();
									individualGenotypes[individualIndex].add(gtCode);
		                        }
		                    });

		                int writtenGenotypeCount = 0;
		                
//		                HashMap<String, String> genotypeStringCache = new HashMap<>();
		        		Map<String, Integer> nbAllel0 = new LinkedHashMap<>();//<population,number>
		        		Map<String, Integer> nbAllel1 = new LinkedHashMap<>();

		                for (String individual : individualPositions.keySet() /* we use this list because it has the proper ordering */) {
		                    int individualIndex = individualPositions.get(individual);
		                    while (writtenGenotypeCount < individualIndex) {
//		                        sb.append(missingGenotype);
		                        writtenGenotypeCount++;
		                    }

		                	String mostFrequentGenotype = null;
		                    LinkedHashMap<Object, Integer> genotypeCounts = AbstractMarkerOrientedExportHandler.sortGenotypesFromMostFound(individualGenotypes[individualPositions.get(individual)]);
		                    if (genotypeCounts.size() == 1 || (genotypeCounts.values().stream().distinct().count() == 2)) {
		                        mostFrequentGenotype = genotypeCounts.keySet().iterator().next().toString();
		                    }
		                    if(mostFrequentGenotype != null) {
		                    	String population = allIndividualPops.get(individual);

		                    	Pattern pattern = Pattern.compile("(\\d)/(\\d)");
		                        Matcher matcher = pattern.matcher(mostFrequentGenotype);

		                        if (matcher.find()) {
		                        	  int gtA1 = Integer.parseInt(matcher.group(1));
			                          int gtA2 = Integer.parseInt(matcher.group(2));
			                          if (gtA1 == 0 && gtA2 == 0) {
//			                        	  int value = nbAllel0.getOrDefault(population, 0);
//			                        	  value += 2;
//			                        	  nbAllel0.put(population, value);
			                              nbAllel0.put(population, nbAllel0.getOrDefault(individual, 0) + 2);
			                          } else if (gtA1 == 0 || gtA2 == 0) {
			                              nbAllel0.put(population, nbAllel0.getOrDefault(population, 0) + 1);
			                              nbAllel1.put(population, nbAllel1.getOrDefault(population, 0) + 1);
			                          } else if (gtA1 == 1 && gtA2 == 1) {
			                              nbAllel1.put(population, nbAllel1.getOrDefault(population, 0) + 2);
			                          }
		                        }
		                        	
		                          

//		                        nbAllel0.put(population, nbAllel0.getOrDefault(population, 0) + 1);
		                    }
		                    	

//		                    String exportedGT = mostFrequentGenotype == null ? missingGenotype : genotypeStringCache.get(mostFrequentGenotype);
//		                    if (exportedGT == null) {
//		                    	exportedGT = StringUtils.join(variant.safelyGetAllelesFromGenotypeCode(mostFrequentGenotype, mongoTemplate), fIsSNP ? "" : "/");
//		                    	genotypeStringCache.put(mostFrequentGenotype, exportedGT);
//		                    }
		                    
		                    
		                    if (genotypeCounts.size() > 1)
		                    	warningFileWriter.write("- Dissimilar genotypes found for variant " + idOfVarToWrite + ", individual " + individual + ". " + (mostFrequentGenotype == null ? "Exporting as missing data" : "Exporting most frequent" ) + "\n");

//		                    sb.append("\t");
//		                    sb.append(exportedGT);
		                    writtenGenotypeCount++;
		                }
		                
		                
		                
		                
		                
		                
		                for (String pop : samplePops.stream().sorted().collect(Collectors.toList())) {
		                	if (!nbAllel0.isEmpty()) {
		                		int nb0 = nbAllel0.getOrDefault(pop, 0);

//		                		int nb0 = nbAllel0.get(pop);
		                		sb.append(nb0);
			                    sb.append("\t");
		                	}
		                	if (!nbAllel1.isEmpty()) {
		                		int nb1 = nbAllel1.getOrDefault(pop, 0);

//		                		int nb1 = nbAllel1.get(pop);
			                    sb.append(nb1);
			                    sb.append("\t");

		                	}
	
		                    
		                    
		                }
		                sb.append(LINE_SEPARATOR);

		                
		                
		                
		                
		                
		                
		                
		                
//		                while (writtenGenotypeCount < individualPositions.size()) {
//		                    sb.append(missingGenotype);
//		                    writtenGenotypeCount++;
//		                }
//		                sb.append(LINE_SEPARATOR);
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
		ExportManager exportManager = new ExportManager(mongoTemplate, nAssemblyId, collWithPojoCodec, VariantRunData.class, varQuery, samplesToExport, true, nQueryChunkSize, writingThread, markerCount, warningFileWriter, progress);
		exportManager.readAndWrite();
        zos.closeEntry();
        
        warningFileWriter.close();
        if (warningFile.length() > 0) {
            zos.putNextEntry(new ZipEntry(exportName + "-REMARKS.txt"));
            int nWarningCount = 0;
            BufferedReader in = new BufferedReader(new FileReader(warningFile));
            String sLine;
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
        
    }

    /* (non-Javadoc)
	 * @see fr.cirad.mgdb.exporting.IExportHandler#getStepList()
     */
    @Override
    public List<String> getStepList() {
        return Arrays.asList(new String[]{"Exporting data to HAPMAP format"});
    }
    
	@Override
	public String[] getExportDataFileExtensions() {
		return new String[] {"hapmap"};
	}

    @Override
    public int[] getSupportedPloidyLevels() {
        return null;
    }
}