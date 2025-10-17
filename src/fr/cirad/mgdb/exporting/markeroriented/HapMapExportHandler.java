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

import java.io.File;
import java.io.IOException;
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

import fr.cirad.mgdb.model.mongo.maintypes.*;
import fr.cirad.mgdb.model.mongodao.MgdbDao;
import org.apache.commons.lang.StringUtils;
import org.apache.log4j.Logger;
import org.bson.Document;
import org.springframework.data.mongodb.core.MongoTemplate;
import org.springframework.data.mongodb.core.query.Criteria;
import org.springframework.data.mongodb.core.query.Query;

import com.mongodb.BasicDBList;

import fr.cirad.mgdb.exporting.IExportHandler;
import fr.cirad.mgdb.exporting.tools.ExportManager;
import fr.cirad.mgdb.exporting.tools.ExportManager.ExportOutputs;
import fr.cirad.mgdb.model.mongo.subtypes.AbstractVariantData;
import fr.cirad.mgdb.model.mongo.subtypes.Callset;
import fr.cirad.mgdb.model.mongo.subtypes.ReferencePosition;
import fr.cirad.mgdb.model.mongo.subtypes.SampleGenotype;
import fr.cirad.tools.Helper;
import fr.cirad.tools.ProgressIndicator;
import fr.cirad.tools.mgdb.VariantQueryWrapper;
import fr.cirad.tools.mongo.MongoTemplateManager;
import htsjdk.variant.variantcontext.VariantContext.Type;

/**
 * The Class HapMapExportHandler.
 */
public class HapMapExportHandler extends AbstractMarkerOrientedExportHandler {

    /**
     * The Constant LOG.
     */
    private static final Logger LOG = Logger.getLogger(HapMapExportHandler.class);

    public static final String missingGenotype = "NN";
    
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
        return "HAPMAP";
    }

    /* (non-Javadoc)
	 * @see fr.cirad.mgdb.exporting.IExportHandler#getExportFormatDescription()
     */
    @Override
    public String getExportFormatDescription() {
        return "Exports data in HapMap Format. See <a target='_blank' href='https://bitbucket.org/tasseladmin/tassel-5-source/wiki/UserManual/Load/Load'>https://bitbucket.org/tasseladmin/tassel-5-source/wiki/UserManual/Load/Load</a> for more details";
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
    public void exportData(OutputStream outputStream, String sModule, Integer nAssemblyId, String sExportingUser, ProgressIndicator progress, String tmpVarCollName, VariantQueryWrapper varQueryWrapper, long markerCount, Map<String, String> markerSynonyms, Map<String, Collection<String>> individualsByPop, boolean workWithSamples, Map<String, HashMap<String, Float>> annotationFieldThresholds, Collection<Callset> callSetsToExport, Collection<String> individualMetadataFieldsToExport, Map<String, InputStream> readyToExportFiles) throws Exception {
    	exportData(false, false, false, false, outputStream, sModule, nAssemblyId, sExportingUser, progress, tmpVarCollName, varQueryWrapper, markerCount, markerSynonyms, individualsByPop, workWithSamples, annotationFieldThresholds, callSetsToExport, individualMetadataFieldsToExport, readyToExportFiles);
    }

    public void exportData(boolean fSkipHapmapColumns, boolean fWriteAllelesAsIndexes, boolean fShowAlleleSeparator, boolean fEmptyStringForMissingData, OutputStream outputStream, String sModule, Integer nAssemblyId, String sExportingUser, ProgressIndicator progress, String tmpVarCollName, VariantQueryWrapper varQueryWrapper, long markerCount, Map<String, String> markerSynonyms, Map<String, Collection<String>> individualsByPop, boolean workWithSamples, Map<String, HashMap<String, Float>> annotationFieldThresholds, Collection<Callset> callSetsToExport, Collection<String> individualMetadataFieldsToExport, Map<String, InputStream> readyToExportFiles) throws Exception {
		List<String> materialNames = individualsByPop.values().stream().flatMap(Collection::stream).toList();
		
		MongoTemplate mongoTemplate = MongoTemplateManager.get(sModule);
        ZipOutputStream zos = IExportHandler.createArchiveOutputStream(outputStream, readyToExportFiles, null);
        Assembly assembly = mongoTemplate.findOne(new Query(Criteria.where("_id").is(nAssemblyId)), Assembly.class);
//        boolean workWithSamples = callSetsToExport.stream().filter(sp -> sp.isDetached()).count() == callSetsToExport.size();
        String exportName = IExportHandler.buildExportName(sModule, assembly, markerCount, materialNames.size(), workWithSamples);

        if (individualMetadataFieldsToExport == null || !individualMetadataFieldsToExport.isEmpty())
        	IExportHandler.addMetadataEntryIfAny(sModule + "__" + materialNames.size() + (workWithSamples ? "sample" : "individual" ) + "s_metadata.tsv", sModule, sExportingUser, materialNames, individualMetadataFieldsToExport, zos, (workWithSamples ? "sample" : "individual"), workWithSamples);

        zos.putNextEntry(new ZipEntry(exportName + ".hapmap"));

        final Map<Integer, String> callSetIdToIndividualMap = new HashMap<>();
        for (Callset cs : callSetsToExport)
        	callSetIdToIndividualMap.put(cs.getId(), workWithSamples ? cs.getSampleId() : cs.getIndividual());
        
        Collection<BasicDBList> variantDataQueries = varQueryWrapper.getVariantDataQueries();
        Document variantQueryForTargetCollection = variantDataQueries.isEmpty() ? new Document() : (tmpVarCollName == null ? new Document("$and", variantDataQueries.iterator().next()) : (varQueryWrapper.getBareQueries().iterator().hasNext() ? new Document("$and", varQueryWrapper.getBareQueries().iterator().next()) : new Document()));

		File[] warningFiles = writeGenotypeFile(fSkipHapmapColumns, fWriteAllelesAsIndexes, fShowAlleleSeparator, fEmptyStringForMissingData, zos, sModule, assembly, individualsByPop, workWithSamples, callSetIdToIndividualMap, annotationFieldThresholds, progress, tmpVarCollName, variantQueryForTargetCollection, markerCount, markerSynonyms, callSetsToExport).getWarningFiles();
		
        zos.closeEntry();
        
        IExportHandler.writeZipEntryFromChunkFiles(zos, warningFiles, exportName + "-REMARKS.txt");

        zos.finish();
        zos.close();
        progress.setCurrentStepProgress((short) 100);
    }
    
    public ExportOutputs writeGenotypeFile(boolean fSkipHapmapColumns, boolean fWriteAllesAsIndexes, boolean fShowAlleleSeparator, boolean fEmptyStringForMissingData, OutputStream os, String sModule, Assembly assembly, Map<String, Collection<String>> individualsByPop, boolean workWithSamples, Map<Integer, String> callSetIdToIndividualMap, Map<String, HashMap<String, Float>> annotationFieldThresholds, ProgressIndicator progress, String tmpVarCollName, Document variantQuery, long markerCount, Map<String, String> markerSynonyms, Collection<Callset> callSetsToExport) throws Exception {
    	MongoTemplate mongoTemplate = MongoTemplateManager.get(sModule);
		Integer projectId = null;
        for (Callset cs : callSetsToExport) {
            if (projectId == null)
                projectId = cs.getProjectId();
            else if (projectId != cs.getProjectId())
            {
                projectId = 0;
                break;	// more than one project are involved: no header will be written
            }
        }
		
		Map<String, Integer> individualPositions = IExportHandler.buildIndividualPositions(callSetsToExport, workWithSamples);

        os.write(((fSkipHapmapColumns ? "variant" : "rs#") + "\talleles\tchrom\tpos" + (fSkipHapmapColumns ? "" : "\tstrand\tassembly#\tcenter\tprotLSID\tassayLSID\tpanelLSID\tQCcode")).getBytes());
        for (String individual : individualPositions.keySet()) {
            os.write("\t".getBytes());
            os.write(individual.getBytes());
        }
        os.write((LINE_SEPARATOR).getBytes());        
        
		final AtomicInteger initialStringBuilderCapacity = new AtomicInteger();

		int nQueryChunkSize = IExportHandler.computeQueryChunkSize(mongoTemplate, markerCount);
		ExportManager.AbstractExportWriter exportWriter = new ExportManager.AbstractExportWriter(false) {
			@Override
			public void writeChunkRuns(Collection<Collection<VariantRunData>> markerRunsToWrite, List<String> orderedMarkerIDs, OutputStream genotypeOS, OutputStream variantOS, OutputStream warningOS) throws IOException {
	
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

		                boolean fIsSNP = variant.getType().equals(Type.SNP.toString());
		                ReferencePosition rp = variant.getReferencePosition(assembly == null ? null : assembly.getId());
		                sb.append(idOfVarToWrite).append("\t").append(StringUtils.join(variant.safelyGetKnownAlleles(mongoTemplate), "/") + "\t" + (rp == null ? 0 : rp.getSequence()) + "\t" + (rp == null ? 0 : rp.getStartSite()) + (fSkipHapmapColumns ? "" : ("\t" + "+\t" + (assembly == null ? "NA" : assembly.getName()) + "\tNA\tNA\tNA\tNA\tNA")));

		                List<String>[] individualGenotypes = new ArrayList[individualPositions.size()];

		                if (runsToWrite != null)
		                	runsToWrite.forEach( run -> {
		                    	for (Integer callSetId : run.getSampleGenotypes().keySet()) {
	                                String individualId = callSetIdToIndividualMap.get(callSetId);
	                                Integer individualIndex = individualPositions.get(individualId);
	                                if (individualIndex == null)
	                                    continue;   // unwanted sample
	
									SampleGenotype sampleGenotype = run.getSampleGenotypes().get(callSetId);
		                            String gtCode = sampleGenotype.getCode();
		                            
									if (gtCode == null || !VariantData.gtPassesVcfAnnotationFilters(individualId, sampleGenotype, individualsByPop, annotationFieldThresholds))
										continue;	// skip genotype
									
									if (individualGenotypes[individualIndex] == null)
										individualGenotypes[individualIndex] = new ArrayList<String>();
									individualGenotypes[individualIndex].add(gtCode);
		                        }
		                    });

		                int writtenGenotypeCount = 0;
		                
		                HashMap<String, String> genotypeStringCache = new HashMap<>();
		                for (String individual : individualPositions.keySet() /* we use this list because it has the proper ordering */) {
		                    int individualIndex = individualPositions.get(individual);
		                    while (writtenGenotypeCount < individualIndex) {
		                        sb.append(missingGenotype);
		                        writtenGenotypeCount++;
		                    }

		                	String mostFrequentGenotype = null;
		                    LinkedHashMap<Object, Integer> genotypeCounts = AbstractMarkerOrientedExportHandler.sortGenotypesFromMostFound(individualGenotypes[individualPositions.get(individual)]);
                            if (genotypeCounts.size() == 1 || genotypeCounts.values().stream().limit(2).distinct().count() == 2)
                            	mostFrequentGenotype = genotypeCounts.keySet().iterator().next().toString();

		                    String exportedGT = mostFrequentGenotype == null ? (fEmptyStringForMissingData ? "" : missingGenotype) : genotypeStringCache.get(mostFrequentGenotype);
		                    if (exportedGT == null) {
		                    	exportedGT = fWriteAllesAsIndexes ? mostFrequentGenotype : StringUtils.join(variant.safelyGetAllelesFromGenotypeCode(mostFrequentGenotype, mongoTemplate).stream().map(all -> "N".equals(all) ? "D" : ("NN".equals(all) ? "I" : all)).toList(), !fShowAlleleSeparator && fIsSNP ? "" : "/");
		                    	genotypeStringCache.put(mostFrequentGenotype, exportedGT);
		                    }
		                    
		                    if (genotypeCounts.size() > 1)
		                    	warningOS.write(("- Dissimilar genotypes found for variant " + idOfVarToWrite + ", individual " + individual + ". " + (mostFrequentGenotype == null ? "Exporting as missing data" : "Exporting most frequent: " + exportedGT) + "\n").getBytes());

		                    sb.append("\t");
		                    sb.append(exportedGT);
		                    writtenGenotypeCount++;
		                }
		                
		                while (writtenGenotypeCount < individualPositions.size()) {
		                    sb.append(missingGenotype);
		                    writtenGenotypeCount++;
		                }
		                sb.append(LINE_SEPARATOR);
		                if (initialStringBuilderCapacity.get() == 0)
		                    initialStringBuilderCapacity.set(sb.length());
			            genotypeOS.write(sb.toString().getBytes());
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

        String usedCollName = tmpVarCollName != null ? tmpVarCollName : mongoTemplate.getCollectionName(VariantRunData.class);
		ExportManager exportManager = new ExportManager(sModule, assembly == null ? null : assembly.getId(), mongoTemplate.getDb().withCodecRegistry(ExportManager.pojoCodecRegistry).getCollection(usedCollName), VariantRunData.class, variantQuery, callSetsToExport, true, nQueryChunkSize, exportWriter, markerCount, progress);
		if (tmpFolderPath != null)
			exportManager.setTmpExtractionFolder(tmpFolderPath + File.separator + "all_users" + File.separator + Helper.convertToMD5(progress.getProcessId()));
		exportManager.readAndWrite(os);
		return exportManager.getOutputs();
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