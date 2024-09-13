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
import java.util.TreeMap;
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
import fr.cirad.mgdb.model.mongo.maintypes.VariantData;
import fr.cirad.mgdb.model.mongo.maintypes.VariantRunData;
import fr.cirad.mgdb.model.mongo.subtypes.AbstractVariantData;
import fr.cirad.mgdb.model.mongo.subtypes.ReferencePosition;
import fr.cirad.mgdb.model.mongo.subtypes.SampleGenotype;
import fr.cirad.tools.AlphaNumericComparator;
import fr.cirad.tools.ProgressIndicator;
import fr.cirad.tools.mgdb.VariantQueryWrapper;
import fr.cirad.tools.mongo.MongoTemplateManager;
import htsjdk.variant.variantcontext.VariantContext.Type;

/**
 * The Class GFFExportHandler.
 */
public class GFFExportHandler extends AbstractMarkerOrientedExportHandler {

    /**
     * The Constant LOG.
     */
    private static final Logger LOG = Logger.getLogger(GFFExportHandler.class);

    /* (non-Javadoc)
	 * @see fr.cirad.mgdb.exporting.IExportHandler#getExportFormatName()
     */
    @Override
    public String getExportFormatName() {
        return "GFF3";
    }

    /* (non-Javadoc)
	 * @see fr.cirad.mgdb.exporting.IExportHandler#getExportFormatDescription()
     */
    @Override
    public String getExportFormatDescription() {
        return "Exports data in GFF3 Format based on Sequence Ontology. See <a target='_blank' href='https://m.ensembl.org/info/website/upload/gff3.html'>https://m.ensembl.org/info/website/upload/gff3.html</a> and <a target='_blank' href='https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md'>https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md</a>";
    }
    
	@Override
	public String getExportArchiveExtension() {
		return "zip";
	}

    @Override
    public void exportData(OutputStream outputStream, String sModule, Integer nAssemblyId, String sExportingUser, ProgressIndicator progress, String tmpVarCollName, VariantQueryWrapper varQueryWrapper, long markerCount, Map<String, String> markerSynonyms, Map<String, Collection<String>> individuals, Map<String, HashMap<String, Float>> annotationFieldThresholds, List<GenotypingSample> samplesToExport, Collection<String> individualMetadataFieldsToExport, Map<String, InputStream> readyToExportFiles) throws Exception {
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

        zos.putNextEntry(new ZipEntry(exportName + ".gff3"));
        String header = "##gff-version 3" + LINE_SEPARATOR;
        zos.write(header.getBytes());

        TreeMap<String, String> typeToOntology = new TreeMap<>();
        typeToOntology.put(Type.SNP.toString(), "SO:0000694");
        typeToOntology.put(Type.INDEL.toString(), "SO:1000032");
        typeToOntology.put(Type.MIXED.toString(), "SO:0001059");
        typeToOntology.put(Type.SYMBOLIC.toString(), "SO:0000109");
        typeToOntology.put(Type.MNP.toString(), "SO:0001059");

        final Map<Integer, String> sampleIdToIndividualMap = samplesToExport.stream().collect(Collectors.toMap(GenotypingSample::getId, sp -> sp.getIndividual()));
		
        ArrayList<Comparable> unassignedMarkers = new ArrayList<>();
		int nQueryChunkSize = IExportHandler.computeQueryChunkSize(mongoTemplate, markerCount);
		final AtomicInteger initialStringBuilderCapacity = new AtomicInteger();
		
		ExportManager.AbstractExportWriter writingThread = new ExportManager.AbstractExportWriter() {
			@Override
			public void writeChunkRuns(Collection<Collection<VariantRunData>> markerRunsToWrite, List<String> orderedMarkerIDs, OutputStream mainOS, OutputStream warningOS) {	
				final Iterator<String> exportedVariantIterator = orderedMarkerIDs.iterator();
				for (Collection<VariantRunData> runsToWrite : markerRunsToWrite) {
                	String idOfVarToWrite = exportedVariantIterator.next();
					if (progress.isAborted() || progress.getError() != null)
						return;
					
					AbstractVariantData variant = runsToWrite == null || runsToWrite.isEmpty() ? mongoTemplate.findById(idOfVarToWrite, VariantData.class) : runsToWrite.iterator().next();

					List<String> variantDataOrigin = new ArrayList<>();
					StringBuilder sb = new StringBuilder(initialStringBuilderCapacity.get() == 0 ? 3 * individualPositions.size() /* rough estimation */ : initialStringBuilderCapacity.get());
					try
					{
		                if (markerSynonyms != null) {
		                	String syn = markerSynonyms.get(idOfVarToWrite);
		                    if (syn != null)
		                    	idOfVarToWrite = syn;
		                }

		                List<String>[] individualGenotypes = new ArrayList[individualPositions.size()];
		                if (runsToWrite != null)
		                	for (VariantRunData run : runsToWrite) {
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
		                    }

	                	String refAllele;
						try {
							refAllele = variant.getKnownAlleles().iterator().next();
						}
						catch (IndexOutOfBoundsException ioobe) {	// VariantRunData's known-allele-list is not up to date
							variant.setKnownAlleles(mongoTemplate.findById(variant.getVariantId(), VariantData.class).getKnownAlleles());
							mongoTemplate.save(variant);
							refAllele = variant.getKnownAlleles().iterator().next();
						}
						
						ReferencePosition rp = variant.getReferencePosition(nAssemblyId);
		                String chrom = rp == null ? "0" : rp.getSequence();
		                long start = rp == null ? 0 : rp.getStartSite();
		                long end = Type.SNP.toString().equals(variant.getType()) ? start : start + refAllele.length() - 1;
		                sb.append(chrom).append("\t").append(StringUtils.join(variantDataOrigin, ";") /*source*/).append("\t").append(typeToOntology.get(variant.getType())).append("\t").append(start).append("\t").append(end).append("\t.\t+\t.\t");
		                String syn = markerSynonyms == null ? null : markerSynonyms.get(variant.getVariantId());
		                sb.append("ID=").append(idOfVarToWrite).append(";").append((syn != null ? "Name=" + syn + ";" : "")).append("alleles=").append(StringUtils.join(variant.getKnownAlleles(), "/")).append(";").append("refallele=").append(refAllele).append(";sample_genotypes=");

		                boolean fFirstIndividual = true;
		                for (String individual : individualPositions.keySet() /* we use this list because it has the proper ordering */) {      
		                	String mostFrequentGenotype = null;
		                    LinkedHashMap<Object, Integer> genotypeCounts = AbstractMarkerOrientedExportHandler.sortGenotypesFromMostFound(individualGenotypes[individualPositions.get(individual)]);
                            if (genotypeCounts.size() == 1 || genotypeCounts.values().stream().limit(2).distinct().count() == 2)
                            	mostFrequentGenotype = genotypeCounts.keySet().iterator().next().toString();

		                    String exportedGT = mostFrequentGenotype == null ? "." : StringUtils.join(variant.safelyGetAllelesFromGenotypeCode(mostFrequentGenotype, mongoTemplate).stream().map(all -> "N".equals(all) ? "D" : ("NN".equals(all) ? "I" : all)).toList(), "/");

		                    if (genotypeCounts.size() > 1)
		                    	warningOS.write(("- Dissimilar genotypes found for variant " + idOfVarToWrite + ", individual " + individual + ". " + (mostFrequentGenotype == null ? "Exporting as missing data" : "Exporting most frequent: " + exportedGT) + "\n").getBytes());
		                    
		                    sb.append(fFirstIndividual ? "" : ",").append(individual).append(":").append(exportedGT);
		                    fFirstIndividual = false;
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
				}
			}
		};
		
		Collection<BasicDBList> variantRunDataQueries = varQueryWrapper.getVariantRunDataQueries();
		ExportManager exportManager = new ExportManager(sModule, nAssemblyId, collWithPojoCodec, VariantRunData.class, !variantRunDataQueries.isEmpty() ? variantRunDataQueries.iterator().next() : new BasicDBList(), samplesToExport, true, nQueryChunkSize, writingThread, markerCount, progress);
		File[] warningFiles = exportManager.readAndWrite(zos);	
        zos.closeEntry();
        
        if (unassignedMarkers.size() > 0)
        	LOG.info("No chromosomal position found for " + unassignedMarkers.size() + " markers " + StringUtils.join(unassignedMarkers, ", "));

        IExportHandler.writeWarnings(zos, warningFiles, exportName);

        zos.finish();
        zos.close();
        progress.setCurrentStepProgress((short) 100);
    }


    /* (non-Javadoc)
     * @see fr.cirad.mgdb.exporting.IExportHandler#getStepList()
     */
    @Override
    public List<String> getStepList() {
        return Arrays.asList(new String[]{"Exporting data to GFF3 format"});
    }
    
	@Override
	public String[] getExportDataFileExtensions() {
		return new String[] {"gff3"};
	}

    @Override
    public int[] getSupportedPloidyLevels() {
        return null;
    }
}
