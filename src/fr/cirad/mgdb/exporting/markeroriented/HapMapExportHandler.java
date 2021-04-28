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
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.stream.Collectors;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;

import org.apache.commons.lang.StringUtils;
import org.apache.log4j.Logger;
import org.bson.Document;
import org.springframework.data.mongodb.core.MongoTemplate;

import com.mongodb.client.MongoCollection;
import com.mongodb.client.MongoCursor;

import fr.cirad.mgdb.exporting.IExportHandler;
import fr.cirad.mgdb.model.mongo.maintypes.GenotypingSample;
import fr.cirad.mgdb.model.mongo.maintypes.VariantData;
import fr.cirad.mgdb.model.mongo.maintypes.VariantRunData;
import fr.cirad.mgdb.model.mongo.maintypes.VariantRunDataV2;
import fr.cirad.mgdb.model.mongo.subtypes.AbstractVariantData;
import fr.cirad.mgdb.model.mongo.subtypes.AbstractVariantDataV2;
import fr.cirad.mgdb.model.mongo.subtypes.ReferencePosition;
import fr.cirad.mgdb.model.mongo.subtypes.SampleGenotype;
import fr.cirad.mgdb.model.mongodao.MgdbDao;
import fr.cirad.tools.AlphaNumericComparator;
import fr.cirad.tools.Helper;
import fr.cirad.tools.ProgressIndicator;
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

    public static final String missingGenotype = "\tNN";
    
    /**
     * The supported variant types.
     */
    private static List<String> supportedVariantTypes;

    static {
        supportedVariantTypes = new ArrayList<String>();
        supportedVariantTypes.add(Type.SNP.toString());
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
        return "Exports data in HapMap Format. See <a target='_blank' href='http://heidi.chnebu.ch/doku.php?id=hapmap'>http://heidi.chnebu.ch/doku.php?id=hapmap</a> for more details";
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

    @Override
    public void exportData(OutputStream outputStream, String sModule, Integer nAssemblyId, Collection<GenotypingSample> samples1, Collection<GenotypingSample> samples2, ProgressIndicator progress, String tmpVarCollName, Document varQuery, long variantCount, Map<String, String> markerSynonyms, HashMap<String, Float> annotationFieldThresholds, HashMap<String, Float> annotationFieldThresholds2, List<GenotypingSample> samplesToExport, Map<String, InputStream> readyToExportFiles) throws Exception {
        
		List<String> individuals1 = MgdbDao.getIndividualsFromSamples(sModule, samples1).stream().map(ind -> ind.getId()).collect(Collectors.toList());	
		List<String> individuals2 = MgdbDao.getIndividualsFromSamples(sModule, samples2).stream().map(ind -> ind.getId()).collect(Collectors.toList());

		List<String> sortedIndividuals = samplesToExport.stream().map(gs -> gs.getIndividual()).distinct().sorted(new AlphaNumericComparator<String>()).collect(Collectors.toList());
		
		MongoTemplate mongoTemplate = MongoTemplateManager.get(sModule);
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

		boolean fV2Model = nAssemblyId == null || nAssemblyId < 0;
		MongoCollection collWithPojoCodec = mongoTemplate.getDb().withCodecRegistry(pojoCodecRegistry).getCollection(tmpVarCollName != null ? tmpVarCollName : mongoTemplate.getCollectionName(fV2Model ? VariantRunDataV2.class : VariantRunData.class));

		long markerCount = collWithPojoCodec.countDocuments(varQuery);
        String exportName = sModule + "__" + markerCount + "variants__" + sortedIndividuals.size() + "individuals";
        zos.putNextEntry(new ZipEntry(exportName + ".hapmap"));
        String header = "rs#" + "\t" + "alleles" + "\t" + "chrom" + "\t" + "pos" + "\t" + "strand" + "\t" + "assembly#" + "\t" + "center" + "\t" + "protLSID" + "\t" + "assayLSID" + "\t" + "panelLSID" + "\t" + "QCcode";
        zos.write(header.getBytes());
        for (String individual : sortedIndividuals)
            zos.write(("\t" + individual).getBytes());
        zos.write((LINE_SEPARATOR).getBytes());

		final Map<Integer, String> sampleIdToIndividualMap = new HashMap<>();
		for (GenotypingSample gs : samplesToExport)
			sampleIdToIndividualMap.put(gs.getId(), gs.getIndividual());
		
		Class runVersionType = fV2Model ? VariantRunDataV2.class : VariantRunData.class;

		Document projection = new Document();
		projection.append(!fV2Model ? (AbstractVariantData.FIELDNAME_REFERENCE_POSITION + "." + nAssemblyId) : AbstractVariantDataV2.FIELDNAME_REFERENCE_POSITION, 1);
		projection.append(!fV2Model ? AbstractVariantData.FIELDNAME_KNOWN_ALLELE_LIST : AbstractVariantDataV2.FIELDNAME_KNOWN_ALLELE_LIST, 1);
		projection.append(!fV2Model ? AbstractVariantData.FIELDNAME_TYPE : AbstractVariantDataV2.FIELDNAME_TYPE, 1);
		projection.append(!fV2Model ? AbstractVariantData.FIELDNAME_SYNONYMS : VariantRunDataV2.FIELDNAME_SYNONYMS, 1);
		projection.append(!fV2Model ? AbstractVariantData.FIELDNAME_ANALYSIS_METHODS : VariantRunDataV2.FIELDNAME_ANALYSIS_METHODS, 1);
		for (GenotypingSample sp : samplesToExport) {
			if (fV2Model)
				projection.append(VariantRunDataV2.FIELDNAME_SAMPLEGENOTYPES + "." + sp.getId() + "." + SampleGenotype.FIELDNAME_GENOTYPECODE, 1);
			else
				projection.append(VariantRunData.FIELDNAME_GENOTYPES + "." + sp.getId(), 1);
		}
		
		Number avgObjSize = (Number) mongoTemplate.getDb().runCommand(new Document("collStats", mongoTemplate.getCollectionName(VariantRunData.class))).get("avgObjSize");
		int nQueryChunkSize = (int) Math.max(1, (nMaxChunkSizeInMb*1024*1024 / avgObjSize.doubleValue()));
		
//		long b4 = System.currentTimeMillis();
		MongoCursor markerCursor = IExportHandler.getVariantCursorSortedWithCollation(mongoTemplate, collWithPojoCodec, runVersionType, varQuery, samplesToExport, false, nAssemblyId, nQueryChunkSize);
//		LOG.debug("cursor obtained in " + (System.currentTimeMillis() - b4) + "ms");
		LinkedHashMap<String, List<Object>> markerRunsToWrite = new LinkedHashMap<>();

		Thread writingThread = new Thread() {
			public void run() {
//				long b4 = System.currentTimeMillis();
//				String id = markerRunsToWrite.keySet().iterator().next();
//				System.out.println("writing " + id);
				
                HashMap<Object, Integer> genotypeCounts = new HashMap<Object, Integer>();	// will help us to keep track of missing genotypes
				for (String idOfVarToWrite : markerRunsToWrite.keySet()) {
					List<Object> runsToWrite = markerRunsToWrite.get(idOfVarToWrite);
					if (runsToWrite.isEmpty())
						continue;

					StringBuffer sb = new StringBuffer();
					try
					{
		                if (markerSynonyms != null) {
		                	String syn = markerSynonyms.get(idOfVarToWrite);
		                    if (syn != null)
		                    	idOfVarToWrite = syn;
		                }
	
		                Object variant = runsToWrite.get(0);
		                boolean fIsSNP = (nAssemblyId == null ? ((VariantRunDataV2) variant).getType() : ((VariantRunData) variant).getType()).equals(Type.SNP.toString());
	
		                ReferencePosition rp = nAssemblyId == null ? ((VariantRunDataV2) variant).getReferencePosition() : ((VariantRunData) variant).getReferencePosition(nAssemblyId);
		                sb.append(/*(variantId == null ? variant.getId() : */idOfVarToWrite/*)*/ + "\t" + StringUtils.join((nAssemblyId == null ? ((VariantRunDataV2) variant).getKnownAlleleList() : ((VariantRunData) variant).getKnownAlleleList()), "/") + "\t" + (rp == null ? 0 : rp.getSequence()) + "\t" + (rp == null ? 0 : rp.getStartSite()) + "\t" + "+\tNA\tNA\tNA\tNA\tNA\tNA");
	
		                Map<String, List<String>> individualGenotypes = new TreeMap<String, List<String>>(new AlphaNumericComparator<String>());
		                if (fV2Model) {
		                	for (Object vrd : runsToWrite) {
		                    	VariantRunDataV2 run = (VariantRunDataV2) vrd;
		                    	for (Integer sampleId : run.getSampleGenotypes().keySet()) {
									SampleGenotype sampleGenotype = run.getSampleGenotypes().get(sampleId);
									String individualId = sampleIdToIndividualMap.get(sampleId);
		                            
									if (!VariantData.gtPassesVcfAnnotationFiltersV2(individualId, sampleGenotype, individuals1, annotationFieldThresholds, individuals2, annotationFieldThresholds2))
										continue;	// skip genotype
									
		                            String gtCode = sampleGenotype.getCode().intern();
		                            List<String> storedIndividualGenotypes = individualGenotypes.get(individualId);
		                            if (storedIndividualGenotypes == null) {
		                                storedIndividualGenotypes = new ArrayList<String>();
		                                individualGenotypes.put(individualId, storedIndividualGenotypes);
		                            }
		                            storedIndividualGenotypes.add(gtCode);
		                        }
		                    }
		                }
		                else {
		                    for (Object vrd : runsToWrite) {
		                    	VariantRunData run = (VariantRunData) vrd;
								for (Integer sampleId : run.getGenotypes().keySet()) {
									String individualId = sampleIdToIndividualMap.get(sampleId);
		                            
									if (!VariantData.gtPassesVcfAnnotationFilters(individualId, sampleId, run.getMetadata(), individuals1, annotationFieldThresholds, individuals2, annotationFieldThresholds2))
										continue;	// skip genotype
									
		                            List<String> storedIndividualGenotypes = individualGenotypes.get(individualId);
		                            if (storedIndividualGenotypes == null) {
		                                storedIndividualGenotypes = new ArrayList<String>();
		                                individualGenotypes.put(individualId, storedIndividualGenotypes);
		                            }
		                            storedIndividualGenotypes.add(run.getGenotypes().get(sampleId).intern());
		                        }
		                    }
			            }
	
		                int writtenGenotypeCount = 0;
		                for (String individual : individualGenotypes.keySet() /* we use this list because it has the proper ordering */) {
		                    int individualIndex = sortedIndividuals.indexOf(individual);
		                    while (writtenGenotypeCount < individualIndex) {
		                        sb.append(missingGenotype);
		                        writtenGenotypeCount++;
		                    }
	
		                    List<String> genotypes = individualGenotypes.get(individual);
		                    genotypeCounts.clear();
		                    int highestGenotypeCount = 0;
		                    String mostFrequentGenotype = null;
		                    if (genotypes != null) {
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
		                    }
	
		                    String exportedGT = mostFrequentGenotype == null ? missingGenotype : ("\t" + StringUtils.join(nAssemblyId == null ? ((VariantRunDataV2) variant).getAllelesFromGenotypeCode(mostFrequentGenotype) : ((VariantRunData) variant).getAllelesFromGenotypeCode(mostFrequentGenotype), fIsSNP ? "" : "/"));
		                    sb.append(exportedGT);
		                    writtenGenotypeCount++;
	
		                    if (genotypeCounts.size() > 1)
		                        warningFileWriter.write("- Dissimilar genotypes found for variant " + /*(variantId == null ? variant.getId() : */idOfVarToWrite/*)*/ + ", individual " + individual + ". Exporting most frequent: " + new String(exportedGT) + "\n");
		                }
	
		                while (writtenGenotypeCount < sortedIndividuals.size()) {
		                    sb.append(missingGenotype);
		                    writtenGenotypeCount++;
		                }
		                sb.append(LINE_SEPARATOR);
	
			            zos.write(sb.toString().getBytes());
	                }
					catch (Exception e)
					{
						if (progress.getError() == null)	// only log this once
							LOG.debug("Unable to export " + idOfVarToWrite, e);
						progress.setError("Unable to export " + idOfVarToWrite + ": " + e.getMessage());
					}
				}
				markerRunsToWrite.clear();
//				long duration = System.currentTimeMillis() - b4;
//				System.out.println("wrote " + id + " in " + duration + "ms");
			}
		};

		IExportHandler.readAndWrite(markerCursor, progress, warningFileWriter, fV2Model, nQueryChunkSize, markerRunsToWrite, writingThread, variantCount);

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
}
