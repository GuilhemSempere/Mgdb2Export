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
import org.bson.Document;
import org.springframework.data.mongodb.core.MongoTemplate;

import com.mongodb.client.MongoCollection;
import com.mongodb.client.MongoCursor;

import fr.cirad.mgdb.exporting.tools.AsyncExportTool;
import fr.cirad.mgdb.exporting.tools.AsyncExportTool.AbstractDataOutputHandler;
import fr.cirad.mgdb.model.mongo.maintypes.GenotypingSample;
import fr.cirad.mgdb.model.mongo.maintypes.VariantData;
import fr.cirad.mgdb.model.mongo.maintypes.VariantRunData;
import fr.cirad.mgdb.model.mongo.subtypes.ReferencePosition;
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
    public void exportData(OutputStream outputStream, String sModule, int nAssemblyId, Collection<GenotypingSample> samples1, Collection<GenotypingSample> samples2, ProgressIndicator progress, MongoCollection<Document> varColl, Document varQuery, Map<String, String> markerSynonyms, HashMap<String, Float> annotationFieldThresholds, HashMap<String, Float> annotationFieldThresholds2, List<GenotypingSample> samplesToExport, Map<String, InputStream> readyToExportFiles) throws Exception {
        
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

		long markerCount = varColl.countDocuments(varQuery);
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
		
		AbstractDataOutputHandler<Integer, LinkedHashMap<VariantData, Collection<VariantRunData>>> dataOutputHandler = new AbstractDataOutputHandler<Integer, LinkedHashMap<VariantData, Collection<VariantRunData>>>() {				
			@Override
			public Void call() {
				StringBuffer sb = new StringBuffer();
				for (VariantData variant : variantDataChunkMap.keySet())
					try
					{
						String variantId = variant.getId();
		                if (markerSynonyms != null) {
		                	String syn = markerSynonyms.get(variantId);
		                    if (syn != null) {
		                        variantId = syn;
		                    }
		                }
	
		                boolean fIsSNP = variant.getType().equals(Type.SNP.toString());
	
		                ReferencePosition rp = variant.getReferencePosition(nAssemblyId);
		                sb.append((variantId == null ? variant.getId() : variantId) + "\t" + StringUtils.join(variant.getKnownAlleleList(), "/") + "\t" + (rp == null ? 0 : rp.getSequence()) + "\t" + (rp == null ? 0 : rp.getStartSite()) + "\t" + "+\tNA\tNA\tNA\tNA\tNA\tNA");
	
		                Map<String, List<String>> individualGenotypes = new TreeMap<String, List<String>>(new AlphaNumericComparator<String>());
		                Collection<VariantRunData> runs = variantDataChunkMap.get(variant);
		                if (runs != null) {
		                    for (VariantRunData run : runs) {
								for (Integer sampleId : run.getGenotypes().keySet()) {
									String individualId = sampleIdToIndividualMap.get(sampleId);
		                            
									if (!VariantData.gtPassesVcfAnnotationFilters(individualId, sampleId, run.getMetadata(), individuals1, annotationFieldThresholds, individuals2, annotationFieldThresholds2))
										continue;	// skip genotype
									
		                            List<String> storedIndividualGenotypes = individualGenotypes.get(individualId);
		                            if (storedIndividualGenotypes == null) {
		                                storedIndividualGenotypes = new ArrayList<String>();
		                                individualGenotypes.put(individualId, storedIndividualGenotypes);
		                            }
		                            storedIndividualGenotypes.add(run.getGenotypes().get(sampleId));
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
		                    HashMap<Object, Integer> genotypeCounts = new HashMap<Object, Integer>();	// will help us to keep track of missing genotypes
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
	
		                    String exportedGT = mostFrequentGenotype == null ? missingGenotype : ("\t" + StringUtils.join(variant.getAllelesFromGenotypeCode(mostFrequentGenotype), fIsSNP ? "" : "/"));
		                    sb.append(exportedGT);
		                    writtenGenotypeCount++;
	
		                    if (genotypeCounts.size() > 1) {
		                        warningFileWriter.write("- Dissimilar genotypes found for variant " + (variantId == null ? variant.getId() : variantId) + ", individual " + individual + ". Exporting most frequent: " + new String(exportedGT) + "\n");
		                    }
		                }
	
		                while (writtenGenotypeCount < sortedIndividuals.size()) {
		                    sb.append(missingGenotype);
		                    writtenGenotypeCount++;
		                }
		                sb.append(LINE_SEPARATOR);
	
			            if (progress.isAborted()) {
			                warningFileWriter.close();
			                return null;
			            }
	                }
					catch (Exception e)
					{
						if (progress.getError() == null)	// only log this once
							LOG.debug("Unable to export " + variant.getId(), e);
						progress.setError("Unable to export " + variant.getId() + ": " + e.getMessage());
					}
				
                try
                {
    				zos.write(sb.toString().getBytes());
				}
                catch (IOException ioe)
                {
                	progress.setError("Unable to export data for " + variantDataChunkMap.keySet() + ": " + ioe.getMessage());
                }
				return null;
			}
		};
	
		Number avgObjSize = (Number) mongoTemplate.getDb().runCommand(new Document("collStats", mongoTemplate.getCollectionName(VariantRunData.class))).get("avgObjSize");
		int nQueryChunkSize = (int) Math.max(1, (nMaxChunkSizeInMb*1024*1024 / avgObjSize.doubleValue()) / AsyncExportTool.WRITING_QUEUE_CAPACITY);
		try (MongoCursor<Document> markerCursor = varColl.find(varQuery).projection(projectionDoc(nAssemblyId)).sort(sortDoc(nAssemblyId)).noCursorTimeout(true).collation(collationObj).batchSize(nQueryChunkSize).iterator()) {
			AsyncExportTool syncExportTool = new AsyncExportTool(markerCursor, markerCount, nQueryChunkSize, mongoTemplate, samplesToExport, dataOutputHandler, progress);
			syncExportTool.launch();
	
			while (progress.getCurrentStepProgress() < 100 && !progress.isAborted())
				Thread.sleep(500);
		}
		
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
