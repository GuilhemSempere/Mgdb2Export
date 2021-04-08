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
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
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
import fr.cirad.mgdb.exporting.tools.AsyncExportTool;
import fr.cirad.mgdb.exporting.tools.AsyncExportTool.AbstractDataOutputHandler;
import fr.cirad.mgdb.model.mongo.maintypes.GenotypingSample;
import fr.cirad.mgdb.model.mongo.maintypes.VariantData;
import fr.cirad.mgdb.model.mongo.maintypes.VariantDataV2;
import fr.cirad.mgdb.model.mongo.maintypes.VariantRunData;
import fr.cirad.mgdb.model.mongo.maintypes.VariantRunDataV2;
import fr.cirad.mgdb.model.mongo.subtypes.ReferencePosition;
import fr.cirad.mgdb.model.mongo.subtypes.SampleGenotype;
import fr.cirad.mgdb.model.mongodao.MgdbDao;
import fr.cirad.tools.AlphaNumericComparator;
import fr.cirad.tools.Helper;
import fr.cirad.tools.ProgressIndicator;
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
        return "Exports data in GFF3 Format based on Sequence Ontology. See <a target='_blank' href='http://rice.bio.indiana.edu:7082/annot/gff3.html'>http://rice.bio.indiana.edu:7082/annot/gff3.html</a> and <a target='_blank' href='http://www.sequenceontology.org/resources/gff3.html'>http://www.sequenceontology.org/resources/gff3.html</a>";
    }
    
	@Override
	public String getExportArchiveExtension() {
		return "zip";
	}

    /* (non-Javadoc)
	 * @see fr.cirad.mgdb.exporting.markeroriented.AbstractMarkerOrientedExportHandler#exportData(java.io.OutputStream, java.lang.String, java.util.List, fr.cirad.tools.ProgressIndicator, com.mongodb.DBCursor, java.util.Map, int, int, java.util.Map)
     */
    @Override
    public void exportData(OutputStream outputStream, String sModule, Integer nAssemblyId,Collection<GenotypingSample> samples1, Collection<GenotypingSample> samples2, ProgressIndicator progress, MongoCollection<Document> varColl, Document varQuery, Map<String, String> markerSynonyms, HashMap<String, Float> annotationFieldThresholds, HashMap<String, Float> annotationFieldThresholds2, List<GenotypingSample> samplesToExport, Map<String, InputStream> readyToExportFiles) throws Exception 
    {
		List<String> individuals1 = MgdbDao.getIndividualsFromSamples(sModule, samples1).stream().map(ind -> ind.getId()).collect(Collectors.toList());	
		List<String> individuals2 = MgdbDao.getIndividualsFromSamples(sModule, samples2).stream().map(ind -> ind.getId()).collect(Collectors.toList());

		List<String> sortedIndividuals = samplesToExport.stream().map(gs -> gs.getIndividual()).distinct().sorted(new AlphaNumericComparator<String>()).collect(Collectors.toList());
	
        MongoTemplate mongoTemplate = MongoTemplateManager.get(sModule);
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

        File warningFile = File.createTempFile("export_warnings_", "");
        FileWriter warningFileWriter = new FileWriter(warningFile);

		long markerCount = varColl.countDocuments(varQuery);
        String exportName = sModule + "__" + markerCount + "variants__" + sortedIndividuals.size() + "individuals";
        zos.putNextEntry(new ZipEntry(exportName + ".gff3"));
        String header = "##gff-version 3" + LINE_SEPARATOR;
        zos.write(header.getBytes());

        TreeMap<String, String> typeToOntology = new TreeMap<>();
        typeToOntology.put(Type.SNP.toString(), "SO:0000694");
        typeToOntology.put(Type.INDEL.toString(), "SO:1000032");
        typeToOntology.put(Type.MIXED.toString(), "SO:0001059");
        typeToOntology.put(Type.SYMBOLIC.toString(), "SO:0000109");
        typeToOntology.put(Type.MNP.toString(), "SO:0001059");

		final Map<Integer, String> sampleIdToIndividualMap = new HashMap<>();
		for (GenotypingSample gs : samplesToExport)
			sampleIdToIndividualMap.put(gs.getId(), gs.getIndividual());
		
        ArrayList<Comparable> unassignedMarkers = new ArrayList<>();

		AbstractDataOutputHandler<Integer, LinkedHashMap> dataOutputHandler = new AbstractDataOutputHandler<Integer, LinkedHashMap>() {				
			@Override
			public Void call() {
				StringBuffer sb = new StringBuffer();
				for (Object variant : variantDataChunkMap.keySet()) {
					String variantId = null;
					try
					{
						variantId = nAssemblyId == null ? ((VariantDataV2) variant).getId() : ((VariantData) variant).getId();
		                if (markerSynonyms != null) {
		                	String syn = markerSynonyms.get(variantId);
		                    if (syn != null)
		                        variantId = syn;
		                }
		                List<String> variantDataOrigin = new ArrayList<>();
		                Map<String, List<String>> individualGenotypes = new TreeMap<>(new AlphaNumericComparator<String>());
		                if ((nAssemblyId == null && ((VariantDataV2) variant).getReferencePosition() == null) || (nAssemblyId != null && ((VariantData) variant).getReferencePosition(nAssemblyId) == null))
		                	unassignedMarkers.add(variantId);
		                
		                if (markerSynonyms != null) {
		                	String syn = markerSynonyms.get(variantId);
		                    if (syn != null)
		                        variantId = syn;
		                }

		                if (nAssemblyId == null) {
			                Collection<VariantRunDataV2> runs = (Collection<VariantRunDataV2>) variantDataChunkMap.get((VariantDataV2) variant);
			                if (runs != null) {
			                    for (VariantRunDataV2 run : runs) {
			                    	for (Integer sampleId : run.getSampleGenotypes().keySet()) {
										SampleGenotype sampleGenotype = run.getSampleGenotypes().get(sampleId);
										String individualId = sampleIdToIndividualMap.get(sampleId);
			                            
										if (!VariantData.gtPassesVcfAnnotationFiltersV2(individualId, sampleGenotype, individuals1, annotationFieldThresholds, individuals2, annotationFieldThresholds2))
											continue;	// skip genotype
										
			                            String gtCode = sampleGenotype.getCode();
			                            List<String> storedIndividualGenotypes = individualGenotypes.get(individualId);
			                            if (storedIndividualGenotypes == null) {
			                                storedIndividualGenotypes = new ArrayList<String>();
			                                individualGenotypes.put(individualId, storedIndividualGenotypes);
			                            }
			                            storedIndividualGenotypes.add(gtCode);
			                        }
			                    }
			                }
		                }
		                else {
			                Collection<VariantRunData> runs = (Collection<VariantRunData>) variantDataChunkMap.get((VariantData) variant);
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
		                }

		                ReferencePosition rp = nAssemblyId == null ? ((VariantDataV2) variant).getReferencePosition() : ((VariantData) variant).getReferencePosition(nAssemblyId);
		                List<String> knownAlleleList = nAssemblyId == null ? ((VariantDataV2) variant).getKnownAlleleList() : ((VariantData) variant).getKnownAlleleList();
		                String refAllele = knownAlleleList.get(0);
		                String chrom = rp == null ? "0" : rp.getSequence();
		                long start = rp == null ? 0 : rp.getStartSite();
		                String variantType = nAssemblyId == null ? ((VariantDataV2) variant).getType() : ((VariantData) variant).getType();
		                long end = Type.SNP.equals(variantType) ? start : start + refAllele.length() - 1;
		                sb.append(chrom + "\t" + StringUtils.join(variantDataOrigin, ";") /*source*/ + "\t" + typeToOntology.get(variantType) + "\t" + start + "\t" + end + "\t" + "." + "\t" + "+" + "\t" + "." + "\t");
		                Comparable syn = markerSynonyms == null ? null : markerSynonyms.get(variantId);
		                sb.append("ID=" + variantId + ";" + (syn != null ? "Name=" + syn + ";" : "") + "alleles=" + StringUtils.join(knownAlleleList, "/") + ";" + "refallele=" + refAllele + ";");

		                for (String individualId : individualGenotypes.keySet() /* we use this object because it has the proper ordering*/) {
		                    NumberFormat nf = NumberFormat.getInstance(Locale.US);
		                    nf.setMaximumFractionDigits(4);
		                    HashMap<String, Integer> compt1 = new HashMap<>();
		                    int highestGenotypeCount = 0;
		                    int sum = 0;

		                    List<String> genotypes = individualGenotypes.get(individualId);
		                    HashMap<Object, Integer> genotypeCounts = new HashMap<>(); // will help us to keep track of missing genotypes

		                    String mostFrequentGenotype = null;
		                    if (genotypes != null) {
		                        for (String genotype : genotypes) {
		                            if (genotype == null) {
		                                continue; /* skip missing genotypes */
		                            }

		                            int count = 0;
		                            List<String> alleles = nAssemblyId == null ? ((VariantDataV2) variant).getAllelesFromGenotypeCode(genotype) : ((VariantData) variant).getAllelesFromGenotypeCode(genotype);
		                            for (String t : alleles) {
		                                for (String t1 : knownAlleleList) {
		                                    if (t.equals(t1) && !(compt1.containsKey(t1))) {
		                                        count++;
		                                        compt1.put(t1, count);
		                                    } else if (t.equals(t1) && compt1.containsKey(t1)) {
		                                        if (compt1.get(t1) != 0) {
		                                            count++;
		                                            compt1.put(t1, count);
		                                        } else {
		                                            compt1.put(t1, count);
		                                        }
		                                    } else if (!(compt1.containsKey(t1))) {
		                                        compt1.put(t1, 0);
		                                    }
		                                }
		                            }
		                            for (int countValue : compt1.values()) {
		                                sum += countValue;
		                            }

		                            int gtCount = 1 + Helper.getCountForKey(genotypeCounts, genotype);
		                            if (gtCount > highestGenotypeCount) {
		                                highestGenotypeCount = gtCount;
		                                mostFrequentGenotype = genotype;
		                            }
		                            genotypeCounts.put(genotype, gtCount);
		                        }
		                    }

		                    List<String> alleles = mostFrequentGenotype == null ? new ArrayList<>() : (nAssemblyId == null ? ((VariantDataV2) variant).getAllelesFromGenotypeCode(mostFrequentGenotype) : ((VariantData) variant).getAllelesFromGenotypeCode(mostFrequentGenotype));

		                    if (!alleles.isEmpty()) {
		                        sb.append("acounts=" + individualId + ":");

		                        for (String knowAllelesCompt : compt1.keySet()) {
		                            sb.append(knowAllelesCompt + " " + nf.format(compt1.get(knowAllelesCompt) / (float) sum) + " " + compt1.get(knowAllelesCompt) + " ");
		                        }
		                        sb.append(alleles.size() + ";");
		                    }
		                    if (genotypeCounts.size() > 1) {
		                        Comparable sVariantId = markerSynonyms != null ? markerSynonyms.get(variantId) : variantId;
		                        warningFileWriter.write("- Dissimilar genotypes found for variant " + (sVariantId == null ? variantId : sVariantId) + ", individual " + individualId + ". Exporting most frequent: " + StringUtils.join(alleles, ",") + "\n");
		                    }
		                }
		                sb.append(LINE_SEPARATOR);		            
		            }
					catch (Exception e)
					{
						if (progress.getError() == null)	// only log this once
							LOG.debug("Unable to export " + variantId, e);
						progress.setError("Unable to export " + variantId + ": " + e.getMessage());
					}
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
		try (MongoCursor<Document> markerCursor = IExportHandler.getMarkerCursorWithCorrectCollation(varColl, varQuery, nAssemblyId, nQueryChunkSize)) {
			AsyncExportTool asyncExportTool = new AsyncExportTool(markerCursor, markerCount, nQueryChunkSize, mongoTemplate, samplesToExport, dataOutputHandler, progress);
			asyncExportTool.launch();
	
			while (progress.getCurrentStepProgress() < 100 && !progress.isAborted())
				Thread.sleep(500);
		}		
        zos.closeEntry();
        
        if (unassignedMarkers.size() > 0)
        	LOG.info("No chromosomal position found for " + unassignedMarkers.size() + " markers " + StringUtils.join(unassignedMarkers, ", "));

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
        return Arrays.asList(new String[]{"Exporting data to GFF3 format"});
    }
    
	@Override
	public String[] getExportDataFileExtensions() {
		return new String[] {"gff3"};
	}
}
