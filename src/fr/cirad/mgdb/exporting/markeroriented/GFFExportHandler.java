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

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
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
import org.springframework.data.mongodb.core.MongoTemplate;

import com.mongodb.DBCursor;

import fr.cirad.mgdb.exporting.tools.AsyncExportTool;
import fr.cirad.mgdb.exporting.tools.AsyncExportTool.AbstractDataOutputHandler;
import fr.cirad.mgdb.model.mongo.maintypes.GenotypingSample;
import fr.cirad.mgdb.model.mongo.maintypes.VariantData;
import fr.cirad.mgdb.model.mongo.maintypes.VariantRunData;
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
    public void exportData(OutputStream outputStream, String sModule, Collection<GenotypingSample> samples1, Collection<GenotypingSample> samples2, ProgressIndicator progress, DBCursor markerCursor, Map<String, String> markerSynonyms, HashMap<String, Integer> annotationFieldThresholds, HashMap<String, Integer> annotationFieldThresholds2, List<GenotypingSample> samplesToExport, Map<String, InputStream> readyToExportFiles) throws Exception 
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

        int markerCount = markerCursor.count();

        String exportName = sModule + "_" + markerCount + "variants_" + sortedIndividuals.size() + "individuals";
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

		AbstractDataOutputHandler<Integer, LinkedHashMap<VariantData, Collection<VariantRunData>>> dataOutputHandler = new AbstractDataOutputHandler<Integer, LinkedHashMap<VariantData, Collection<VariantRunData>>>() {				
			@Override
			public Void call() {
				for (VariantData variant : variantDataChunkMap.keySet())
					try
					{
						String variantId = variant.getId();
		                List<String> variantDataOrigin = new ArrayList<>();
		                Map<String, List<String>> individualGenotypes = new TreeMap<>(new AlphaNumericComparator<String>());
		                if (variant.getReferencePosition() == null)
		                	unassignedMarkers.add(variantId);
		                
		                if (markerSynonyms != null) {
		                	String syn = markerSynonyms.get(variantId);
		                    if (syn != null)
		                        variantId = syn;
		                }

		                Collection<VariantRunData> runs = variantDataChunkMap.get(variant);
		                if (runs != null) {
		                    for (VariantRunData run : runs) {
								for (Integer sampleId : run.getSampleGenotypes().keySet()) {
									SampleGenotype sampleGenotype = run.getSampleGenotypes().get(sampleId);
		                            String individualId = sampleIdToIndividualMap.get(sampleId);
									if (!VariantData.gtPassesVcfAnnotationFilters(individualId, sampleGenotype, individuals1, annotationFieldThresholds, individuals2, annotationFieldThresholds2))
										continue;	// skip genotype

		                            String gtCode = sampleGenotype.getCode();
		                            List<String> storedIndividualGenotypes = individualGenotypes.get(individualId);
		                            if (storedIndividualGenotypes == null) {
		                                storedIndividualGenotypes = new ArrayList<>();
		                                individualGenotypes.put(individualId, storedIndividualGenotypes);
		                            }
		                            storedIndividualGenotypes.add(gtCode);
		                        }
		                    }
		                }

		                String refAllele = variant.getKnownAlleleList().get(0);
		                String chrom = variant.getReferencePosition() == null ? "0" : variant.getReferencePosition().getSequence();
		                long start = variant.getReferencePosition() == null ? 0 : variant.getReferencePosition().getStartSite();
		                long end = Type.SNP.equals(variant.getType()) ? start : start + refAllele.length() - 1;
		                zos.write((chrom + "\t" + StringUtils.join(variantDataOrigin, ";") /*source*/ + "\t" + typeToOntology.get(variant.getType()) + "\t" + start + "\t" + end + "\t" + "." + "\t" + "+" + "\t" + "." + "\t").getBytes());
		                Comparable syn = markerSynonyms == null ? null : markerSynonyms.get(variant.getId());
		                zos.write(("ID=" + variant.getId() + ";" + (syn != null ? "Name=" + syn + ";" : "") + "alleles=" + StringUtils.join(variant.getKnownAlleleList(), "/") + ";" + "refallele=" + refAllele + ";").getBytes());

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
		                            for (String t : variant.getAllelesFromGenotypeCode(genotype)) {
		                                for (String t1 : variant.getKnownAlleleList()) {
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

		                    List<String> alleles = mostFrequentGenotype == null ? new ArrayList<>() : variant.getAllelesFromGenotypeCode(mostFrequentGenotype);

		                    if (!alleles.isEmpty()) {
		                        zos.write(("acounts=" + individualId + ":").getBytes());

		                        for (String knowAllelesCompt : compt1.keySet()) {
		                            zos.write((knowAllelesCompt + " " + nf.format(compt1.get(knowAllelesCompt) / (float) sum) + " " + compt1.get(knowAllelesCompt) + " ").getBytes());
		                        }
		                        zos.write((alleles.size() + ";").getBytes());
		                    }
		                    if (genotypeCounts.size() > 1) {
		                        Comparable sVariantId = markerSynonyms != null ? markerSynonyms.get(variant.getId()) : variant.getId();
		                        warningFileWriter.write("- Dissimilar genotypes found for variant " + (sVariantId == null ? variant.getId() : sVariantId) + ", individual " + individualId + ". Exporting most frequent: " + StringUtils.join(alleles, ",") + "\n");
		                    }
		                }
		                zos.write((LINE_SEPARATOR).getBytes());		            
		            }
					catch (Exception e)
					{
						if (progress.getError() == null)	// only log this once
							LOG.debug("Unable to export " + variant.getId(), e);
						progress.setError("Unable to export " + variant.getId() + ": " + e.getMessage());
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
                in.readLine();
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
