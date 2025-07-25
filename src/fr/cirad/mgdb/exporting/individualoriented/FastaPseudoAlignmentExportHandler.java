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

import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.Scanner;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;

import org.apache.log4j.Logger;
import org.bson.Document;
import org.springframework.data.mongodb.core.MongoTemplate;
import org.springframework.data.mongodb.core.query.Criteria;
import org.springframework.data.mongodb.core.query.Query;

import com.mongodb.BasicDBList;
import com.mongodb.client.MongoCursor;

import fr.cirad.mgdb.exporting.IExportHandler;
import fr.cirad.mgdb.exporting.tools.ExportManager.ExportOutputs;
import fr.cirad.mgdb.model.mongo.maintypes.Assembly;
import fr.cirad.mgdb.model.mongo.maintypes.VariantData;
import fr.cirad.mgdb.model.mongo.subtypes.ReferencePosition;
import fr.cirad.tools.Helper;
import fr.cirad.tools.ProgressIndicator;
import fr.cirad.tools.mgdb.VariantQueryWrapper;
import fr.cirad.tools.mongo.MongoTemplateManager;
import htsjdk.variant.variantcontext.VariantContext.Type;

/**
 * The Class FastaAlignmentExportHandler.
 */
public class FastaPseudoAlignmentExportHandler extends AbstractIndividualOrientedExportHandler {

    /**
     * The Constant LOG.
     */
    private static final Logger LOG = Logger.getLogger(FastaPseudoAlignmentExportHandler.class);

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
        return "FASTA";
    }

    /* (non-Javadoc)
	 * @see fr.cirad.mgdb.exporting.IExportHandler#getExportFormatDescription()
     */
    @Override
    public String getExportFormatDescription() {
    	return "Exports a zipped FASTA file containing a pseudo-alignment consisting in the concatenation of SNP alleles, compatible with tree construction tools like FastTree, ClearCut or RapidNJ. An additional PLINK-style map file is added for reference.";
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
    public void exportData(OutputStream outputStream, String sModule, Integer nAssemblyId, String sExportingUser, ExportOutputs exportOutputs, boolean fDeleteSampleExportFilesOnExit, ProgressIndicator progress, String tmpVarCollName, VariantQueryWrapper varQueryWrapper, long markerCount, Map<String, String> markerSynonyms, Collection<String> individualMetadataFieldsToExport, Map<String, String> individualPopulations, Map<String, InputStream> readyToExportFiles) throws Exception {
		MongoTemplate mongoTemplate = MongoTemplateManager.get(sModule);
        int nQueryChunkSize = IExportHandler.computeQueryChunkSize(mongoTemplate, markerCount);

        // save existing warnings into a temp file so we can append to it
        File warningFile = File.createTempFile("export_warnings_", "");
        try {
	        FileOutputStream warningOS = new FileOutputStream(warningFile);
	        for (File f : exportOutputs.getWarningFiles()) {
		    	if (f != null && f.length() > 0) {
		            BufferedReader in = new BufferedReader(new FileReader(f));
		            String sLine;
		            while ((sLine = in.readLine()) != null)
		            	warningOS.write((sLine + "\n").getBytes());
		            in.close();
			    	f.delete();
		    	}
	        }
	
	        ZipOutputStream zos = IExportHandler.createArchiveOutputStream(outputStream, readyToExportFiles);
			Assembly assembly = mongoTemplate.findOne(new Query(Criteria.where("_id").is(nAssemblyId)), Assembly.class);
			String exportName = sModule + (assembly != null && assembly.getName() != null ? "__" + assembly.getName() : "") + "__" + markerCount + "variants__" + exportOutputs.getGenotypeFiles().length + "individuals";
	        
	        ArrayList<String> exportedIndividuals = new ArrayList<>();
	        for (File indFile : exportOutputs.getGenotypeFiles())
	        	try (Scanner scanner = new Scanner(indFile)) {
	        		exportedIndividuals.add(scanner.nextLine());
	        	}
	
	        if (individualMetadataFieldsToExport == null || !individualMetadataFieldsToExport.isEmpty())
	        	IExportHandler.addMetadataEntryIfAny(sModule + "__" + exportOutputs.getGenotypeFiles().length + "individuals_metadata.tsv", sModule, sExportingUser, exportedIndividuals, individualMetadataFieldsToExport, zos, "individual");
	        
	        Collection<BasicDBList> variantDataQueries = varQueryWrapper.getVariantDataQueries();
	        Document variantQueryForTargetCollection = variantDataQueries.isEmpty() ? new Document() : (new Document("$and", tmpVarCollName == null ? variantDataQueries.iterator().next() : (varQueryWrapper.getBareQueries().iterator().hasNext() ? varQueryWrapper.getBareQueries().iterator().next() : new BasicDBList())));
	
	        zos.putNextEntry(new ZipEntry(exportName + "." + getExportDataFileExtensions()[0]));
	        zos.write(getHeaderlines(exportOutputs.getGenotypeFiles().length, (int) markerCount).getBytes());
	        writeGenotypeFile(zos, sModule, exportedIndividuals, nQueryChunkSize, markerSynonyms, exportOutputs.getGenotypeFiles(), warningOS, progress);
	        zos.write(getFooterlines().getBytes());
	    	zos.closeEntry();
	    	
	        zos.putNextEntry(new ZipEntry(exportName + ".map"));
	        String refPosPath = Assembly.getVariantRefPosPath(nAssemblyId);
	        int nMarkerIndex = 0;
	        ArrayList<Comparable> unassignedMarkers = new ArrayList<>();
	    	String refPosPathWithTrailingDot = Assembly.getThreadBoundVariantRefPosPath() + ".";
	    	Document projectionAndSortDoc = new Document(refPosPathWithTrailingDot + ReferencePosition.FIELDNAME_SEQUENCE, 1).append(refPosPathWithTrailingDot + ReferencePosition.FIELDNAME_START_SITE, 1);
	    	try (MongoCursor<Document> markerCursor = IExportHandler.getMarkerCursorWithCorrectCollation(mongoTemplate.getCollection(tmpVarCollName != null ? tmpVarCollName : mongoTemplate.getCollectionName(VariantData.class)), variantQueryForTargetCollection, projectionAndSortDoc, nQueryChunkSize)) {
	            progress.addStep("Writing map file");
	            progress.moveToNextStep();
		        while (markerCursor.hasNext()) {
		            Document exportVariant = markerCursor.next();
		            Document refPos = (Document) Helper.readPossiblyNestedField(exportVariant, refPosPath, ";", null);
		            Long pos = refPos == null ? null : ((Number) refPos.get(ReferencePosition.FIELDNAME_START_SITE)).longValue();
		            String chrom = refPos == null ? null : (String) refPos.get(ReferencePosition.FIELDNAME_SEQUENCE);
	                String markerId = (String) exportVariant.get("_id");
		            if (chrom == null)
		            	unassignedMarkers.add(markerId);
		            String exportedId = markerSynonyms == null ? markerId : markerSynonyms.get(markerId);
		            zos.write(((chrom == null ? "0" : chrom) + " " + exportedId + " " + 0 + " " + (pos == null ? 0 : pos) + LINE_SEPARATOR).getBytes());
	
	                progress.setCurrentStepProgress(nMarkerIndex++ * 100 / markerCount);
		        }
			}
	        zos.closeEntry();
	
	        warningOS.close();
	        if (warningFile.length() > 0) {
	            progress.addStep("Adding lines to warning file");
	            progress.moveToNextStep();
	            progress.setPercentageEnabled(false);
	            zos.putNextEntry(new ZipEntry(exportName + "-REMARKS.txt"));
	            int nWarningCount = 0;
	            BufferedReader in = new BufferedReader(new FileReader(warningFile));
	            String sLine;
	            while ((sLine = in.readLine()) != null) {
	                zos.write((sLine + "\n").getBytes());
	                progress.setCurrentStepProgress(nWarningCount++);
	            }
	            LOG.info("Number of Warnings for export (" + exportName + "): " + nWarningCount);
	            in.close();
	            zos.closeEntry();
	        }
	        zos.finish();
	        zos.close();
	    }
        finally {
            warningFile.delete();
        }
        progress.setPercentageEnabled(true);
        progress.setCurrentStepProgress((short) 100);
    }

	public void writeGenotypeFile(OutputStream os, String sModule, Collection<String> individualsToExport, int nQueryChunkSize, Map<String, String> markerSynonyms, File[] individualExportFiles, OutputStream warningOS, ProgressIndicator progress) throws IOException, InterruptedException {
        short nProgress = 0, nPreviousProgress = 0;
        
        int i = 0, nNConcurrentThreads = Math.max(1, Runtime.getRuntime().availableProcessors());	// use multiple threads so we can prepare several lines at once
        StringBuilder[] individualLines = new StringBuilder[nNConcurrentThreads];
        final ArrayList<Thread> threadsToWaitFor = new ArrayList<>(nNConcurrentThreads);
        final AtomicInteger initialStringBuilderCapacity = new AtomicInteger();

        try
        {
            int nWrittenIndividualCount = 0;
	        for (final File f : individualExportFiles) {
	            if (progress.isAborted() || progress.getError() != null)
	            	return;

	        	final int nThreadIndex = i % nNConcurrentThreads;	        	
	            Thread thread = new Thread() {
	                @Override
	                public void run() {
                        StringBuilder indLine = individualLines[nThreadIndex];
                        if (indLine == null) {
                            indLine = new StringBuilder((int) f.length() / 3 /* rough estimation */);
                            individualLines[nThreadIndex] = indLine;
                        }

	                	BufferedReader in = null;
	    	            try {
		    	            in = new BufferedReader(new FileReader(f));
			                String individualId, line = in.readLine();	// read sample id
			                if (line != null) {
			                    individualId = line;
			                    indLine.append(getLinePrefix()).append(individualId).append(getIndividualToSequenceSeparator());
			                } else {
			                    throw new Exception("Unable to read first line of temp export file " + f.getName());
			                }
			
			                int nMarkerIndex = 0;
			                while ((line = in.readLine()) != null) {
                                String mostFrequentGenotype = findOutMostFrequentGenotype(line, warningOS, nMarkerIndex, individualId);
			                    String[] alleles = mostFrequentGenotype == null ? new String[0] : mostFrequentGenotype.replaceAll("\\*", "-").split(" ");
			                    if (alleles.length > 2 && warningOS != null)
			                    	warningOS.write(("- More than 2 alleles found for variant n. " + nMarkerIndex + ", individual " + individualId + ". Exporting only the first 2 alleles.\n").getBytes());
			
			                    String all1 = alleles.length == 0 ? getMissingAlleleString() : alleles[0];
			                    String all2 = alleles.length == 0 ? getMissingAlleleString() : alleles[alleles.length == 1 ? 0 : 1];
		                        indLine.append(all1).append(all2);
			
			                    nMarkerIndex++;
			                }
	    	            } catch (Exception e) {
	    	                LOG.error("Error exporting data", e);
	    	                progress.setError("Error exporting data: " + e.getClass().getSimpleName() + (e.getMessage() != null ? " - " + e.getMessage() : ""));
	    	            } finally {
	    	            	if (in != null)
								try {
									in.close();
								} catch (IOException e) {
									LOG.warn(e);
								}               
	    	            }

	    	            indLine.append(LINE_SEPARATOR);
	                }
	            };
	            
	            thread.start();
	            threadsToWaitFor.add(thread);

	            if (++i % nNConcurrentThreads == 0 || i == individualExportFiles.length) {
	                for (Thread t : threadsToWaitFor) // wait for all previously launched async threads
	               		t.join();
	                
                    for (int j=0; j<nNConcurrentThreads && nWrittenIndividualCount++ < individualExportFiles.length; j++) {
                        StringBuilder indLine = individualLines[j];
                        if (indLine == null || indLine.length() == 0)
                            LOG.warn("No line to export for individual " + j);
                        else {
                            os.write(indLine.toString().getBytes());
                            individualLines[j] = new StringBuilder(initialStringBuilderCapacity.get());
                        }
                    }

    	            nProgress = (short) (i * 100 / individualExportFiles.length);
    	            if (nProgress > nPreviousProgress) {
    	                progress.setCurrentStepProgress(nProgress);
    	                nPreviousProgress = nProgress;
    	            }
	                threadsToWaitFor.clear();
	            }
	        }
        }
        finally
        {
        	for (File f : individualExportFiles)
                if (!f.delete()) {
                    f.deleteOnExit();
                    LOG.info("Unable to delete tmp export file " + f.getAbsolutePath());
                }
        }
    }

    protected String getMissingAlleleString() {
		return "-";
	}
    
    protected String getLinePrefix() {
		return ">";
	}
    
    protected String getIndividualToSequenceSeparator() {
		return "\n";
	}
    
    protected String getHeaderlines(int nIndividualCount, int nMarkerCount) {
		return "";
	}

	protected String getFooterlines() {
		return "";
	}
	
	/* (non-Javadoc)
	 * @see fr.cirad.mgdb.exporting.IExportHandler#getStepList()
     */
    @Override
    public List<String> getStepList() {
        return Arrays.asList(new String[]{"Exporting data to FASTA format"});
    }

	@Override
	public String[] getExportDataFileExtensions() {
		return new String[] {"fasta"};
	}

    @Override
    public int[] getSupportedPloidyLevels() {
        return new int[] {2};
    }
}