/*******************************************************************************
 * MGDB - Mongo Genotype DataBase
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

import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.PipedInputStream;
import java.io.PipedOutputStream;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Modifier;
import java.nio.charset.StandardCharsets;
import java.util.Collection;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.concurrent.atomic.AtomicReference;
import java.util.stream.Collectors;

import org.apache.log4j.Logger;
import org.bson.Document;
import org.springframework.beans.factory.config.BeanDefinition;
import org.springframework.context.annotation.ClassPathScanningCandidateComponentProvider;
import org.springframework.core.type.filter.AssignableTypeFilter;
import org.springframework.data.mongodb.core.query.Criteria;
import org.springframework.data.mongodb.core.query.Query;

import fr.cirad.mgdb.exporting.IExportHandler;
import fr.cirad.mgdb.exporting.markeroriented.AbstractMarkerOrientedExportHandler;
import fr.cirad.mgdb.exporting.markeroriented.HapMapExportHandler;
import fr.cirad.mgdb.exporting.tools.ExportManager.ExportOutputs;
import fr.cirad.mgdb.model.mongo.maintypes.Assembly;
import fr.cirad.mgdb.model.mongo.maintypes.CallSet;
import fr.cirad.mgdb.model.mongo.maintypes.GenotypingSample;
import fr.cirad.tools.Helper;
import fr.cirad.tools.ProgressIndicator;
import fr.cirad.tools.mgdb.VariantQueryWrapper;
import fr.cirad.tools.mongo.MongoTemplateManager;

/**
 * The Class AbstractIndividualOrientedExportHandler.
 */
public abstract class AbstractIndividualOrientedExportHandler implements IExportHandler
{
	
	/** The Constant LOG. */
	private static final Logger LOG = Logger.getLogger(AbstractIndividualOrientedExportHandler.class);
	
	/** The individual oriented export handlers. */
	static private TreeMap<String, AbstractIndividualOrientedExportHandler> individualOrientedExportHandlers = null;

	protected String tmpFolderPath;
		
	/**
	 * Export data.
	 *
	 * @param outputStream the output stream
	 * @param sModule the module
     * @param nAssemblyId ID of the assembly to work with
	 * @param exportOutputs contains individual export files, variant files (if any), warning files, all chunked
	 * @param fDeleteSampleExportFilesOnExit whether or not to delete sample export files on exit
	 * @param progress the progress
	 * @param tmpVarCollName the variant collection name (null if not temporary)
	 * @param varQueryWrapper variant query wrapper
	 * @param markerCount number of variants to export
	 * @param markerSynonyms the marker synonyms
	 * @param metadataPopField metadata field to use as population String (overriding "fixed" individual-population field if exists) 
	 * @param readyToExportFiles the ready to export files
	 * @throws Exception the exception
	 */
	
	abstract public void exportData(OutputStream outputStream, String sModule, Integer nAssemblyId, ExportOutputs exportOutputs, boolean fDeleteSampleExportFilesOnExit, ProgressIndicator progress, String tmpVarCollName, VariantQueryWrapper varQueryWrapper, long markerCount, Map<String, String> markerSynonyms, Map<String, String> individualPopulations, Map<String, InputStream> readyToExportFiles) throws Exception;
	/**
	 * Creates the export files.
	 *
	 * @param sModule the module
     * @param nAssemblyId ID of the assembly to work with
	 * @param sExportingUser the user who launched the export 
	 * @param tmpVarCollName the variant collection name (null if not temporary)
	 * @param variantQuery the query to apply to the variant collection
	 * @param markerCount number of variants to export
	 * @param exportID the export id
	 * @param individualsByPop List of the individuals in each group
	 * @param workWithSamples if true export sample names, otherwise individual names
	 * @param annotationFieldThresholds the annotation field thresholds for each group
	 * @param callSetsToExport
	 * @param metadataFieldsToExport metadata fields to export for individuals
	 * @param progress the progress
	 * @return a map providing one File per individual
	 * @throws Exception the exception
	 */
	public ExportOutputs createExportFiles(String sModule, Integer nAssemblyId, String sExportingUser, String tmpVarCollName, Document variantQuery, long markerCount, String exportID, Map<String, Collection<String>> individualsByPop, boolean workWithSamples, Map<String, HashMap<String, Float>> annotationFieldThresholds, List<CallSet> callSetsToExport, Collection<String> metadataFieldsToExport, final ProgressIndicator progress) throws Exception
	{
		long before = System.currentTimeMillis();

		Map<String, Integer> individualPositions = IExportHandler.buildIndividualPositions(callSetsToExport, workWithSamples);
		
		File[] indFiles = new File[individualPositions.size()];
		BufferedOutputStream[] indOS = new BufferedOutputStream[individualPositions.size()];
		int i = 0;
		for (String individual : individualPositions.keySet()) {
			indFiles[i] = File.createTempFile(exportID.replaceAll("\\|", "&curren;") +  "-" + individual + "-", ".tsv");
			if (i == 0)
				LOG.debug("First individual file for export " + exportID + ": " + indFiles[i].getPath());
			indOS[i] = new BufferedOutputStream(new FileOutputStream(indFiles[i]));
			indOS[i++].write((individual + LINE_SEPARATOR).getBytes());
		}

        final Map<Integer, String> callSetIdToIndividualMap = new HashMap<>();
        for (CallSet cs : callSetsToExport)
        	callSetIdToIndividualMap.put(cs.getId(), workWithSamples ? cs.getSampleId() : cs.getIndividual());
		PipedOutputStream pos = new PipedOutputStream();
		AtomicReference<ExportOutputs> exportOutputs = new AtomicReference<>();
		
		Integer nAssemblyID= Assembly.getThreadBoundAssembly();
//        boolean workWithSamples = callSetsToExport.stream().filter(sp -> sp.isDetached()).count() == callSetsToExport.size();

		// Run data reading in a thread so that we can immediately wait for the streamed output (avoids the need for a global temporary file)
		Thread hapMapExportThread = new Thread(() -> {
		    try {
				Assembly.setThreadAssembly(nAssemblyID);
		        progress.addStep("Extracting genotypes");
		        progress.moveToNextStep();

		        HapMapExportHandler heh = (HapMapExportHandler) AbstractMarkerOrientedExportHandler.getMarkerOrientedExportHandlers().get("HAPMAP");
		        exportOutputs.set(heh.writeGenotypeFile(true, false, true, true, pos, sModule, MongoTemplateManager.get(sModule).findOne(new Query(Criteria.where("_id").is(nAssemblyID)), Assembly.class), individualsByPop, workWithSamples, callSetIdToIndividualMap, annotationFieldThresholds, progress, tmpVarCollName, variantQuery, markerCount, null, callSetsToExport));
		        exportOutputs.get().setWorkWithSamples(workWithSamples);
		    } catch (Exception e) {
		        LOG.error("Error reading genotypes for export", e);
		        progress.setError(e.getMessage());
		    } finally {
		        Assembly.cleanupThreadAssembly();
		    }
		});
		hapMapExportThread.start();

		if (progress.isAborted() || progress.getError() != null)
			return exportOutputs.get();
		
        int nLinesProcessed = -1;
    	Integer[] indPos = individualPositions.values().toArray(new Integer[individualPositions.size()]);
        String line;
        try (BufferedReader reader = new BufferedReader(new InputStreamReader(new PipedInputStream(pos), StandardCharsets.UTF_8))) {

	        while (nLinesProcessed < markerCount && (line = reader.readLine()) != null) {
	        	if (nLinesProcessed++ == -1)
	        		continue;	// skip header

	        	String[] splitLine = line.split("\t", -1);

	        	// write genotypes collected in this chunk to each individual's file
	        	for (int j=4; j<splitLine.length; j++) {
					indOS[indPos[j - 4]].write((splitLine[j].replace("/", " ")).getBytes(StandardCharsets.UTF_8));
					indOS[indPos[j - 4]].write(LINE_SEPARATOR.getBytes(StandardCharsets.UTF_8));
	        	}
	        }
	        for (BufferedOutputStream ios : indOS)
	        	ios.close();
	        
    		hapMapExportThread.join();
	        exportOutputs.get().setGenotypeFiles(indFiles);
        }
        catch (Exception e) {
        	if (!(e instanceof IOException || e instanceof FileNotFoundException) || !progress.isAborted())
        		throw e;
        }

        if (metadataFieldsToExport != null && !metadataFieldsToExport.isEmpty()) {
            String exportName = IExportHandler.buildExportName(sModule, MongoTemplateManager.get(sModule).findOne(new Query(Criteria.where("_id").is(nAssemblyId)), Assembly.class), markerCount, individualPositions.size(), workWithSamples);        
        	String initialContents = getMetadataContentsPrefix() != null ? getMetadataContentsPrefix() : (workWithSamples ? "sample" : "individual");
        	String metadataContents = IExportHandler.buildMetadataFile(sModule, sExportingUser, individualPositions.keySet(), metadataFieldsToExport, workWithSamples, initialContents);
            if (metadataContents.length() != initialContents.length())
            	exportOutputs.get().setMetadata(exportName + "." + getMetadataFileExtension(), metadataContents);
        }

	 	if (!progress.isAborted())
	 		LOG.info("createExportFiles took " + (System.currentTimeMillis() - before)/1000d + "s to process " + markerCount + " variants and " + indFiles.length + " individuals");
		
		return exportOutputs.get();
	}
	
	/**
	 * Gets the individual oriented export handlers.
	 *
	 * @return the individual oriented export handlers
	 * @throws ClassNotFoundException the class not found exception
	 * @throws InstantiationException the instantiation exception
	 * @throws IllegalAccessException the illegal access exception
	 * @throws IllegalArgumentException the illegal argument exception
	 * @throws InvocationTargetException the invocation target exception
	 * @throws NoSuchMethodException the no such method exception
	 * @throws SecurityException the security exception
	 */
	public static TreeMap<String, AbstractIndividualOrientedExportHandler> getIndividualOrientedExportHandlers() throws ClassNotFoundException, InstantiationException, IllegalAccessException, IllegalArgumentException, InvocationTargetException, NoSuchMethodException, SecurityException
	{
		if (individualOrientedExportHandlers == null)
		{
			individualOrientedExportHandlers = new TreeMap<String, AbstractIndividualOrientedExportHandler>();
			ClassPathScanningCandidateComponentProvider provider = new ClassPathScanningCandidateComponentProvider(false);
			provider.addIncludeFilter(new AssignableTypeFilter(AbstractIndividualOrientedExportHandler.class));
			try
			{
				for (BeanDefinition component : provider.findCandidateComponents("fr.cirad"))
				{
				    Class cls = Class.forName(component.getBeanClassName());
				    if (!Modifier.isAbstract(cls.getModifiers()))
				    {
						AbstractIndividualOrientedExportHandler exportHandler = (AbstractIndividualOrientedExportHandler) cls.getConstructor().newInstance();
						String sFormat = exportHandler.getExportFormatName();
						AbstractIndividualOrientedExportHandler previouslyFoundExportHandler = individualOrientedExportHandlers.get(sFormat);
						if (previouslyFoundExportHandler != null)
						{
							if (exportHandler.getClass().isAssignableFrom(previouslyFoundExportHandler.getClass()))
							{
								LOG.debug(previouslyFoundExportHandler.getClass().getName() + " implementation was preferred to " + exportHandler.getClass().getName() + " to handle exporting to '" + sFormat + " format");
								continue;	// skip adding the current exportHandler because we already have a "better" one
							}
							else if (previouslyFoundExportHandler.getClass().isAssignableFrom(exportHandler.getClass()))
								LOG.debug(exportHandler.getClass().getName() + " implementation was preferred to " + previouslyFoundExportHandler.getClass().getName() + " to handle exporting to " + sFormat + " format");
							else
								LOG.warn("Unable to choose between " + previouslyFoundExportHandler.getClass().getName() + " and " + exportHandler.getClass().getName() + ". Keeping first found: " + previouslyFoundExportHandler.getClass().getName());
						}
				    	individualOrientedExportHandlers.put(sFormat, exportHandler);
				    }
				}
			}
			catch (Exception e)
			{
				LOG.warn("Error scanning export handlers", e);
			}
		}
		return individualOrientedExportHandlers;
	}
	
	protected String findOutMostFrequentGenotype(String line, OutputStream warningOS, int nMarkerIndex, String individualId) throws IOException {
        String mostFrequentGenotype = null;
        if (!line.isEmpty()) {
            List<String> genotypes = Helper.split(line, "|");
            if (genotypes.size() == 1)
                mostFrequentGenotype = genotypes.get(0);
            else {
                HashMap<Object, Integer> genotypeCounts = new HashMap<Object, Integer>();   // will help us to keep track of missing genotypes
                int highestGenotypeCount = 0;

                for (String genotype : genotypes) {
                    if (genotype == null)
                        continue;   /* skip missing genotypes */

                    int gtCount = 1 + Helper.getCountForKey(genotypeCounts, genotype);
                    if (gtCount > highestGenotypeCount) {
                        highestGenotypeCount = gtCount;
                        mostFrequentGenotype = genotype;
                    }
                    genotypeCounts.put(genotype, gtCount);
                }

                if (genotypeCounts.size() > 1) {
                    List<Integer> reverseSortedGtCounts = genotypeCounts.values().stream().sorted(Comparator.reverseOrder()).collect(Collectors.toList());
                    if (reverseSortedGtCounts.get(0) == reverseSortedGtCounts.get(1))
                        mostFrequentGenotype = null;
                    if (warningOS != null)
                    	warningOS.write(("- Dissimilar genotypes found for variant n. " + nMarkerIndex + ", individual " + individualId + ". " + (mostFrequentGenotype == null ? "Exporting as missing data" : "Exporting most frequent: " + mostFrequentGenotype.replaceAll(" ", "/")) + "\n").getBytes());
                }
            }
        }
        return mostFrequentGenotype;
	}

	/* (non-Javadoc)
	 * @see fr.cirad.mgdb.exporting.IExportHandler#getSupportedVariantTypes()
	 */
	@Override
	public List<String> getSupportedVariantTypes()
	{
		return null;	// means any type
	}
	
	@Override
	public String getExportContentType() {
		return "application/zip";
	}
	
	@Override
	public void setTmpFolder(String tmpFolderPath) {
		this.tmpFolderPath = tmpFolderPath;	
	}

	@Override
	public String getMetadataFileExtension() {
		return "tsv";
	}

	@Override
	public String getMetadataContentsPrefix() {
		return null;
	}
}