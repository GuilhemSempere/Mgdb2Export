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
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.stream.Collectors;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;

import fr.cirad.mgdb.model.mongo.maintypes.*;
import org.apache.log4j.Logger;
import org.bson.Document;
import org.springframework.data.mongodb.core.MongoTemplate;
import org.springframework.data.mongodb.core.query.Criteria;
import org.springframework.data.mongodb.core.query.Query;

import com.mongodb.BasicDBList;
import com.mongodb.BasicDBObject;
import com.mongodb.client.MongoCollection;
import com.mongodb.client.MongoCursor;

import fr.cirad.mgdb.exporting.IExportHandler;
import fr.cirad.mgdb.exporting.tools.ExportManager;
import fr.cirad.mgdb.model.mongo.maintypes.DBVCFHeader.VcfHeaderId;
import fr.cirad.mgdb.model.mongo.subtypes.AbstractVariantData;
import fr.cirad.mgdb.model.mongo.subtypes.Callset;
import fr.cirad.mgdb.model.mongo.subtypes.ReferencePosition;
import fr.cirad.mgdb.model.mongodao.MgdbDao;
import fr.cirad.tools.AlphaNumericComparator;
import fr.cirad.tools.Helper;
import fr.cirad.tools.ProgressIndicator;
import fr.cirad.tools.mgdb.VariantQueryWrapper;
import fr.cirad.tools.mongo.MongoTemplateManager;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.CustomVCFWriter;
import htsjdk.variant.vcf.VCFContigHeaderLine;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;

/**
 * The Class VcfExportHandler.
 */
public class VcfExportHandler extends AbstractMarkerOrientedExportHandler {

	/** The Constant LOG. */
	private static final Logger LOG = Logger.getLogger(VcfExportHandler.class);

	/* (non-Javadoc)
	 * @see fr.cirad.mgdb.exporting.IExportHandler#getExportFormatName()
	 */
	@Override
	public String getExportFormatName()
	{
		return "VCF";
	}

	/* (non-Javadoc)
	 * @see fr.cirad.mgdb.exporting.IExportHandler#getExportFormatDescription()
	 */
	@Override
	public String getExportFormatDescription()
	{
		return "Exports data in Variant Call Format. See <a target='_blank' href='http://samtools.github.io/hts-specs/VCFv4.1.pdf'>http://samtools.github.io/hts-specs/VCFv4.1.pdf</a> for more details.";
	}

	/* (non-Javadoc)
	 * @see fr.cirad.mgdb.exporting.IExportHandler#getStepList()
	 */
	@Override
	public List<String> getStepList()
	{
		return Arrays.asList(new String[] {"Creating sequence list", "Exporting data to VCF format"});
	}
	
	@Override
	public String getExportArchiveExtension() {
		return "zip";
	}

	/**
	 * Creates the sam sequence dictionary.
	 *
	 * @param sModule the module
	 * @param sequenceIDs the sequence IDs
	 * @return the SAM sequence dictionary
	 * @throws Exception the exception
	 */
	public SAMSequenceDictionary createSAMSequenceDictionary(String sModule, Collection<String> sequenceIDs) throws Exception
	{
		SAMSequenceDictionary dict1 = new SAMSequenceDictionary();
		MongoTemplate mongoTemplate = MongoTemplateManager.get(sModule);
		if (mongoTemplate.collectionExists(MongoTemplateManager.getMongoCollectionName(Sequence.class))) {
			long before = System.currentTimeMillis();
			Query query = sequenceIDs.isEmpty() ? new Query() : new Query(Criteria.where("_id").in(sequenceIDs));
			query.fields().include(Sequence.FIELDNAME_LENGTH).include(Sequence.FIELDNAME_ASSEMBLY);
			for (Sequence seq : MongoTemplateManager.get(sModule).find(query, Sequence.class))
				dict1.addSequence(new SAMSequenceRecord((String) seq.getId(), (int) seq.getLength()) {{setAssembly(seq.getAssembly());}} );
	    	LOG.info("createSAMSequenceDictionary took " + (System.currentTimeMillis() - before)/1000d + "s to write " + sequenceIDs.size() + " sequences");
		}
		else
			LOG.info("No sequence data was found. No SAMSequenceDictionary could be created.");
	    return dict1;
	}
	
	protected CustomVCFWriter buildVariantContextWriter(OutputStream os, SAMSequenceDictionary dict) {
		return new CustomVCFWriter(null, os, dict, false, false, true);
	}

    @Override
    public void exportData(OutputStream outputStream, String sModule, Integer nAssemblyId, String sExportingUser, ProgressIndicator progress, String tmpVarCollName, VariantQueryWrapper varQueryWrapper, long markerCount, Map<String, String> markerSynonyms, Map<String, Collection<String>> individualsByPop, boolean workWithSamples, Map<String, HashMap<String, Float>> annotationFieldThresholds, Collection<Callset> callSetsToExport, Collection<String> individualMetadataFieldsToExport, Map<String, InputStream> readyToExportFiles) throws Exception {
    	TreeSet<String> indSet = new TreeSet<>(new AlphaNumericComparator<String>());
		for (Callset cs : callSetsToExport)
			indSet.add(workWithSamples ? cs.getSampleId() : cs.getIndividual());
		List<String> sortedIndividuals = new ArrayList<>(indSet);
		
		MongoTemplate mongoTemplate = MongoTemplateManager.get(sModule);
        ZipOutputStream zos = IExportHandler.createArchiveOutputStream(outputStream, readyToExportFiles, null);
        
        Assembly assembly = mongoTemplate.findOne(new Query(Criteria.where("_id").is(nAssemblyId)), Assembly.class);
//        boolean workWithSamples = callSetsToExport.stream().filter(sp -> sp.isDetached()).count() == callSetsToExport.size();
        String exportName = IExportHandler.buildExportName(sModule, assembly, markerCount, indSet.size(), workWithSamples);
        
        if (individualMetadataFieldsToExport == null || !individualMetadataFieldsToExport.isEmpty())
        	IExportHandler.addMetadataEntryIfAny(sModule + "__" + indSet.size() + (workWithSamples ? "sample" : "individual" ) + "s_metadata.tsv", sModule, sExportingUser, indSet, individualMetadataFieldsToExport, zos, (workWithSamples ? "sample" : "individual"), workWithSamples);
        
        String refPosPathWithTrailingDot = Assembly.getVariantRefPosPath(nAssemblyId) + ".";
        
        Collection<BasicDBList> variantDataQueries = varQueryWrapper.getVariantDataQueries();
 
        CustomVCFWriter writer = null;		
		try {
			List<String> distinctSequenceNames = new ArrayList<String>();

//				String sequenceSeqCollName = MongoTemplateManager.getMongoCollectionName(Sequence.class);
//				if (mongoTemplate.collectionExists(sequenceSeqCollName))
//					try (MongoCursor<Document> markerCursor = IExportHandler.getSortedExportCursor(mongoTemplate, collWithPojoCodec, Document.class, varQuery, null, true, nQueryChunkSize)) {
//						while (markerCursor.hasNext())
//						{
//							int nLoadedMarkerCountInLoop = 0;
//							boolean fStartingNewChunk = true;
//							while (markerCursor.hasNext() && (fStartingNewChunk || nLoadedMarkerCountInLoop%nQueryChunkSize != 0)) {
//								Document exportVariant = markerCursor.next();
//								String chr = (String) ((Document) exportVariant.get(VariantData.FIELDNAME_REFERENCE_POSITION)).get(ReferencePosition.FIELDNAME_SEQUENCE);
//								if (chr != null && !distinctSequenceNames.contains(chr))
//									distinctSequenceNames.add(chr);
//							}
//						}
//					}
//					else {
					for (Object chr : mongoTemplate.getCollection(mongoTemplate.getCollectionName(VariantData.class)).distinct(refPosPathWithTrailingDot + ReferencePosition.FIELDNAME_SEQUENCE, !variantDataQueries.isEmpty() ? new BasicDBObject("$and", variantDataQueries.iterator().next()) : new BasicDBObject(), String.class))	// find out distinctSequenceNames by looking at exported variant list
						if (chr != null)
							distinctSequenceNames.add(chr.toString());
//					}

			Collections.sort(distinctSequenceNames, new AlphaNumericComparator());
			SAMSequenceDictionary dict = createSAMSequenceDictionary(sModule, distinctSequenceNames);
			writer = buildVariantContextWriter(zos, dict);
			zos.putNextEntry(new ZipEntry(exportName + "." + getExportDataFileExtensions()[0]));
			Document variantQueryForTargetCollection = variantDataQueries.isEmpty() ? new Document() : (tmpVarCollName == null ? new Document("$and", variantDataQueries.iterator().next()) : (varQueryWrapper.getBareQueries().iterator().hasNext() ? new Document("$and", varQueryWrapper.getBareQueries().iterator().next()) : new Document()));

			File[] warningFiles = writeGenotypeFile(sModule, nAssemblyId, individualsByPop, workWithSamples, annotationFieldThresholds, progress, tmpVarCollName, variantQueryForTargetCollection, markerCount, markerSynonyms, callSetsToExport, sortedIndividuals, distinctSequenceNames, dict, writer);
			zos.closeEntry();

			IExportHandler.writeZipEntryFromChunkFiles(zos, warningFiles, exportName + "-REMARKS.txt");
		}
		catch (Exception e) {
			LOG.error("Error exporting", e);
			progress.setError(e.getMessage());
			return;
		}
		finally {
			if (writer != null)
				try {
					zos.finish();
					writer.close();
				}
				catch (Throwable ignored)
				{}
		}
    }
		
    public File[] writeGenotypeFile(String sModule, Integer nAssemblyId, Map<String, Collection<String>> individualsByPop, boolean workWithSamples, Map<String, HashMap<String, Float>> annotationFieldThresholds, ProgressIndicator progress, String tmpVarCollName, Document variantQuery, long markerCount, Map<String, String> markerSynonyms, Collection<Callset> callSetsToExport, List<String> sortedIndividuals, List<String> distinctSequenceNames, SAMSequenceDictionary dict, CustomVCFWriter writer) throws Exception {
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

//			VariantContextWriterBuilder vcwb = new VariantContextWriterBuilder();
//			vcwb.unsetOption(Options.INDEX_ON_THE_FLY);
//			vcwb.unsetOption(Options.DO_NOT_WRITE_GENOTYPES);
//			vcwb.setOption(Options.USE_ASYNC_IOINDEX_ON_THE_FLY);
//			vcwb.setOption(Options.ALLOW_MISSING_FIELDS_IN_HEADER);
//			vcwb.setReferenceDictionary(dict);
//			writer = vcwb.build();
//			writer = new AsyncVariantContextWriter(writer, 3000);

		progress.moveToNextStep();	// done with dictionary
		MongoCollection<org.bson.Document> vcfHeaderColl = mongoTemplate.getCollection(MongoTemplateManager.getMongoCollectionName(DBVCFHeader.class));
		Document vcfHeaderQuery = new Document("_id." + VcfHeaderId.FIELDNAME_PROJECT, projectId);
		long nHeaderCount = vcfHeaderColl.countDocuments(vcfHeaderQuery);
		MongoCursor<Document> headerCursor = vcfHeaderColl.find(vcfHeaderQuery).iterator();
		Set<VCFHeaderLine> headerLines = new HashSet<VCFHeaderLine>();
		boolean fWriteCommandLine = true, fWriteEngineHeaders = true;	// default values

		while (headerCursor.hasNext()) {
			DBVCFHeader dbVcfHeader = DBVCFHeader.fromDocument(headerCursor.next());
			headerLines.addAll(dbVcfHeader.getHeaderLines());

			fWriteCommandLine = nHeaderCount == 1 && dbVcfHeader.getWriteCommandLine();	// wouldn't make sense to include command lines for several runs
			if (!dbVcfHeader.getWriteEngineHeaders())
				fWriteEngineHeaders = false;
		}
		headerCursor.close();

		if (headerLines.size() == 0)
			headerLines.add(new VCFFormatHeaderLine("GT", 1, VCFHeaderLineType.String, "Genotype"));	// minimum required
		
		// Add sequence header lines (not stored in our vcf header collection)
		int nSequenceIndex = 0;
		String sequenceInfoCollName = MongoTemplateManager.getMongoCollectionName(SequenceStats.class);
		boolean fCollectionExists = mongoTemplate.collectionExists(sequenceInfoCollName);
		for (String sequenceName : distinctSequenceNames) {
			Map<String, String> sequenceLineData = new LinkedHashMap<String, String>();
			if (fCollectionExists) { /* this should not happen anymore (relates to Gigwa V1's contig management section) */
				Document record = mongoTemplate.getCollection(sequenceInfoCollName).find(new Query(Criteria.where("_id").is(sequenceName)).getQueryObject()).projection(new Document(SequenceStats.FIELDNAME_SEQUENCE_LENGTH, true)).first();
				if (record == null) {
					LOG.warn("Sequence '" + sequenceName + "' not found in collection " + sequenceInfoCollName);
					continue;
				}
				sequenceLineData.put("ID", (String) record.get("_id"));
				sequenceLineData.put("length", ((Number) record.get(SequenceStats.FIELDNAME_SEQUENCE_LENGTH)).toString());
			}
			else {
				sequenceLineData.put("ID", sequenceName);
				SAMSequenceRecord dictSeq = dict.getSequence(sequenceName); 
				if (dictSeq != null) {
				    if (dictSeq.getSequenceLength() > 0)
	                    sequenceLineData.put("length", "" + dict.getSequence(sequenceName).getSequenceLength());
				    String sAssembly = dictSeq.getAssembly();
                    if (sAssembly != null)
                        sequenceLineData.put("assembly", sAssembly);
				}
			}
			headerLines.add(new VCFContigHeaderLine(sequenceLineData, nSequenceIndex++));
		}

		VCFHeader header = new VCFHeader(headerLines, sortedIndividuals);
		header.setWriteCommandLine(fWriteCommandLine);
		header.setWriteEngineHeaders(fWriteEngineHeaders);
		writer.writeHeader(header);

		HashMap<Integer, Object /*phID*/> phasingIDsByCallSet = new HashMap<>();
		ExportManager.AbstractExportWriter writingThread = new ExportManager.AbstractExportWriter() {
			public void writeChunkRuns(Collection<Collection<VariantRunData>> markerRunsToWrite, List<String> orderedMarkerIDs, OutputStream genotypeOS, OutputStream variantOS, OutputStream warningOS) throws IOException {
				if (markerRunsToWrite.isEmpty())
					return;
				
				CustomVCFWriter chunkWriter = new CustomVCFWriter(null, genotypeOS, null, false, false, true);
				chunkWriter.accountForHeaderWithoutWritingIt(header);

			    ArrayList<ArrayList<Collection<VariantRunData>>> splitVrdColls = Helper.evenlySplitCollection(markerRunsToWrite, Runtime.getRuntime().availableProcessors() - 1);
			    ArrayList<ArrayList<String>> splitVariantIdColls = Helper.evenlySplitCollection(orderedMarkerIDs, splitVrdColls.size());
			    
                VariantContext[][] vcChunks = new VariantContext[splitVrdColls.size()][];
		        ExecutorService executor = Executors.newFixedThreadPool(vcChunks.length);
		        for (int i=0; i<vcChunks.length; i++) {
		            final ArrayList<Collection<VariantRunData>> vrdCollChunk = splitVrdColls.get(i);
		            final ArrayList<String> variantIdCollChunk = splitVariantIdColls.get(i);
		            final int nChunkIndex = i;
		            Thread t = new Thread() {
		                public void run() {
                            vcChunks[nChunkIndex] = new VariantContext[vrdCollChunk.size()];		                    
                            int nVariantIndex = 0;
		                    for (Collection<VariantRunData> runsToWrite : vrdCollChunk) {
		                    	String idOfVarToWrite = variantIdCollChunk.get(nVariantIndex);
		    					if (progress.isAborted() || progress.getError() != null)
		    						return;
		    					
		    					AbstractVariantData variant = runsToWrite == null || runsToWrite.isEmpty() ? mongoTemplate.findById(idOfVarToWrite, VariantData.class) : runsToWrite.iterator().next();
            					try
            					{
            		                vcChunks[nChunkIndex][nVariantIndex++] = variant.toVariantContext(mongoTemplate, runsToWrite, nAssemblyId, !MgdbDao.idLooksGenerated(idOfVarToWrite), callSetsToExport, individualPositions, individualsByPop, workWithSamples, annotationFieldThresholds, phasingIDsByCallSet, warningOS, markerSynonyms == null ? idOfVarToWrite : markerSynonyms.get(idOfVarToWrite));
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
		            executor.execute(t);
		        }
		        executor.shutdown();
		        try {
                    executor.awaitTermination(Integer.MAX_VALUE, TimeUnit.DAYS);
                } catch (InterruptedException e) {
                    progress.setError(e.getMessage());
                    LOG.error("Error exporting VCF", e);
                }
		        
		        for (VariantContext[] vcChunk : vcChunks)
		            for (VariantContext vc : vcChunk)
		            	chunkWriter.add(vc);
		        chunkWriter.close();
			}
		};

		int nQueryChunkSize = IExportHandler.computeQueryChunkSize(mongoTemplate, markerCount);
        String usedCollName = tmpVarCollName != null ? tmpVarCollName : mongoTemplate.getCollectionName(VariantRunData.class);
		MongoCollection collWithPojoCodec = mongoTemplate.getDb().withCodecRegistry(ExportManager.pojoCodecRegistry).getCollection(usedCollName);
		ExportManager exportManager = new ExportManager(sModule, nAssemblyId, collWithPojoCodec, VariantRunData.class, variantQuery, callSetsToExport, true, nQueryChunkSize, writingThread, markerCount, progress);
		if (tmpFolderPath != null)
			exportManager.setTmpExtractionFolder(tmpFolderPath + File.separator + Helper.convertToMD5(progress.getProcessId()));
		exportManager.readAndWrite(writer.getOutputStream());
		return exportManager.getOutputs().getWarningFiles();
	}

    @Override
	public String[] getExportDataFileExtensions() {
		return new String[] {"vcf"};
	}

    @Override
    public int[] getSupportedPloidyLevels() {
        return null;
    }
}