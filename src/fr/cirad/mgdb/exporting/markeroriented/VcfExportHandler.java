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
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.atomic.AtomicLong;
import java.util.stream.Collectors;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;

import org.apache.log4j.Logger;
import org.bson.Document;
import org.springframework.data.mongodb.core.MongoTemplate;
import org.springframework.data.mongodb.core.mapreduce.MapReduceResults;
import org.springframework.data.mongodb.core.query.Criteria;
import org.springframework.data.mongodb.core.query.Query;

import com.mongodb.client.MongoCollection;
import com.mongodb.client.MongoCursor;

import fr.cirad.mgdb.exporting.AbstractExportWritingThread;
import fr.cirad.mgdb.exporting.IExportHandler;
import fr.cirad.mgdb.exporting.tools.ExportManager;
import fr.cirad.mgdb.model.mongo.maintypes.DBVCFHeader;
import fr.cirad.mgdb.model.mongo.maintypes.DBVCFHeader.VcfHeaderId;
import fr.cirad.mgdb.model.mongo.maintypes.GenotypingSample;
import fr.cirad.mgdb.model.mongo.maintypes.Sequence;
import fr.cirad.mgdb.model.mongo.maintypes.SequenceStats;
import fr.cirad.mgdb.model.mongo.maintypes.VariantData;
import fr.cirad.mgdb.model.mongo.maintypes.VariantRunData;
import fr.cirad.mgdb.model.mongo.subtypes.ReferencePosition;
import fr.cirad.mgdb.model.mongodao.MgdbDao;
import fr.cirad.tools.AlphaNumericComparator;
import fr.cirad.tools.ProgressIndicator;
import fr.cirad.tools.mongo.MongoTemplateManager;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.CustomVCFWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
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
	 * @param sequenceIDs the sequence i ds
	 * @return the SAM sequence dictionary
	 * @throws Exception the exception
	 */
	public SAMSequenceDictionary createSAMSequenceDictionary(String sModule, Collection<String> sequenceIDs) throws Exception
	{
		SAMSequenceDictionary dict1 = new SAMSequenceDictionary();
		MongoTemplate mongoTemplate = MongoTemplateManager.get(sModule);
		String sequenceSeqCollName = MongoTemplateManager.getMongoCollectionName(Sequence.class);
		if (mongoTemplate.collectionExists(sequenceSeqCollName) && sequenceIDs.size() > 1)
		{
			long before = System.currentTimeMillis();
			Query query = new Query(Criteria.where("_id").in(sequenceIDs));
			String mapFunction = "function() { emit(this._id, this." + Sequence.FIELDNAME_SEQUENCE + ".length); }";
			String reduceFunction = "function() { }";
			MapReduceResults<Map> rs = MongoTemplateManager.get(sModule).mapReduce(query, sequenceSeqCollName, mapFunction, reduceFunction, Map.class);
			Iterator<Map> it = rs.iterator();
			while (it.hasNext())
			{
				Map map = it.next();
				dict1.addSequence(new SAMSequenceRecord((String) map.get("_id"), ((Double) map.get("value")).intValue()));
			}
	    	LOG.info("createSAMSequenceDictionary took " + (System.currentTimeMillis() - before)/1000d + "s to write " + sequenceIDs.size() + " sequences");
		}
		else
			LOG.info("No sequence data was found. No SAMSequenceDictionary could be created.");
	    return dict1;
	}

    @Override
    public void exportData(OutputStream outputStream, String sModule, Collection<String> individuals1, Collection<String> individuals2, ProgressIndicator progress, String tmpVarCollName, Document varQuery, long markerCount, Map<String, String> markerSynonyms, HashMap<String, Float> annotationFieldThresholds, HashMap<String, Float> annotationFieldThresholds2, List<GenotypingSample> samplesToExport, Map<String, InputStream> readyToExportFiles) throws Exception {
		List<String> sortedIndividuals = samplesToExport.stream().map(gs -> gs.getIndividual()).distinct().sorted(new AlphaNumericComparator<String>()).collect(Collectors.toList());

		MongoTemplate mongoTemplate = MongoTemplateManager.get(sModule);
        File warningFile = File.createTempFile("export_warnings_", "");
        FileWriter warningFileWriter = new FileWriter(warningFile);
        ZipOutputStream zos = IExportHandler.createArchiveOutputStream(outputStream, readyToExportFiles);
		String exportName = sModule + "__" + markerCount + "variants__" + sortedIndividuals.size() + "individuals";

    	VariantContextWriter writer = null;		
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
					for (Object chr : mongoTemplate.getCollection(mongoTemplate.getCollectionName(VariantData.class)).distinct(VariantData.FIELDNAME_REFERENCE_POSITION + "." + ReferencePosition.FIELDNAME_SEQUENCE, varQuery, String.class))	// find out distinctSequenceNames by looking at exported variant list
						if (chr != null)
							distinctSequenceNames.add(chr.toString());
//					}

			Collections.sort(distinctSequenceNames, new AlphaNumericComparator());
			SAMSequenceDictionary dict = createSAMSequenceDictionary(sModule, distinctSequenceNames);
			writer = new CustomVCFWriter(null, zos, dict, false, false, true);
			zos.putNextEntry(new ZipEntry(exportName + ".vcf"));
			writeGenotypeFile(sModule, individuals1, individuals2, progress, tmpVarCollName, varQuery, markerCount, markerSynonyms, annotationFieldThresholds, annotationFieldThresholds2, samplesToExport, sortedIndividuals, distinctSequenceNames, warningFileWriter, writer);
		}
		catch (Exception e)
		{
			LOG.error("Error exporting", e);
			progress.setError(e.getMessage());
			return;
		}
		finally
		{
			zos.closeEntry();

			warningFileWriter.close();
			if (warningFile.length() > 0) {
				zos.putNextEntry(new ZipEntry(exportName + "-REMARKS.txt"));
				int nWarningCount = 0;
				BufferedReader in = new BufferedReader(new FileReader(warningFile));
				String sLine;
				while ((sLine = in.readLine()) != null) {
					outputStream.write((sLine + "\n").getBytes());
					nWarningCount++;
				}
				LOG.info("Number of Warnings for export (" + exportName + "): " + nWarningCount);
				in.close();
				zos.closeEntry();
			}
			warningFile.delete();

			if (writer != null)
				try
				{
					zos.finish();
					writer.close();
				}
				catch (Throwable ignored)
				{}
		}
    }
		
    public void writeGenotypeFile(String sModule, Collection<String> individuals1, Collection<String> individuals2, ProgressIndicator progress, String tmpVarCollName, Document varQuery, long markerCount, Map<String, String> markerSynonyms, HashMap<String, Float> annotationFieldThresholds, HashMap<String, Float> annotationFieldThresholds2, List<GenotypingSample> samplesToExport, List<String> sortedIndividuals, List<String> distinctSequenceNames, FileWriter warningFileWriter, VariantContextWriter writer) throws Exception {
    	MongoTemplate mongoTemplate = MongoTemplateManager.get(sModule);
		Integer projectId = null;
		for (GenotypingSample sample : samplesToExport)
		{
			if (projectId == null)
				projectId = sample.getProjectId();
			else if (projectId != sample.getProjectId())
			{
				projectId = 0;
				break;	// more than one project are involved: no header will be written
			}
		}

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
			if (fCollectionExists) {
				Document record = mongoTemplate.getCollection(sequenceInfoCollName).find(new Query(Criteria.where("_id").is(sequenceName)).getQueryObject()).projection(new Document(SequenceStats.FIELDNAME_SEQUENCE_LENGTH, true)).first();
				if (record == null) {
					LOG.warn("Sequence '" + sequenceName + "' not found in collection " + sequenceInfoCollName);
					continue;
				}

				sequenceLineData.put("ID", (String) record.get("_id"));
				sequenceLineData.put("length", ((Number) record.get(SequenceStats.FIELDNAME_SEQUENCE_LENGTH)).toString());
			}
			else
				sequenceLineData.put("ID", sequenceName);
			headerLines.add(new VCFContigHeaderLine(sequenceLineData, nSequenceIndex++));
		}

		VCFHeader header = new VCFHeader(headerLines, sortedIndividuals);
		header.setWriteCommandLine(fWriteCommandLine);
		header.setWriteEngineHeaders(fWriteEngineHeaders);
		writer.writeHeader(header);

		HashMap<Integer, Object /*phID*/> phasingIDsBySample = new HashMap<>();
		
		final VariantContextWriter finalVariantContextWriter = writer;
//			AtomicLong timeConverting = new AtomicLong(0), timeWriting = new AtomicLong(0);
		AbstractExportWritingThread writingThread = new AbstractExportWritingThread() {
			public void run() {
				for (String idOfVarToWrite : markerRunsToWrite.keySet()) {
					if (progress.isAborted() || progress.getError() != null)
						return;

					List<VariantRunData> runsToWrite = markerRunsToWrite.get(idOfVarToWrite);
					if (runsToWrite.isEmpty())
						continue;
					
					String variantId = null;
					try
					{
//							long b4 = System.currentTimeMillis();
//			                if (markerSynonyms != null) {
//			                	String syn = markerSynonyms.get(idOfVarToWrite);
//			                    if (syn != null)
//			                    	idOfVarToWrite = syn;
//			                }
	
		                VariantRunData vrd = runsToWrite.get(0);
						variantId = vrd.getId().getVariantId();
						
		                if (markerSynonyms != null) {
		                	String syn = markerSynonyms.get(variantId);
		                    if (syn != null)
		                        variantId = syn;
		                }

						VariantContext vc = vrd.toVariantContext(markerRunsToWrite.get(idOfVarToWrite), !MgdbDao.idLooksGenerated(variantId.toString()), samplesToExport, individuals1, individuals2, phasingIDsBySample, annotationFieldThresholds, annotationFieldThresholds2, warningFileWriter, markerSynonyms == null ? variantId : markerSynonyms.get(variantId));
//							timeConverting.addAndGet(System.currentTimeMillis() - b4);
						
//							b4 = System.currentTimeMillis();
						finalVariantContextWriter.add(vc);
//							timeWriting.addAndGet(System.currentTimeMillis() - b4);
					}
					catch (Exception e)
					{
						if (progress.getError() == null)	// only log this once
							LOG.debug("Unable to export " + idOfVarToWrite, e);
						progress.setError("Unable to export " + idOfVarToWrite + ": " + e.getMessage());
					}
				}
				markerRunsToWrite.clear();
			}
		};
		
		int nQueryChunkSize = IExportHandler.computeQueryChunkSize(mongoTemplate, markerCount);
        String usedCollName = tmpVarCollName != null ? tmpVarCollName : mongoTemplate.getCollectionName(VariantRunData.class);
		MongoCollection collWithPojoCodec = mongoTemplate.getDb().withCodecRegistry(ExportManager.pojoCodecRegistry).getCollection(usedCollName);
		ExportManager exportManager = new ExportManager(mongoTemplate, collWithPojoCodec, VariantRunData.class, varQuery, samplesToExport, true, nQueryChunkSize, writingThread, markerCount, warningFileWriter, progress);
		exportManager.readAndWrite();

//			System.out.println("time spent converting: " + timeConverting.get());
//			System.out.println("time spent writing: " + timeWriting.get());
	}
    
	@Override
	public String[] getExportDataFileExtensions() {
		return new String[] {"vcf"};
	}
}