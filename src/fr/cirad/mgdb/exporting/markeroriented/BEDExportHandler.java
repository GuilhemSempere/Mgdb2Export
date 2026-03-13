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
import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;

import org.apache.log4j.Logger;
import org.bson.Document;
import org.springframework.data.mongodb.core.MongoTemplate;
import org.springframework.data.mongodb.core.query.Criteria;
import org.springframework.data.mongodb.core.query.Query;

import com.mongodb.BasicDBList;
import com.mongodb.BasicDBObject;
import com.mongodb.client.MongoCursor;

import fr.cirad.mgdb.exporting.IExportHandler;
import fr.cirad.mgdb.exporting.tools.ExportManager;
import fr.cirad.mgdb.exporting.tools.ExportManager.ExportOutputs;
import fr.cirad.mgdb.model.mongo.maintypes.Assembly;
import fr.cirad.mgdb.model.mongo.maintypes.DBVCFHeader;
import fr.cirad.mgdb.model.mongo.maintypes.VariantData;
import fr.cirad.mgdb.model.mongo.subtypes.AbstractVariantData;
import fr.cirad.mgdb.model.mongo.subtypes.Callset;
import fr.cirad.mgdb.model.mongo.subtypes.ReferencePosition;
import fr.cirad.tools.Helper;
import fr.cirad.tools.ProgressIndicator;
import fr.cirad.tools.mgdb.VariantQueryWrapper;
import fr.cirad.tools.mongo.MongoTemplateManager;

/**
 * The Class BEDExportHandler.
 */
public class BEDExportHandler extends AbstractMarkerOrientedExportHandler
{
	
	/** The Constant LOG. */
	static final Logger LOG = Logger.getLogger(BEDExportHandler.class);
	
	/* (non-Javadoc)
	 * @see fr.cirad.mgdb.exporting.IExportHandler#getExportFormatName()
	 */
	@Override
	public String getExportFormatName() {
		return "BED";
	}

	/* (non-Javadoc)
	 * @see fr.cirad.mgdb.exporting.IExportHandler#getExportFormatDescription()
	 */
	@Override
	public String getExportFormatDescription() {
		return "Exports data in BED Format. See <a target='_blank' href='http://genome.ucsc.edu/FAQ/FAQformat.html#format1'>http://genome.ucsc.edu/FAQ/FAQformat.html#format1</a> for more details";
	}
	
	@Override
	public String getExportArchiveExtension() {
		return "zip";
	}
	
	/* (non-Javadoc)
	 * @see fr.cirad.mgdb.exporting.IExportHandler#getStepList()
	 */
	@Override
	public List<String> getStepList() {
		return Arrays.asList(new String[] {"Exporting data to BED format"});
	}

	@Override
	public void exportData(OutputStream outputStream, String sModule, Integer nAssemblyId, String sExportingUser, ProgressIndicator progress, String tmpVarCollName, VariantQueryWrapper varQueryWrapper, long markerCount, Map<String, String> markerSynonyms, Map<String, Collection<String>> individualsByPop, boolean workWithSamples, Map<String, HashMap<String, Float>> annotationFieldThresholds, Collection<Callset> callSetsToExport, Collection<String> individualMetadataFieldsToExport, Map<String, InputStream> readyToExportFiles) throws Exception {
        ZipOutputStream zos = IExportHandler.createArchiveOutputStream(outputStream, readyToExportFiles, null);
        MongoTemplate mongoTemplate = MongoTemplateManager.get(sModule);

        Assembly assembly = mongoTemplate.findOne(new Query(Criteria.where("_id").is(nAssemblyId)), Assembly.class);
		String exportName = sModule + (assembly != null && assembly.getName() != null ? "__" + assembly.getName() : "") + "__" + markerCount + "variants";
		
        Collection<BasicDBList> variantDataQueries = varQueryWrapper.getVariantDataQueries();
        Document variantQueryForTargetCollection = variantDataQueries.isEmpty() ? new Document() : (tmpVarCollName == null ? new Document("$and", variantDataQueries.iterator().next()) : (varQueryWrapper.getBareQueries().iterator().hasNext() ? new Document("$and", varQueryWrapper.getBareQueries().iterator().next()) : new Document()));
		
        BasicDBObject varAnnField = new BasicDBObject();
        varAnnField.put(AbstractVariantData.VCF_CONSTANT_INFO_META_DATA, 1);
        varAnnField.put(AbstractVariantData.VCF_CONSTANT_INFO_FORMAT_META_DATA, 1);

        if (null != MongoTemplateManager.get(sModule).getCollection(MongoTemplateManager.getMongoCollectionName(DBVCFHeader.class)).find().projection(varAnnField).first()) {
        	// we have annotations to export so let's fake a hapmap export and build BED format from that 
            ByteArrayOutputStream baos = new ByteArrayOutputStream();
            HapMapExportHandler heh = (HapMapExportHandler) AbstractMarkerOrientedExportHandler.getMarkerOrientedExportHandlers().get("HAPMAP");
            ExportOutputs outputs = heh.writeGenotypeFile(true, true, true, true, baos, sModule, mongoTemplate.findOne(new Query(Criteria.where("_id").is(Assembly.getThreadBoundAssembly())), Assembly.class), new HashMap<>(), workWithSamples, new HashMap<>(), new HashMap<>(), progress, tmpVarCollName != null ? tmpVarCollName : null, variantQueryForTargetCollection, markerCount, null, callSetsToExport);
            zos.putNextEntry(new ZipEntry(exportName + ".bed"));
            try (BufferedReader reader = new BufferedReader(new InputStreamReader(new ByteArrayInputStream(baos.toByteArray())))) {
            	int nLinesRead = 0;
                String line;
                while ((line = reader.readLine()) != null) {
                	if (nLinesRead > 0)	{ // skip header
                		String[] splitLine = line.split("\t");
                		int start = Integer.parseInt(splitLine[3]) - 1;
                		int end = start == -1 ? 0 : (start + splitLine[1].split("/")[0].length() - 1);
                		zos.write((splitLine[2] + "\t" + Math.max(0, start) + "\t" + end + "\t" + splitLine[0] + "\t0\t+\n").getBytes());
                	}
                	nLinesRead++;
                }
            }
			zos.closeEntry();

            IExportHandler.writeZipEntryFromChunkFiles(zos, outputs.getAnnotationFiles(), exportName + ".ann", IExportHandler.VEP_LIKE_HEADER_LINE);
			zos.finish();
	        zos.close();
			progress.setCurrentStepProgress((short) 100);
        }
        else {	// quick mode: don't use VariantRunData at all
			short nProgress = 0, nPreviousProgress = 0;
			int nQueryChunkSize = (int) Math.min(2000, markerCount), nLoadedMarkerCount = 0;
			String refPosPath = Assembly.getVariantRefPosPath(nAssemblyId);
	    	String refPosPathWithTrailingDot = Assembly.getThreadBoundVariantRefPosPath() + ".";
	    	Document projectionAndSortDoc = new Document(refPosPathWithTrailingDot + ReferencePosition.FIELDNAME_SEQUENCE, 1).append(refPosPathWithTrailingDot + ReferencePosition.FIELDNAME_START_SITE, 1);
	    	
			zos.putNextEntry(new ZipEntry(exportName + ".bed"));
	    	try (MongoCursor<Document> markerCursor = IExportHandler.getMarkerCursorWithCorrectCollation(mongoTemplate.getDb().withCodecRegistry(ExportManager.pojoCodecRegistry).getCollection(tmpVarCollName != null ? tmpVarCollName : mongoTemplate.getCollectionName(VariantData.class)), variantQueryForTargetCollection, projectionAndSortDoc, nQueryChunkSize)) {
				while (markerCursor.hasNext())
				{
					Document exportVariant = markerCursor.next();
					Document refPos = (Document) Helper.readPossiblyNestedField(exportVariant, refPosPath, ";", null);
					String chrom = refPos == null ? null : (String) refPos.get(ReferencePosition.FIELDNAME_SEQUENCE);
					if (chrom != null) {
						Long pos = chrom == null ? null : ((Number) refPos.get(ReferencePosition.FIELDNAME_START_SITE)).longValue();
						Number end = chrom == null ? null : (Number) refPos.get(ReferencePosition.FIELDNAME_END_SITE);
						zos.write((chrom + "\t" + (pos - 1) + "\t" + ((end == null ? pos : end.longValue()) - 1) + "\t" + exportVariant.get("_id") + "\t" + "0" + "\t" + "+").getBytes());
					}
					else
						zos.write(("0\t0\t0\t" + exportVariant.get("_id") + "\t" + "0" + "\t" + "+").getBytes());
					zos.write((LINE_SEPARATOR).getBytes());
					
					nLoadedMarkerCount++;
					nProgress = (short) (nLoadedMarkerCount * 100 / markerCount);
					if (nProgress > nPreviousProgress) {
						progress.setCurrentStepProgress(nProgress);
						nPreviousProgress = nProgress;
					}
				}
	
	            if (progress.isAborted())
	            	return;
			}
			zos.closeEntry();
			zos.finish();
	        zos.close();
			progress.setCurrentStepProgress((short) 100);
        }
	}
	
	@Override
	public String[] getExportDataFileExtensions() {
		return new String[] {"bed"};
	}

    @Override
    public int[] getSupportedPloidyLevels() {
        return null;
    }
}