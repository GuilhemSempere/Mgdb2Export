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

import java.io.InputStream;
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
import fr.cirad.mgdb.model.mongo.maintypes.Assembly;
import fr.cirad.mgdb.model.mongo.maintypes.GenotypingSample;
import fr.cirad.mgdb.model.mongo.maintypes.VariantData;
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
	public void exportData(OutputStream outputStream, String sModule, Integer nAssemblyId, String sExportingUser, ProgressIndicator progress, String tmpVarCollName, VariantQueryWrapper varQueryWrapper, long markerCount, Map<String, String> markerSynonyms, Map<String, Collection<String>> individuals, Map<String, HashMap<String, Float>> annotationFieldThresholds, Collection<GenotypingSample> samplesToExport, Collection<String> individualMetadataFieldsToExport, Map<String, InputStream> readyToExportFiles) throws Exception {
        ZipOutputStream zos = IExportHandler.createArchiveOutputStream(outputStream, readyToExportFiles);
        MongoTemplate mongoTemplate = MongoTemplateManager.get(sModule);

        Assembly assembly = mongoTemplate.findOne(new Query(Criteria.where("_id").is(nAssemblyId)), Assembly.class);
		String exportName = sModule + (assembly != null && assembly.getName() != null ? "__" + assembly.getName() : "") + "__" + markerCount + "variants";
		zos.putNextEntry(new ZipEntry(exportName + ".bed"));
		
        Collection<BasicDBList> variantDataQueries = varQueryWrapper.getVariantDataQueries();
        BasicDBObject varQuery = !variantDataQueries.isEmpty() ? new BasicDBObject("$and", variantDataQueries.iterator().next()) : new BasicDBObject();
		
		short nProgress = 0, nPreviousProgress = 0;
		int nQueryChunkSize = (int) Math.min(2000, markerCount), nLoadedMarkerCount = 0;
		String refPosPath = Assembly.getVariantRefPosPath(nAssemblyId);
    	String refPosPathWithTrailingDot = Assembly.getThreadBoundVariantRefPosPath() + ".";
    	Document projectionAndSortDoc = new Document(refPosPathWithTrailingDot + ReferencePosition.FIELDNAME_SEQUENCE, 1).append(refPosPathWithTrailingDot + ReferencePosition.FIELDNAME_START_SITE, 1).append(refPosPathWithTrailingDot + ReferencePosition.FIELDNAME_END_SITE, 1);
    	try (MongoCursor<Document> markerCursor = IExportHandler.getMarkerCursorWithCorrectCollation(mongoTemplate.getDb().withCodecRegistry(ExportManager.pojoCodecRegistry).getCollection(tmpVarCollName != null ? tmpVarCollName : mongoTemplate.getCollectionName(VariantData.class)), tmpVarCollName != null ? new Document() : new Document(varQuery), projectionAndSortDoc, nQueryChunkSize)) {
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
	
	@Override
	public String[] getExportDataFileExtensions() {
		return new String[] {"bed"};
	}

    @Override
    public int[] getSupportedPloidyLevels() {
        return null;
    }
}