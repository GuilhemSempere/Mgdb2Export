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
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;

import org.apache.log4j.Logger;
import org.bson.Document;

import com.mongodb.client.MongoCollection;
import com.mongodb.client.MongoCursor;

import fr.cirad.mgdb.exporting.IExportHandler;
import fr.cirad.mgdb.model.mongo.maintypes.GenotypingSample;
import fr.cirad.mgdb.model.mongo.maintypes.VariantData;
import fr.cirad.mgdb.model.mongo.maintypes.VariantDataV2;
import fr.cirad.mgdb.model.mongo.subtypes.ReferencePosition;
import fr.cirad.tools.Helper;
import fr.cirad.tools.ProgressIndicator;

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

	/* (non-Javadoc)
	 * @see fr.cirad.mgdb.exporting.markeroriented.AbstractMarkerOrientedExportHandler#exportData(java.io.OutputStream, java.lang.String, java.util.List, fr.cirad.tools.ProgressIndicator, com.mongodb.DBCursor, java.util.Map, int, int, java.util.Map)
	 */
	@Override
	public void exportData(OutputStream outputStream, String sModule, Integer nAssemblyId, Collection<GenotypingSample> samples1, Collection<GenotypingSample> samples2, ProgressIndicator progress, MongoCollection<Document> varColl, Document varQuery, Map<String, String> markerSynonyms, HashMap<String, Float> annotationFieldThresholds, HashMap<String, Float> annotationFieldThresholds2, List<GenotypingSample> samplesToExport, Map<String, InputStream> readyToExportFiles) throws Exception {
		ZipOutputStream zos = new ZipOutputStream(outputStream);
		
		if (readyToExportFiles != null)
			for (String readyToExportFile : readyToExportFiles.keySet())
			{
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

		long markerCount = varColl.countDocuments(varQuery);
		String exportName = sModule + "__" + markerCount + "variants";
		zos.putNextEntry(new ZipEntry(exportName + ".bed"));
		
		short nProgress = 0, nPreviousProgress = 0;	
		int nQueryChunkSize = (int) Math.min(2000, markerCount), nLoadedMarkerCount = 0;
		try (MongoCursor<Document> markerCursor = IExportHandler.getMarkerCursorWithCorrectCollation(varColl, varQuery, nAssemblyId, nQueryChunkSize)) {
			while (markerCursor.hasNext())
			{
				int nLoadedMarkerCountInLoop = 0;
				Map<String, String> markerChromosomalPositions = new LinkedHashMap<>();
				boolean fStartingNewChunk = true;
				while (markerCursor.hasNext() && (fStartingNewChunk || nLoadedMarkerCountInLoop%nQueryChunkSize != 0)) {
					Document exportVariant = markerCursor.next();
					Document refPos = (Document) Helper.readPossiblyNestedField(exportVariant, (nAssemblyId != null ? VariantData.FIELDNAME_REFERENCE_POSITION + "." + nAssemblyId : VariantDataV2.FIELDNAME_REFERENCE_POSITION), "; ", null);
					markerChromosomalPositions.put((String) exportVariant.get("_id"), refPos == null ? null : (refPos.get(ReferencePosition.FIELDNAME_SEQUENCE) + ":" + refPos.get(ReferencePosition.FIELDNAME_START_SITE)));
					nLoadedMarkerCountInLoop++;
					fStartingNewChunk = false;
				}
	
				for (String variantId : markerChromosomalPositions.keySet())
				{
					String refPos = markerChromosomalPositions.get(variantId);
					if (refPos != null)
					{
						String[] chromAndPos = refPos.split(":");
						zos.write((chromAndPos[0] + "\t" + (Long.parseLong(chromAndPos[1])-1)  + "\t" + (Long.parseLong(chromAndPos[1])-1) + "\t" + variantId + "\t" + "0" + "\t" + "+").getBytes());
					}
					else
						zos.write(("0\t0\t0\t" + variantId + "\t" + "0" + "\t" + "+").getBytes());
					zos.write((LINE_SEPARATOR).getBytes());
				}
				
	            if (progress.isAborted())
	            	return;
	
	            nLoadedMarkerCount += nLoadedMarkerCountInLoop;
				nProgress = (short) (nLoadedMarkerCount * 100 / markerCount);
				if (nProgress > nPreviousProgress)
				{
					progress.setCurrentStepProgress(nProgress);
					nPreviousProgress = nProgress;
				}	
	        }
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
}