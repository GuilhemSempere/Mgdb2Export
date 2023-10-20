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
import java.io.OutputStream;
import java.util.Arrays;
import java.util.List;
import java.util.zip.Deflater;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.samtools.util.zip.DeflaterFactory;
import htsjdk.variant.variantcontext.writer.CustomVCFWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;

/**
 * The Class VcfExportHandler.
 */
public class VcfGzExportHandler extends VcfExportHandler {

	/* (non-Javadoc)
	 * @see fr.cirad.mgdb.exporting.IExportHandler#getExportFormatName()
	 */
	@Override
	public String getExportFormatName()
	{
		return "VCF.gz";
	}

	/* (non-Javadoc)
	 * @see fr.cirad.mgdb.exporting.IExportHandler#getExportFormatDescription()
	 */
	@Override
	public String getExportFormatDescription()
	{
		return "Exports data in bgzipped Variant Call Format. See <a target='_blank' href='http://samtools.github.io/hts-specs/VCFv4.1.pdf'>http://samtools.github.io/hts-specs/VCFv4.1.pdf</a> and <a target='_blank' href='http://www.htslib.org/doc/bgzip.html'>http://www.htslib.org/doc/bgzip.html</a> for more details.";
	}

	/* (non-Javadoc)
	 * @see fr.cirad.mgdb.exporting.IExportHandler#getStepList()
	 */
	@Override
	public List<String> getStepList()
	{
		return Arrays.asList(new String[] {"Creating sequence list", "Exporting data to VCF.gz format"});
	}
	
	@Override
	protected VariantContextWriter buildVariantContextWriter(OutputStream os, SAMSequenceDictionary dict) {
		return new CustomVCFWriter(null, new BlockCompressedOutputStream(os, (File) null, Deflater.BEST_COMPRESSION, new DeflaterFactory()), dict, false, false, true);
	}

    @Override
	public String[] getExportDataFileExtensions() {
		return new String[] {"vcf.gz"};
	}
}