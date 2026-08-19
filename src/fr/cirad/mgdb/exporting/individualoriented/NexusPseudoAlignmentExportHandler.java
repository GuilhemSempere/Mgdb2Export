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

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import htsjdk.variant.variantcontext.VariantContext.Type;

/**
 * The Class NexusPseudoAlignmentExportHandler.
 */
public class NexusPseudoAlignmentExportHandler extends PhylipPseudoAlignmentExportHandler {

    /**
     * The Constant LOG.
     */
    private static final Logger LOG = LoggerFactory.getLogger(NexusPseudoAlignmentExportHandler.class);

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
        return "NEXUS";
    }

    /* (non-Javadoc)
	 * @see fr.cirad.mgdb.exporting.IExportHandler#getExportFormatDescription()
     */
    @Override
    public String getExportFormatDescription() {
    	return "Exports a zipped NEXUS file containing a pseudo-alignment consisting in the concatenation of SNP alleles, compatible with tree construction tools like IQ-TREE. An additional PLINK-style map file is added for reference.";
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
    protected String getMissingAlleleString() {
		return "?";
	}
    
    @Override
    protected String getLinePrefix() {
		return "";
	}
    
    @Override
    protected String getIndividualToSequenceSeparator() {
		return "\t";
	}
    
    @Override
    protected String getHeaderlines(int nIndividualCount, int nMarkerCount) {
    	return "#NEXUS\nbegin data;\ndimensions ntax=" + nIndividualCount + " nchar=" + 2*nMarkerCount + ";\nformat datatype=dna missing=? gap=-;\nmatrix\n";
	}

    @Override
    protected String getFooterlines() {
    	return ";\nend;";
	}

    /* (non-Javadoc)
	 * @see fr.cirad.mgdb.exporting.IExportHandler#getStepList()
     */
    @Override
    public List<String> getStepList() {
        return Arrays.asList(new String[]{"Exporting data to NEXUS format"});
    }

	@Override
	public String[] getExportDataFileExtensions() {
		return new String[] {"nxs"};
	}

    @Override
    public int[] getSupportedPloidyLevels() {
        return new int[] {2};
    }
}