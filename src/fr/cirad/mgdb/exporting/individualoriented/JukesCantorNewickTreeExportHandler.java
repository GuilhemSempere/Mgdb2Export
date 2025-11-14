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
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;

import com.traviswheeler.ninja.TreeBuilderBinHeap;
import com.traviswheeler.ninja.TreeNode;

import fr.cirad.tools.ProgressIndicator;
import htsjdk.variant.variantcontext.VariantContext.Type;

/**
 * The Class JukesCantorNewickTreeExportHandler.
 */
public class JukesCantorNewickTreeExportHandler extends JukesCantorDistanceMatrixExportHandler {

    public JukesCantorNewickTreeExportHandler() {
    }
    
    public JukesCantorNewickTreeExportHandler(int nMaxMissingDataPercentageForIndividuals) {
    	super(nMaxMissingDataPercentageForIndividuals);
    }
    
    /* (non-Javadoc)
	 * @see fr.cirad.mgdb.exporting.IExportHandler#getExportFormatName()
     */
    @Override
    public String getExportFormatName() {
        return "NJ-NEWICK";
    }

    /* (non-Javadoc)
	 * @see fr.cirad.mgdb.exporting.IExportHandler#getExportFormatDescription()
     */
    @Override
    public String getExportFormatDescription() {
    	return "Exports a zipped <a target='_blank' href='https://en.wikipedia.org/wiki/Newick_format'>Newick</a> file containing a <a target='_blank' href='https://www.sglp.uzh.ch/apps/static/MLS/stemmatology/Jukes-Cantor-model_229150204.html'>Jukes-Cantor</a> neighbour-joining tree based on a pseudo-alignment consisting in the concatenation of SNP alleles. Individuals with more than 50% missing data are automatically excluded.";
    }

    @Override
    protected void finalizeExportUsingDistanceMatrix(List<String> sequenceNames, ArrayList<String> exportedIndividuals, String exportName, double[][] distanceMatrix, ZipOutputStream zos, ProgressIndicator progress) throws Exception {
    	progress.moveToNextStep();
        int[][] intDistanceMatrix = new int[distanceMatrix.length][];
        for (int i = 0; i < distanceMatrix.length; i++) {
        	intDistanceMatrix[i] = new int[distanceMatrix.length - 1 - i];
            for (int j = i + 1; j < distanceMatrix.length; j++)
            	intDistanceMatrix[i][j - 1 - i] = (int) (distanceMatrix[i][j] * 100000000);
        }

		StringBuffer sb = new StringBuffer(); 

		TreeBuilderBinHeap tb = new TreeBuilderBinHeap(sequenceNames.toArray(new String[sequenceNames.size()]), intDistanceMatrix);
		TreeNode[] nodes = tb.build();
//    		nodes[nodes.length-1].rootTreeAt("CR1062");
		nodes[nodes.length-1].buildTreeString(sb);

    	
        zos.putNextEntry(new ZipEntry(exportName + "." + getExportDataFileExtensions()[0]));
		String treeString = sb.toString() + ";";
    	zos.write(treeString.getBytes());
		
//            NewickTreeRerooter rerooter = new NewickTreeRerooter();
//            TreeNode root = rerooter.parseNewick(treeString);
//            TreeNode newRoot = rerooter.reroot(root, "CX280");
//            String newNewick = rerooter.toNewick(newRoot);
//        	zos.write(newNewick.getBytes());
        
//        	zos.write(NewickRooter.rootTree(treeString, /* "IRIS_313-10729" */ "CX280").getBytes());
//    		StringBuffer rootedTreeSB = new StringBuffer();
//    		rootTree(treeString, /* "IRIS_313-10729" */ "CX280").buildTreeString(rootedTreeSB);
//        	zos.write(rootedTreeSB.toString().getBytes());
        zos.closeEntry();
	}

	/* (non-Javadoc)
	 * @see fr.cirad.mgdb.exporting.IExportHandler#getStepList()
     */
    @Override
    public List<String> getStepList() {
        return new ArrayList<>(super.getStepList()) {{ add("Building Jukes-Cantor Newick tree"); }};
    }

	@Override
	public String[] getExportDataFileExtensions() {
		return new String[] {"nwk"};
	}
}