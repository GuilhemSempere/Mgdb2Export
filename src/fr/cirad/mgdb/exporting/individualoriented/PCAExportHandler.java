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

import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.InputStream;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;

import org.apache.log4j.Logger;
import org.ejml.data.FMatrixRMaj;
import org.springframework.data.mongodb.core.MongoTemplate;
import org.springframework.data.mongodb.core.query.Criteria;
import org.springframework.data.mongodb.core.query.Query;

import com.mongodb.BasicDBList;

import fr.cirad.mgdb.exporting.IExportHandler;
import fr.cirad.mgdb.exporting.markeroriented.EigenstratExportHandler;
import fr.cirad.mgdb.exporting.tools.ExportManager.ExportOutputs;
import fr.cirad.mgdb.exporting.tools.pca.ParallelPCACalculator;
import fr.cirad.mgdb.model.mongo.maintypes.Assembly;
import fr.cirad.mgdb.model.mongo.maintypes.GenotypingSample;
import fr.cirad.tools.AlphaNumericComparator;
import fr.cirad.tools.ExperimentalFeature;
import fr.cirad.tools.ProgressIndicator;
import fr.cirad.tools.mgdb.VariantQueryWrapper;
import fr.cirad.tools.mongo.MongoTemplateManager;
import htsjdk.variant.variantcontext.VariantContext.Type;

/**
 * The Class PCAExportHandler.
 */
public class PCAExportHandler extends EigenstratExportHandler implements ExperimentalFeature {

    /**
     * The Constant LOG.
     */
    private static final Logger LOG = Logger.getLogger(PCAExportHandler.class);

    private static final int numoutevec = 30;

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
        return "PCA";
    }

    /* (non-Javadoc)
	 * @see fr.cirad.mgdb.exporting.IExportHandler#getExportFormatDescription()
     */
    @Override
    public String getExportFormatDescription() {
    	return "Exports a zipped PCA based on an Eigenstrat format export. Variants with more than 50% missing data are automatically excluded.";
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
    public void exportData(OutputStream outputStream, String sModule, Integer nAssemblyId, String sExportingUser, ProgressIndicator progress, String tmpVarCollName, VariantQueryWrapper varQueryWrapper, long markerCount, Map<String, String> markerSynonyms, Map<String, Collection<String>> individuals, Map<String, HashMap<String, Float>> annotationFieldThresholds, Collection<GenotypingSample> samplesToExport, Collection<String> individualMetadataFieldsToExport, Map<String, InputStream> readyToExportFiles) throws Exception {
		List<String> sortedIndividuals = samplesToExport.stream().map(gs -> gs.getIndividual()).distinct().sorted(new AlphaNumericComparator<String>()).collect(Collectors.toList());

	    // Check if there is enough memory
//		System.gc();
//	    Runtime runtime = Runtime.getRuntime();
//	    long availableRAM = runtime.maxMemory() - (runtime.totalMemory() - runtime.freeMemory());
//	    long matrixSize = markerCount * sortedIndividuals.size();
//	    if (matrixSize * 3 > availableRAM * 0.01)	// Empirically, we found that if matrix sizer is < 1% of the available memory, SVD calculation is possible. But we account for concurrency by making sure we have at least 3x more RAM available
//	        throw new Exception("Not enough memory to process so much data. Please reduce matrix size to under " + (int) (availableRAM * 0.01 / 3) + " genotypes!");
		
		MongoTemplate mongoTemplate = MongoTemplateManager.get(sModule);
        ZipOutputStream zos = IExportHandler.createArchiveOutputStream(outputStream, readyToExportFiles);
		Assembly assembly = mongoTemplate.findOne(new Query(Criteria.where("_id").is(nAssemblyId)), Assembly.class);
        

//      FIXME: remove me (only for testing)
        if (individualMetadataFieldsToExport == null || !individualMetadataFieldsToExport.isEmpty()) {
        	if (!IExportHandler.addMetadataEntryIfAny(sModule + "__" + sortedIndividuals.size() + "individuals_metadata.tsv", sModule, sExportingUser, sortedIndividuals, individualMetadataFieldsToExport, zos, "individual")) {
		    	zos.putNextEntry(new ZipEntry(sModule + "__" + sortedIndividuals.size() + "individuals_metadata.tsv"));
		    	zos.write("individual\tpopulation".getBytes());
		    	for (String ind : sortedIndividuals) {
		    		zos.write(("\n" + ind + "\t" + ind.split("_")[0]).getBytes());
		    	}
		    	zos.closeEntry();
	        }
        }
        
        Collection<BasicDBList> variantDataQueries = varQueryWrapper.getVariantDataQueries();

		Map<String, Integer> individualPositions = new LinkedHashMap<>();
		for (String ind : samplesToExport.stream().map(gs -> gs.getIndividual()).distinct().sorted(new AlphaNumericComparator<String>()).collect(Collectors.toList()))
			individualPositions.put(ind, individualPositions.size());

        if (progress.isAborted())
            return;

        File tempEigenstratFile = File.createTempFile("eigenstrat_for_pca__", "");
        File warningFile = File.createTempFile("export_warnings_", ""); // save existing warnings into a temp file so we can append to it
        try {
    		OutputStream eigenstratGenoOS = new FileOutputStream(tempEigenstratFile);   	
            ExportOutputs exportOutputs = writeGenotypeFile(eigenstratGenoOS, sModule, nAssemblyId, individuals, annotationFieldThresholds, progress, tmpVarCollName, !variantDataQueries.isEmpty() ? variantDataQueries.iterator().next() : new BasicDBList(), markerCount, markerSynonyms, samplesToExport, individualPositions);
            eigenstratGenoOS.close();

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

            progress.moveToNextStep();
            ParallelPCACalculator pcaCalc = new ParallelPCACalculator();            
            float[][] data = pcaCalc.readAndTransposeEigenstratGenoString(new BufferedInputStream(new FileInputStream(tempEigenstratFile)), warningOS, progress, (int) markerCount);

	        if (progress.isAborted())
	            return;
	        
            progress.moveToNextStep();
            data = pcaCalc.imputeMissingValues(data);
	        if (progress.isAborted())
	            return;
	        
            progress.moveToNextStep();
            data = pcaCalc.centerAndScaleData(data);
	        if (progress.isAborted())
	            return;
	        
	        progress.moveToNextStep();

//	        // compute with doubles
//	        double[][] doubleData = new double[data.length][data.length == 0 ? 0 : data[0].length];
//	        for (int i = 0; i < data.length; i++)
//	            for (int j = 0; j < (data.length == 0 ? 0 : data[0].length); j++)
//	                doubleData[i][j] = (double) data[i][j];
//	        ParallelPCACalculator.PCAResult pcaResult = pcaCalc.performPCA(doubleData);
//	        DoubleMatrix2D eigenVectors = pcaResult.getEigenVectors();
//	        double[][] dataMatrix = pcaCalc.transformData(doubleData, eigenVectors, numoutevec < individualPositions.size() ? numoutevec : null);
//	        double[] eigenValues = pcaResult.getEigenValues();
	        
	        // compute with floats
	        ParallelPCACalculator.FloatPCAResult floatPcaResult = pcaCalc.performPCA_float_disk(data);
	        FMatrixRMaj floatEigenVectors = floatPcaResult.getEigenVectors();
	        float[][] dataMatrix = pcaCalc.transformData(data, floatEigenVectors, numoutevec < individualPositions.size() ? numoutevec : null);
	        float[] eigenValues = floatPcaResult.getEigenValues();
	        
			String exportName = sModule + (assembly != null && assembly.getName() != null ? "__" + assembly.getName() : "") + "__" + dataMatrix[0].length + "variants__" + sortedIndividuals.size() + "individuals";
	        zos.putNextEntry(new ZipEntry(exportName + "." + getExportDataFileExtensions()[0]));
	        zos.write("#eigvals:".getBytes());
	    	for (int c=0; c<eigenValues.length; c++)
	    		zos.write(("\t" + eigenValues[c]).getBytes());
	        for (int r=0; r<sortedIndividuals.size(); r++) {
	        	zos.write(("\n" + sortedIndividuals.get(r)).getBytes());
	        	for (int c=0; c<dataMatrix[r].length; c++)
	        		zos.write(("\t" + dataMatrix[r][c]).getBytes());
	        }

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
        }
        catch (OutOfMemoryError oome) {
	    	LOG.error("Not enough RAM to transpose matrix", oome);
        	progress.setError("Unable to compute PCA (dataset is too large): " + oome.getMessage());
        }
        finally {
            warningFile.delete();
            tempEigenstratFile.delete();
        }

        zos.finish();
        zos.close();
        progress.setPercentageEnabled(true);
        progress.setCurrentStepProgress((short) 100);
    }

	/* (non-Javadoc)
	 * @see fr.cirad.mgdb.exporting.IExportHandler#getStepList()
     */
    @Override
    public List<String> getStepList() {
        return Arrays.asList(new String[]{"Extracting data to Eigenstrat format", "Transposing genotype matrix", "Imputing missing values", "Centering and scaling data", "Computing PCA"});
    }

	@Override
	public String[] getExportDataFileExtensions() {
		return new String[] {"pca"};
	}

    @Override
    public int[] getSupportedPloidyLevels() {
        return new int[] {2};
    }
}