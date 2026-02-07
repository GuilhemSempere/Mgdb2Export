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
import java.util.List;
import java.util.Map;
import java.util.TreeSet;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;

import org.apache.log4j.Logger;
import org.bson.Document;
import org.ejml.data.FMatrixRMaj;
import org.springframework.data.mongodb.core.MongoTemplate;
import org.springframework.data.mongodb.core.query.Criteria;
import org.springframework.data.mongodb.core.query.Query;

import com.mongodb.BasicDBList;

import fr.cirad.mgdb.exporting.IExportHandler;
import fr.cirad.mgdb.exporting.markeroriented.EigenstratExportHandler;
import fr.cirad.mgdb.exporting.tools.ExportManager.ExportOutputs;
import fr.cirad.mgdb.exporting.tools.pca.PCACalculator;
import fr.cirad.mgdb.exporting.tools.pca.PCACalculator.FloatPCAResult;
import fr.cirad.mgdb.exporting.tools.pca.PCACalculator.OptimizedDiskFloatMatrix;
import fr.cirad.mgdb.model.mongo.maintypes.Assembly;
import fr.cirad.mgdb.model.mongo.subtypes.Callset;
import fr.cirad.tools.AlphaNumericComparator;
import fr.cirad.tools.ProgressIndicator;
import fr.cirad.tools.mgdb.VariantQueryWrapper;
import fr.cirad.tools.mongo.MongoTemplateManager;
import htsjdk.variant.variantcontext.VariantContext.Type;

/**
 * The Class PCAExportHandler.
 */
public class PCAExportHandler extends EigenstratExportHandler {

    /**
     * The Constant LOG.
     */
    private static final Logger LOG = Logger.getLogger(PCAExportHandler.class);

    private static final int numoutevec = 20;
    
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
    	return "Exports a zipped PCA based on an <a target='_blank' href='https://github.com/argriffing/eigensoft/blob/master/CONVERTF/README'>Eigenstrat</a> format export. Variants with more than 50% missing data are automatically excluded.";
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
    public void exportData(OutputStream outputStream, String sModule, Integer nAssemblyId, String sExportingUser, ProgressIndicator progress, String tmpVarCollName, VariantQueryWrapper varQueryWrapper, long markerCount, Map<String, String> markerSynonyms, Map<String, Collection<String>> individualsByPop, boolean workWithSamples, Map<String, HashMap<String, Float>> annotationFieldThresholds, Collection<Callset> callSetsToExport, Collection<String> individualMetadataFieldsToExport, Map<String, InputStream> readyToExportFiles) throws Exception {
    	TreeSet<String> indSet = new TreeSet<>(new AlphaNumericComparator<String>());
		for (Callset cs : callSetsToExport)
			indSet.add(workWithSamples ? cs.getSampleId() : cs.getIndividual());
		List<String> sortedIndividuals = new ArrayList<>(indSet);

		MongoTemplate mongoTemplate = MongoTemplateManager.get(sModule);
        ZipOutputStream zos = IExportHandler.createArchiveOutputStream(outputStream, readyToExportFiles, null);
		Assembly assembly = mongoTemplate.findOne(new Query(Criteria.where("_id").is(nAssemblyId)), Assembly.class);
        
		Collection<BasicDBList> variantDataQueries = varQueryWrapper.getVariantDataQueries();
        Document variantQueryForTargetCollection = variantDataQueries.isEmpty() ? new Document() : (tmpVarCollName == null ? new Document("$and", variantDataQueries.iterator().next()) : (varQueryWrapper.getBareQueries().iterator().hasNext() ? new Document("$and", varQueryWrapper.getBareQueries().iterator().next()) : new Document()));

        Map<String, Integer> individualPositions = IExportHandler.buildIndividualPositions(callSetsToExport, workWithSamples);

        if (progress.isAborted())
            return;

        File tempEigenstratFile = File.createTempFile("eigenstrat_for_pca__", "");
        File warningFile = File.createTempFile("export_warnings_", ""); // save existing warnings into a temp file so we can append to it
        try {
    		OutputStream eigenstratGenoOS = new FileOutputStream(tempEigenstratFile);   	
    	    ExportOutputs exportOutputs = writeGenotypeFile(eigenstratGenoOS, sModule, nAssemblyId, workWithSamples, annotationFieldThresholds, progress, tmpVarCollName, variantQueryForTargetCollection, markerCount, markerSynonyms, callSetsToExport, individualsByPop);
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
	        
//	    	long before = System.currentTimeMillis();

            progress.moveToNextStep();
            
			String exportName = sModule + (assembly != null && assembly.getName() != null ? "__" + assembly.getName() : "") + "__" + markerCount + "variants__" + sortedIndividuals.size() + "individuals";
	        if (individualMetadataFieldsToExport == null || !individualMetadataFieldsToExport.isEmpty())
	        	IExportHandler.addMetadataEntryIfAny(sModule + "__" + indSet.size() + (workWithSamples ? "sample" : "individual" ) + "s_metadata.tsv", sModule, sExportingUser, indSet, individualMetadataFieldsToExport, zos, (workWithSamples ? "sample" : "individual"), workWithSamples);
 
	        PCACalculator pcaCalc = new PCACalculator();            
            OptimizedDiskFloatMatrix diskMatrix = pcaCalc.readAndTransposeDirectToDisk(new BufferedInputStream(new FileInputStream(tempEigenstratFile)), warningOS, progress, (int) markerCount);

            if (progress.isAborted() || progress.getError() != null) {
                diskMatrix.close();
                return;
            }

            progress.moveToNextStep();
            pcaCalc.imputeMissingDataOnDisk(diskMatrix);
            if (progress.isAborted() || progress.getError() != null) {
                diskMatrix.close();
                return;
            }

            progress.moveToNextStep();
            pcaCalc.centerAndScaleDataOnDisk(diskMatrix);
            if (progress.isAborted() || progress.getError() != null) {
                diskMatrix.close();
                return;
            }

            ExecutorService executor = Executors.newSingleThreadExecutor();
            Future<FloatPCAResult> future = executor.submit(() -> {
                return pcaCalc.performSmartSVD(diskMatrix, null, warningOS, progress);
            });
            while (!future.isDone()) {
                if (progress.isAborted()) {
                    future.cancel(true);
                    LOG.debug("PCA aborted by user");
                    break;
                }
                Thread.sleep(500);
            }
            
            if (progress.isAborted() || progress.getError() != null) {
                diskMatrix.close();
                return;
            }

            FloatPCAResult floatPcaResult = future.get();
            FMatrixRMaj floatEigenVectors = floatPcaResult.getEigenVectors();

            float[][] transformedData = pcaCalc.transformDataBlockwise(diskMatrix, floatEigenVectors, numoutevec < individualPositions.size() ? numoutevec : null);
            float[] eigenValues = floatPcaResult.getEigenValues();
            
//	        LOG.debug("pca calculation took " + (System.currentTimeMillis() - before) / 1000d + "s");

            // Write results with memory efficiency
            zos.putNextEntry(new ZipEntry(exportName + "." + getExportDataFileExtensions()[0]));

            // Write eigenvalues header
            StringBuilder eigenHeader = new StringBuilder("#eigvals:");
            for (int c = 0; c < eigenValues.length; c++) {
                eigenHeader.append("\t").append(eigenValues[c]);
            }
            zos.write(eigenHeader.toString().getBytes());

            // Write data rows efficiently
            for (int r = 0; r < sortedIndividuals.size(); r++) {
                StringBuilder rowBuilder = new StringBuilder();
                rowBuilder.append("\n").append(sortedIndividuals.get(r));
                
                for (int c = 0; c < transformedData[r].length; c++)
                    rowBuilder.append("\t").append(transformedData[r][c]);
                
                zos.write(rowBuilder.toString().getBytes());
                
                // Help GC for very large datasets
                if (sortedIndividuals.size() > 10000 && r % 5000 == 0) {
                    transformedData[r] = null;
                    System.gc();
                }
            }

            zos.closeEntry();
            diskMatrix.close();
            transformedData = null;
            System.gc();

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
        return Arrays.asList(new String[]{"Extracting data to Eigenstrat format", "Transposing genotype matrix", "Imputing missing values", "Centering and scaling data"});
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