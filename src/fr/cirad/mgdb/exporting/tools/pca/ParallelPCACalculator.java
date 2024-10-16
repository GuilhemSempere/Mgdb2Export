package fr.cirad.mgdb.exporting.tools.pca;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.util.HashSet;
import java.util.Scanner;

import org.apache.log4j.Logger;
import org.ejml.data.FMatrixRBlock;
import org.ejml.data.FMatrixRMaj;
import org.ejml.dense.block.MatrixOps_FDRB;
import org.ejml.dense.row.CommonOps_FDRM;
import org.ejml.dense.row.decomposition.svd.SvdImplicitQrDecompose_FDRM;

import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.colt.matrix.tdouble.algo.decomposition.DenseDoubleSingularValueDecomposition;
import cern.colt.matrix.tdouble.impl.DenseDoubleMatrix2D;
import fr.cirad.tools.ProgressIndicator;

public class ParallelPCACalculator {
	
	static private final Logger LOG = Logger.getLogger(ParallelPCACalculator.class);
	
	private int nMaxMissingDataPercentageForIndividuals = 50;
	
    public ParallelPCACalculator() {
    }
    
    public ParallelPCACalculator(int nMaxMissingDataPercentageForIndividuals) {
    	super();
    	this.nMaxMissingDataPercentageForIndividuals = nMaxMissingDataPercentageForIndividuals;
    }

    public float[][] readAndTransposeEigenstratGenoString(InputStream eigenstratContents, OutputStream warningOS, ProgressIndicator progress, int nVariantCount) throws IOException {
	    try (Scanner scanner = new Scanner(eigenstratContents)) {
	        if (!scanner.hasNextLine())
	            return new float[0][0];
	
	        int row = 0;
	        float[][] transposedData = null;
	        HashSet<Integer> variantsToIgnore = new HashSet<>();
	        
	        while (scanner.hasNextLine()) {
	            String line = scanner.nextLine();
	            if (row == 0)
	            	transposedData = new float[line.length()][nVariantCount];
	
	            int nMissingDataCount = 0;
	            for (int col = 0; col < transposedData.length; col++) {
	                float value = Character.getNumericValue(line.charAt(col));
	                if (value == 9) {
	                    transposedData[col][row] = Float.NaN;  // Handle missing data with NaN
	                    nMissingDataCount++;
	                } else
	                	transposedData[col][row] = value;
	            }
	            if (nMissingDataCount * 100f / line.length() > nMaxMissingDataPercentageForIndividuals) {
	            	warningOS.write(("- Excluding variant #" + row + " from PCA export, it has too much missing data: " + (nMissingDataCount * 100 / line.length()) + "%\n").getBytes());
	            	variantsToIgnore.add(row);
	            }
	            row++;
	            progress.setCurrentStepProgress((int) (row * 100f / nVariantCount));
	        }
	        
	        if (variantsToIgnore.isEmpty())
	        	return transposedData;

        	float[][] finalResult = new float[transposedData.length][];
        	for (row = 0; row < transposedData.length; row++) {
        		finalResult[row] = new float[nVariantCount - variantsToIgnore.size()];	// allocate memory row by row to consume less
        		int j = 0;
        		for (int i = 0; i < transposedData[row].length; i++)
	        		if (!variantsToIgnore.contains(i))
	        			finalResult[row][j++] = transposedData[row][i];
        		transposedData[row] = null;	//free up memory
        	}
        	return finalResult;
	    }
    }

    public float[][] imputeMissingValues(float[][] data) {
        int numRows = data.length;
        int numCols = data[0].length;

        float[] means = new float[numCols];

        // Calculate column means
        for (int col = 0; col < numCols; col++) {
            float sum = 0;
            int validCount = 0;
            for (int row = 0; row < numRows; row++) {
                if (!Float.isNaN(data[row][col])) {
                    sum += data[row][col];
                    validCount++;
                }
            }
            means[col] = sum / validCount;
        }

        // Impute missing values with column means
        float[][] imputedData = new float[numRows][];
        for (int row = 0; row < numRows; row++) {
        	imputedData[row] = new float[numCols];	// allocate memory row by row to consume less
            for (int col = 0; col < numCols; col++) {
                if (Float.isNaN(data[row][col])) {
                    imputedData[row][col] = means[col];
                } else {
                    imputedData[row][col] = data[row][col];
                }
            }
            data[row] = null;	//free up memory
        }

        return imputedData;
    }

    public float[][] centerAndScaleData(float[][] data) {
        int numRows = data.length;
        int numCols = data[0].length;

        float[] means = new float[numCols];
        float[] stdDevs = new float[numCols];

        // Calculate means
        for (int col = 0; col < numCols; col++) {
            float sum = 0;
            for (int row = 0; row < numRows; row++) {
                sum += data[row][col];
            }
            means[col] = sum / numRows;
        }

        // Calculate standard deviations
        for (int col = 0; col < numCols; col++) {
            float sum = 0;
            for (int row = 0; row < numRows; row++) {
                sum += Math.pow(data[row][col] - means[col], 2);
            }
            stdDevs[col] = sum == 0 ? 1e-20f : (float) Math.sqrt(sum / (numRows - 1));
        }

        // Center and scale the data
        float[][] centeredData = new float[numRows][];
        for (int row = 0; row < numRows; row++) {
        	centeredData[row] = new float[numCols];	// allocate memory row by row to consume less
            for (int col = 0; col < numCols; col++)
                centeredData[row][col] = (data[row][col] - means[col]) / stdDevs[col];
            data[row] = null;	//free up memory
        }

        return centeredData;
    }

    public FloatPCAResult performPCA_float_ram(float[][] data) {
    	LOG.debug("Performing in-memory PCA using floats");
    	int blockSize = 45000;
    	
        if (hasInvalidValues(data))
        	LOG.warn("Matrix contains NaN or Infinite values.");
        else
        	LOG.info("Matrix contains no NaN or Infinite values.");

        // Step 1: Convert the input data to a block matrix (FMatrixRBlock)
        FMatrixRBlock matrixBlock = new FMatrixRBlock(data.length, data[0].length, blockSize);
        for (int i = 0; i < data.length; i++) {
            for (int j = 0; j < data[0].length; j++) {
                matrixBlock.set(i, j, data[i][j]);
            }
        }

        // Step 2: Convert FMatrixRBlock to FMatrixRMaj
        FMatrixRMaj matrix = new FMatrixRMaj(matrixBlock.numRows, matrixBlock.numCols);
        MatrixOps_FDRB.convert(matrixBlock, matrix);

        // Step 3: Perform SVD on the FMatrixRMaj using SvdImplicitQrDecompose_FDRM
        SvdImplicitQrDecompose_FDRM svd = new SvdImplicitQrDecompose_FDRM(true, false, true, true);
        if (!svd.decompose(matrix)) {
            throw new RuntimeException("SVD failed to decompose the matrix.");
        }

        // Step 4: Get the V matrix (right singular vectors, i.e., eigenvectors for PCA)
        FMatrixRMaj V = svd.getV(null, false);

        // Step 5: Get the singular values
        float[] singularValues = svd.getSingularValues();

        // Step 6: Calculate eigenvalues from singular values
        int numSamples = matrix.numRows;
        float[] eigenValues = new float[singularValues.length];
        for (int i = 0; i < singularValues.length; i++) {
            eigenValues[i] = (singularValues[i] * singularValues[i]) / (numSamples - 1);
        }

        // Step 7: Return the PCA results
        return new FloatPCAResult(V, eigenValues, singularValues);
    }

    public FloatPCAResult performPCA_float_disk(float[][] data) throws IOException {
    	LOG.debug("Performing memory-mapped-file-based PCA using floats");
    	int blockSize = 45000;
    	
        if (hasInvalidValues(data))
        	LOG.warn("Matrix contains NaN or Infinite values.");
        else
        	LOG.info("Matrix contains no NaN or Infinite values.");

        // Create a memory-mapped file for the data matrix
        File matrixFile = File.createTempFile("FileBackedDenseColumnFloatMatrix2D_", "");
        FileBackedDenseColumnFloatMatrix2D fileDataMatrix = new FileBackedDenseColumnFloatMatrix2D(data, matrixFile);

        int numSamples = data.length;
        data = null;  // free up memory
        
        try {
            // Convert the memory-mapped matrix to an FMatrixRBlock for large-scale SVD
            FMatrixRBlock blockMatrix = new FMatrixRBlock(fileDataMatrix.numRows(), fileDataMatrix.numCols(), blockSize);
            for (int i = 0; i < fileDataMatrix.numRows(); i++)
                for (int j = 0; j < fileDataMatrix.numCols(); j++)
                    blockMatrix.set(i, j, fileDataMatrix.get(i, j));  // Reading from memory-mapped file

            // Perform Singular Value Decomposition (SVD) using block decomposition
            FMatrixRMaj matrix = new FMatrixRMaj(blockMatrix.numRows, blockMatrix.numCols);
            MatrixOps_FDRB.convert(blockMatrix, matrix);
            blockMatrix = null; 	// free up memory
            SvdImplicitQrDecompose_FDRM svd = new SvdImplicitQrDecompose_FDRM(true, false, true, true);
            if (!svd.decompose(matrix)) {
                throw new RuntimeException("SVD failed to decompose the matrix.");
            }
            
            // Get the singular values
            float[] singularValues = svd.getSingularValues();

            // Calculate eigenvalues from singular values
            float[] eigenValues = new float[singularValues.length];
            for (int i = 0; i < singularValues.length; i++)
                eigenValues[i] = (singularValues[i] * singularValues[i]) / (numSamples - 1);

            return new FloatPCAResult(svd.getV(null, false), eigenValues, singularValues);
        }
        catch (OutOfMemoryError oome) {
        	LOG.error("Not enough RAM to compute PCA", oome);
        	throw oome;
        }
    	finally {
            // Clean up the memory-mapped file
            matrixFile.delete();
        }
    }

    private boolean hasInvalidValues(float[][] data) {
        for (float[] row : data) {
            for (float val : row) {
                if (Float.isNaN(val) || Float.isInfinite(val)) {
                    return true;
                }
            }
        }
        return false;
    }

    public PCAResult performPCA_float_double(double[][] data) throws IOException {
    	LOG.debug("Performing memory-mapped-file-based PCA using doubles");
    	
        if (hasInvalidValues(data))
        	LOG.warn("Matrix contains NaN or Infinite values.");
        else
        	LOG.info("Matrix contains no NaN or Infinite values.");
        
        int numSamples = data.length;
        
//      DenseColumnDoubleMatrix2D dataMatrix = new DenseColumnDoubleMatrix2D(data);
        File matrixFile = File.createTempFile("FileBackedDenseColumnDoubleMatrix2D_", "");
        FileBackedDenseColumnDoubleMatrix2D dataMatrix = new FileBackedDenseColumnDoubleMatrix2D(data, matrixFile);

        data = null;	// freeing memory (won't be able to call hasInvalidValues after that...
        try {
	        DenseDoubleSingularValueDecomposition svd = new DenseDoubleSingularValueDecomposition(dataMatrix, true, false);
	        double[] singularValues = svd.getSingularValues();
	        DoubleMatrix2D eigenVectors = svd.getV();	// FIXME: this causes Java heap space on big datasets!!
	
	        // Calculate eigenvalues from singular values
	        double[] eigenValues = new double[singularValues.length];
	        for (int i = 0; i < singularValues.length; i++) {
	            eigenValues[i] = (singularValues[i] * singularValues[i]) / (numSamples - 1);
	        }
	
	        return new PCAResult(eigenVectors, eigenValues, singularValues);
        }
        finally {
        	matrixFile.delete();
        }
    }
    
    public boolean hasInvalidValues(double[][] data) {
        for (int i=0; i<data.length; i++) {
        	for (int j=0; j<data[i].length; j++) {
                if (Double.isNaN(data[i][j]) || Double.isInfinite(data[i][j])) {
                    return true;
                }
            }
        }
        return false;
    }
    
    public double[][] transformData(double[][] data, DoubleMatrix2D pcaMatrix, Integer numEigenvectors) {
        DoubleMatrix2D dataMatrix = new DenseDoubleMatrix2D(data);
        DoubleMatrix2D t = numEigenvectors == null ? pcaMatrix : pcaMatrix.viewPart(0, 0, data[0].length, numEigenvectors);
        DoubleMatrix2D transformedData = dataMatrix.zMult(t, null);
        return transformedData.toArray();
    }
    
    public float[][] transformData(float[][] data, FMatrixRMaj pcaMatrix, Integer numEigenvectors) {
        // Step 1: Convert double[][] data to float[][] and then to FMatrixRMaj
        int rows = data.length;
        int cols = data[0].length;
        FMatrixRMaj dataMatrix = new FMatrixRMaj(rows, cols);
        
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                dataMatrix.set(i, j, (float) data[i][j]);
            }
        }

        // Step 2: Extract part of the PCA matrix if numEigenvectors is specified
        FMatrixRMaj tMatrix;
        if (numEigenvectors != null) {
            tMatrix = new FMatrixRMaj(cols, numEigenvectors);
            CommonOps_FDRM.extract(pcaMatrix, 0, cols, 0, numEigenvectors, tMatrix, 0, 0);
        } else {
            tMatrix = pcaMatrix;
        }

        // Step 3: Perform matrix multiplication (dataMatrix * tMatrix)
        FMatrixRMaj transformedData = new FMatrixRMaj(rows, tMatrix.numCols);
        CommonOps_FDRM.mult(dataMatrix, tMatrix, transformedData);

        // Step 4: Convert the transformed data (FMatrixRMaj) back to float[][]
        float[][] result = new float[transformedData.numRows][transformedData.numCols];
        for (int i = 0; i < transformedData.numRows; i++) {
            for (int j = 0; j < transformedData.numCols; j++) {
                result[i][j] = transformedData.get(i, j);
            }
        }

        // Step 5: Return the transformed data
        return result;
    }

    public class PCAResult {
        private DoubleMatrix2D eigenVectors;
        private double[] eigenValues;
        private double[] singularValues;

        public PCAResult(DoubleMatrix2D eigenVectors, double[] eigenValues, double[] singularValues) {
            this.eigenVectors = eigenVectors;
            this.eigenValues = eigenValues;
            this.singularValues = singularValues;
        }

        public DoubleMatrix2D getEigenVectors() {
            return eigenVectors;
        }
        
        public double[] getEigenValues() {
            return eigenValues;
        }

        public double[] getSingularValues() {
            return singularValues;
        }
    }
    
    public static class FloatPCAResult {
        public final FMatrixRMaj eigenVectors;
        public final float[] eigenValues;
        public final float[] singularValues;

        public FloatPCAResult(FMatrixRMaj eigenVectors, float[] eigenValues, float[] singularValues) {
            this.eigenVectors = eigenVectors;
            this.eigenValues = eigenValues;
            this.singularValues = singularValues;
        }
        
        public FMatrixRMaj getEigenVectors() {
            return eigenVectors;
        }
        
        public float[] getEigenValues() {
            return eigenValues;
        }
    }
}
