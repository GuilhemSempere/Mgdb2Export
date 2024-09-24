package fr.cirad.mgdb.exporting.tools.pca;

import java.io.IOException;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Scanner;

import org.apache.log4j.Logger;

import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.colt.matrix.tdouble.algo.decomposition.DenseDoubleSingularValueDecomposition;
import cern.colt.matrix.tdouble.impl.DenseDoubleMatrix2D;

public class ParallelPCACalculator {
	
	static private final Logger LOG = Logger.getLogger(ParallelPCACalculator.class);
	
	private int nMaxMissingDataPercentageForIndividuals = 50;
	
    public ParallelPCACalculator() {
    }
    
    public ParallelPCACalculator(int nMaxMissingDataPercentageForIndividuals) {
    	super();
    	this.nMaxMissingDataPercentageForIndividuals = nMaxMissingDataPercentageForIndividuals;
    }

    public double[][] readAndTransposeEigenstratGenoString(String fileContents, OutputStream warningOS) throws IOException {
        List<String> lines = new ArrayList<>();
        Scanner scanner = new Scanner(fileContents);
        while (scanner.hasNextLine()) {
            lines.add(scanner.nextLine());
        }
        scanner.close();

        if (lines.isEmpty()) {
            return new double[0][0];
        }

        int variantCount = lines.size();
        int individualCount = lines.get(0).length();
        List<Double>[] transposedData = new List[individualCount];
        ArrayList<Integer> variantsToIgnore = new ArrayList<>();

        for (int row = 0; row < variantCount; row++) {
            String line = lines.get(row);
            int nMissingDataCount = 0;
            for (int col = 0; col < individualCount; col++) {
            	if (row == 0)
            		transposedData[col] = new ArrayList<>();
                double value = Character.getNumericValue(line.charAt(col));
                if (value == 9) {
                    transposedData[col].add(Double.NaN);  // Handle missing data with NaN
                    nMissingDataCount++;
                } else
                	transposedData[col].add(value);
            }
            if (nMissingDataCount * 100 / line.length() > nMaxMissingDataPercentageForIndividuals) {
            	warningOS.write(("- Excluding variant #" + row + " from PCA export, it has too much missing data: " + (nMissingDataCount * 100 / line.length()) + "%\n").getBytes());
            	variantsToIgnore.add(row);
            }
        }
        
        if (!variantsToIgnore.isEmpty())
        	for (List<Double> indLine : transposedData)
        		for (int i = variantsToIgnore.size() - 1; i >= 0; i--)
        			indLine.remove((int) variantsToIgnore.get(i));

        return Arrays.stream(transposedData).map(list -> list.stream().mapToDouble(Double::doubleValue).toArray()).toArray(double[][]::new);
    }

    public double[][] imputeMissingValues(double[][] data) {
        int numRows = data.length;
        int numCols = data[0].length;

        double[] means = new double[numCols];

        // Calculate column means
        for (int col = 0; col < numCols; col++) {
            double sum = 0;
            int validCount = 0;
            for (int row = 0; row < numRows; row++) {
                if (!Double.isNaN(data[row][col])) {
                    sum += data[row][col];
                    validCount++;
                }
            }
            means[col] = sum / validCount;
        }

        // Impute missing values with column means
        double[][] imputedData = new double[numRows][numCols];
        for (int row = 0; row < numRows; row++) {
            for (int col = 0; col < numCols; col++) {
                if (Double.isNaN(data[row][col])) {
                    imputedData[row][col] = means[col];
                } else {
                    imputedData[row][col] = data[row][col];
                }
            }
        }

        return imputedData;
    }

    public double[][] centerAndScaleData(double[][] data) {
        int numRows = data.length;
        int numCols = data[0].length;

        double[] means = new double[numCols];
        double[] stdDevs = new double[numCols];

        // Calculate means
        for (int col = 0; col < numCols; col++) {
            double sum = 0;
            for (int row = 0; row < numRows; row++) {
                sum += data[row][col];
            }
            means[col] = sum / numRows;
        }

        // Calculate standard deviations
        for (int col = 0; col < numCols; col++) {
            double sum = 0;
            for (int row = 0; row < numRows; row++) {
                sum += Math.pow(data[row][col] - means[col], 2);
            }
            stdDevs[col] = sum == 0 ? 1e-20 : Math.sqrt(sum / (numRows - 1));
        }

        // Center and scale the data
        double[][] centeredData = new double[numRows][numCols];
        for (int row = 0; row < numRows; row++) {
            for (int col = 0; col < numCols; col++) {
                centeredData[row][col] = (data[row][col] - means[col]) / stdDevs[col];
            }
        }

        return centeredData;
    }

    public PCAResult performPCA(double[][] data) {
        if (hasInvalidValues(data)) {
        	LOG.warn("Matrix contains NaN or Infinite values.");
        }
        else
        	LOG.info("Matrix contains no NaN or Infinite values.");
        
        DoubleMatrix2D dataMatrix = new DenseDoubleMatrix2D(data);
        try {
	        DenseDoubleSingularValueDecomposition svd = new DenseDoubleSingularValueDecomposition(dataMatrix, true, false);
	        DoubleMatrix2D eigenVectors = svd.getV();
	        double[] singularValues = svd.getSingularValues();
	
	        // Calculate eigenvalues from singular values
	        int numSamples = data.length;
	        double[] eigenValues = new double[singularValues.length];
	        for (int i = 0; i < singularValues.length; i++) {
	            eigenValues[i] = (singularValues[i] * singularValues[i]) / (numSamples - 1);
	        }
	
	        return new PCAResult(eigenVectors, eigenValues, singularValues);
        } catch (IllegalArgumentException e) {
            System.err.println("Error occurred while computing SVD decomposition:");
            e.printStackTrace();
            printMatrixDimensions(data);
            if (hasInvalidValues(data)) {
                System.err.println("Matrix contains NaN or Infinite values.");
            }
            throw e;
        }
    }

    public void printMatrixDimensions(double[][] data) {
        int numRows = data.length;
        int numCols = data[0].length;
        System.out.printf("Matrix dimensions: %d rows x %d columns%n", numRows, numCols);
    }

    public boolean hasInvalidValues(double[][] data) {
//        for (int i=0; i<data.length; i++) {
//        	for (int j=0; i<data[i].length; j++) {
//                if (Double.isNaN(data[i][j]) || Double.isInfinite(data[i][j])) {
//                    return true;
//                }
//            }
//        }
        return false;
    }
    
    public double[][] transformData(double[][] data, DoubleMatrix2D pcaMatrix, Integer numEigenvectors) {
        DoubleMatrix2D dataMatrix = new DenseDoubleMatrix2D(data);
        DoubleMatrix2D t = numEigenvectors == null ? pcaMatrix : pcaMatrix.viewPart(0, 0, data[0].length, numEigenvectors);
        DoubleMatrix2D transformedData = dataMatrix.zMult(t, null);
        return transformedData.toArray();
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
}
