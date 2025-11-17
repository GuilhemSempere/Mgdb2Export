package fr.cirad.mgdb.exporting.tools.pca;

import java.io.BufferedReader;
import java.io.Closeable;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.RandomAccessFile;
import java.nio.FloatBuffer;
import java.nio.MappedByteBuffer;
import java.nio.channels.FileChannel;
import java.nio.file.Files;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Random;
//import java.util.concurrent.Callable;
//import java.util.concurrent.ExecutorService;
//import java.util.concurrent.Executors;
//import java.util.concurrent.Future;
//import java.util.concurrent.ScheduledExecutorService;
//import java.util.concurrent.TimeUnit;

import org.apache.log4j.Logger;
import org.ejml.data.FMatrixRMaj;
import org.ejml.dense.row.CommonOps_FDRM;
import org.ejml.dense.row.decomposition.svd.SvdImplicitQrDecompose_FDRM;

import fr.cirad.tools.ProgressIndicator;

public class PCACalculator {
    
    static private final Logger LOG = Logger.getLogger(PCACalculator.class);
    
    private int nMaxMissingDataPercentageForIndividuals = 50;

    public PCACalculator() {}
    
    public PCACalculator(int nMaxMissingDataPercentageForIndividuals) {
        this.nMaxMissingDataPercentageForIndividuals = nMaxMissingDataPercentageForIndividuals;
    }

    /**
     * Read Eigenstrat data and stream directly to disk matrix without loading into memory
     */
    public OptimizedDiskFloatMatrix readAndTransposeDirectToDisk(
            InputStream eigenstratContents, OutputStream warningOS, 
            ProgressIndicator progress, int nVariantCount) throws IOException {
        
        try (BufferedReader reader = new BufferedReader(new InputStreamReader(eigenstratContents))) {
            String firstLine = reader.readLine();
            if (firstLine == null)
                return new OptimizedDiskFloatMatrix(0, 0);

            int sampleCount = firstLine.length();

            // Create disk matrix with correct dimensions: samples × variants
            OptimizedDiskFloatMatrix diskMatrix = new OptimizedDiskFloatMatrix(sampleCount, nVariantCount);
            HashSet<Integer> variantsToIgnore = new HashSet<>();

            // Process first line
            processVariantLine(firstLine, 0, diskMatrix, variantsToIgnore, warningOS);

            // Process remaining lines
            int variantIndex = 1;
            String line;
            while ((line = reader.readLine()) != null) {
                processVariantLine(line, variantIndex, diskMatrix, variantsToIgnore, warningOS);
                variantIndex++;

                if (progress != null && nVariantCount > 0) {
                    progress.setCurrentStepProgress((int) (variantIndex * 100f / nVariantCount));
                }
            }

            // Handle variants to ignore by creating a new matrix without them
            if (!variantsToIgnore.isEmpty()) {
                return createFilteredDiskMatrix(diskMatrix, variantsToIgnore, sampleCount, nVariantCount);
            }

            return diskMatrix;
        }
    }

    private void processVariantLine(String line, int variantIndex, OptimizedDiskFloatMatrix diskMatrix, HashSet<Integer> variantsToIgnore, OutputStream warningOS) throws IOException {
        int nMissingDataCount = 0;
        for (int sampleIndex = 0; sampleIndex < line.length(); sampleIndex++) {
            float value = Character.getNumericValue(line.charAt(sampleIndex));
            if (value == 9) {
                diskMatrix.set(sampleIndex, variantIndex, Float.NaN);
                nMissingDataCount++;
            } else {
                diskMatrix.set(sampleIndex, variantIndex, value);
            }
        }
        
        if (nMissingDataCount * 100f / line.length() > nMaxMissingDataPercentageForIndividuals) {
            warningOS.write(("- Excluding variant #" + variantIndex + " from PCA export, it has too much missing data: " + 
                           (nMissingDataCount * 100 / line.length()) + "%\n").getBytes());
            variantsToIgnore.add(variantIndex);
        }
    }

    private OptimizedDiskFloatMatrix createFilteredDiskMatrix(OptimizedDiskFloatMatrix original,
                                                            HashSet<Integer> variantsToIgnore,
                                                            int sampleCount, int originalVariantCount) 
                                                            throws IOException {
        int filteredVariantCount = originalVariantCount - variantsToIgnore.size();
        OptimizedDiskFloatMatrix filteredMatrix = new OptimizedDiskFloatMatrix(sampleCount, filteredVariantCount);
        
        for (int sample = 0; sample < sampleCount; sample++) {
            int targetVariant = 0;
            for (int sourceVariant = 0; sourceVariant < originalVariantCount; sourceVariant++) {
                if (!variantsToIgnore.contains(sourceVariant)) {
                    filteredMatrix.set(sample, targetVariant, original.get(sample, sourceVariant));
                    targetVariant++;
                }
            }
        }
        
        original.close(); // Clean up original matrix
        return filteredMatrix;
    }

    /**
     * Impute missing data directly on disk matrix
     * @throws InterruptedException 
     */
    public void imputeMissingDataOnDisk(OptimizedDiskFloatMatrix diskMatrix) throws InterruptedException {
        int numRows = diskMatrix.rows();
        int numCols = diskMatrix.cols();
        float[] means = new float[numCols];
        int[] counts = new int[numCols];

        // Reusable row buffer to read rows contiguously
        float[] rowBuf = new float[numCols];

        // Pass 1: compute sums & counts (row-major => good locality)
        for (int r = 0; r < numRows; r++) {
            diskMatrix.getRow(r, rowBuf); // fills rowBuf[0..numCols-1]
            for (int c = 0; c < numCols; c++) {
                float v = rowBuf[c];
                if (!Float.isNaN(v)) {
                    means[c] += v;
                    counts[c]++;
                }
            }
        }

        // Finalize means
        for (int c = 0; c < numCols; c++) {
            means[c] = counts[c] > 0 ? means[c] / counts[c] : 0f;
        }

        // Pass 2: impute missing values and write back rows
        for (int r = 0; r < numRows; r++) {
            diskMatrix.getRow(r, rowBuf);
            boolean modified = false;
            for (int c = 0; c < numCols; c++) {
                if (Float.isNaN(rowBuf[c])) {
                    rowBuf[c] = means[c];
                    modified = true;
                }
            }
            if (modified) {
                diskMatrix.setRow(r, rowBuf);
            }
        }
    }

    /**
     * Parallel center and scale operations
     */
    public void centerAndScaleDataOnDisk(OptimizedDiskFloatMatrix diskMatrix) throws InterruptedException {
        int numRows = diskMatrix.rows();
        int numCols = diskMatrix.cols();

        float[] means = new float[numCols];
        float[] m2 = new float[numCols]; // for Welford's M2
        int[] counts = new int[numCols];

        float[] rowBuf = new float[numCols];

        // Pass 1: Welford accumulation column-wise by reading rows once
        for (int r = 0; r < numRows; r++) {
            diskMatrix.getRow(r, rowBuf);
            for (int c = 0; c < numCols; c++) {
                float x = rowBuf[c];
                // assuming missing values already imputed (if not, they will be included)
                // Standard Welford update
                counts[c]++;
                float delta = x - means[c];
                means[c] += delta / counts[c];
                float delta2 = x - means[c];
                m2[c] += delta * delta2;
            }
        }

        // Compute stdDevs
        float[] stdDevs = new float[numCols];
        for (int c = 0; c < numCols; c++) {
            if (counts[c] > 1) {
                float variance = m2[c] / (counts[c] - 1);
                stdDevs[c] = variance <= 0f ? 1e-20f : (float) Math.sqrt(variance);
            } else {
                stdDevs[c] = 1e-20f;
            }
        }

        // Pass 2: normalize rows and write back
        for (int r = 0; r < numRows; r++) {
            diskMatrix.getRow(r, rowBuf);
            boolean modified = false;
            for (int c = 0; c < numCols; c++) {
                float centeredScaled = (rowBuf[c] - means[c]) / stdDevs[c];
                // Avoid creating extra NaNs; replace in place
                rowBuf[c] = centeredScaled;
                modified = true;
            }
            if (modified) {
                diskMatrix.setRow(r, rowBuf);
            }
        }
    }

//    public Future<FloatPCAResult> performSmartSVDWithProgressAsync(
//            OptimizedDiskFloatMatrix diskMatrix,
//            Integer numComponents,
//            OutputStream warningOS,
//            ProgressIndicator progress,
//            ExecutorService executor) {
//
//        // Main SVD task (same math, blocking)
//        Callable<FloatPCAResult> svdTask = () -> {
//            try {
//                // Small delay so progress animation has time to start
//                Thread.sleep(150);
//            } catch (InterruptedException e) {}
//
//            return performSmartSVD(diskMatrix, numComponents, warningOS);
//        };
//
//        // Submit SVD job
//        Future<FloatPCAResult> svdFuture = executor.submit(svdTask);
//
//        // Start smooth progress animator
//        ScheduledExecutorService ticker = Executors.newSingleThreadScheduledExecutor();
//        final long startTime = System.currentTimeMillis();
//        final int ESTIMATED_MS = estimateSvdDurationMs(diskMatrix.rows(), diskMatrix.cols());
//
//        ticker.scheduleAtFixedRate(() -> {
//            if (svdFuture.isDone() || svdFuture.isCancelled()) {
//                if (progress != null) {
//                    progress.setCurrentStepProgress(100);
//                }
//                ticker.shutdown();
//                return;
//            }
//
//            if (progress != null) {
//                long elapsed = System.currentTimeMillis() - startTime;
//                float pct = Math.min(95f, (elapsed * 100f / ESTIMATED_MS)); // cap at 95%
//                progress.setCurrentStepProgress((int) pct);
//            }
//        }, 0, 600, TimeUnit.MILLISECONDS);
//
//        return svdFuture;
//    }
//
//    private int estimateSvdDurationMs(int n, int m) {
//        // Heuristic based on expected complexity O(min(n*m^2, m*n^2))
//        double cost = Math.min(
//                (double)n * m * m,
//                (double)m * n * n
//        );
//
//        // Normalize to a baseline (tune as needed)
//        // For example: cost=1e10 -> ~10 seconds
//        double BASE_COST = 1e10;
//        double estimatedSeconds = Math.max(2, cost / BASE_COST * 10.0);
//
//        return (int)(estimatedSeconds * 1000);
//    }


    /**
     * Smart SVD that tries exact method first, falls back to memory-safe methods on failure
     * @throws IOException 
     */
    public FloatPCAResult performSmartSVD(OptimizedDiskFloatMatrix diskMatrix, Integer numComponents, OutputStream warningOS) throws IOException {
        int n = diskMatrix.rows();
        int m = diskMatrix.cols();
        
        if (numComponents == null) {
            numComponents = Math.min(100, Math.min(n, m));
        }
                
        // Try exact SVD first (fastest when it works)
        try {
            LOG.debug("Attempting exact smart SVD for matrix " + n + " × " + m + ", components=" + numComponents);
            return performExactSVDWithMemoryCheck(diskMatrix, numComponents);
        } catch (OutOfMemoryError e) {
            LOG.warn("Exact SVD failed due to memory, falling back to randomized SVD");
            System.gc(); // Clean up after OOM
            
            // Fallback to randomized SVD
            return performRandomizedSVD(diskMatrix, numComponents, 10, warningOS);
        }
    }
    
    public FloatPCAResult performRandomizedSVD(OptimizedDiskFloatMatrix diskMatrix, int k, int p, OutputStream warningOS) throws IOException {
        int n = diskMatrix.rows();     // number of samples
        int m = diskMatrix.cols();     // number of SNPs
        
        if (k <= 0)
        	k = Math.min(100, Math.min(n, m));
        
        // Start with 50% of SNPs, then halve recursively if needed
        int initialTargetSNPs = Math.max(1000, m / 2); // Start with 50%, minimum 1000
        return performRandomizedSVDRecursive(diskMatrix, n, m, k, initialTargetSNPs, warningOS);
    }

	private FloatPCAResult performRandomizedSVDRecursive(OptimizedDiskFloatMatrix diskMatrix, int n, int m, int k, int targetSNPs, OutputStream warningOS) throws IOException {
        int actualTargetSNPs = Math.min(targetSNPs, m);
	    try {
	        int[] selectedIndices = selectRandomSNPs(m, actualTargetSNPs);

	        LOG.info("Performing simplified randomized SVD: " + n + " samples × " + selectedIndices.length + " SNPs (sampled from " + m + " total SNPs), k=" + k);
	        
	        FMatrixRMaj subsampled = extractSubsampledMatrix(diskMatrix, n, selectedIndices);
	        
	        // Perform exact SVD on subsampled matrix
	        SvdImplicitQrDecompose_FDRM svd = new SvdImplicitQrDecompose_FDRM(true, true, true, true);
	        if (!svd.decompose(subsampled))
	            throw new RuntimeException("SVD failed on subsampled matrix");
	        
	        FMatrixRMaj V_small = svd.getV(null, false);
	        float[] singularValues = svd.getSingularValues();
	        
	        // Create full-size eigenvector matrix with zeros for non-sampled SNPs
	        FMatrixRMaj fullEigenVectors = new FMatrixRMaj(m, k);
	        for (int comp = 0; comp < k; comp++)
	            for (int i = 0; i < selectedIndices.length; i++)
	                fullEigenVectors.set(selectedIndices[i], comp, V_small.get(i, comp));

	        float[] eigenValues = new float[k];
	        for (int i = 0; i < k; i++)
	            eigenValues[i] = (singularValues[i] * singularValues[i]) / (n - 1);
	        
	        // Write final warning only once with the actual SNPs used
	        if (actualTargetSNPs < m) {
	            double percentage = (actualTargetSNPs * 100.0) / m;
	            warningOS.write(("- PCA computed using randomized SVD on " + actualTargetSNPs + " of " + m + " SNPs (" + String.format("%.1f", percentage) + "%)\n").getBytes());
	        }
	        
	        return new FloatPCAResult(fullEigenVectors, eigenValues, Arrays.copyOf(singularValues, k));
	        
	    } catch (OutOfMemoryError e) {
	        LOG.warn("Randomized SVD failed with " + actualTargetSNPs + " SNPs, halving to " + (actualTargetSNPs / 2));
	        System.gc();
	        
	        // Halve the number of SNPs and try again (minimum 1000 SNPs)
	        int newTargetSNPs = Math.max(1000, actualTargetSNPs / 2);
	        return performRandomizedSVDRecursive(diskMatrix, n, m, k, newTargetSNPs, warningOS);
	    }
	}

	private int[] selectRandomSNPs(int totalSNPs, int targetSNPs) {
	    if (targetSNPs >= totalSNPs) {
	        int[] all = new int[totalSNPs];
	        for (int i = 0; i < totalSNPs; i++) all[i] = i;
	        return all;
	    }
	    Random rand = new Random(42);
	    int[] reservoir = new int[targetSNPs];
	    // Fill initial reservoir
	    for (int i = 0; i < targetSNPs; i++) reservoir[i] = i;
	    // Replace elements with gradually decreasing probability
	    for (int i = targetSNPs; i < totalSNPs; i++) {
	        int j = rand.nextInt(i + 1);
	        if (j < targetSNPs) reservoir[j] = i;
	    }
	    Arrays.sort(reservoir);
	    return reservoir;
	}

    private FMatrixRMaj extractSubsampledMatrix(OptimizedDiskFloatMatrix diskMatrix, int n, int[] selectedIndices) {
        FMatrixRMaj subsampled = new FMatrixRMaj(n, selectedIndices.length);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < selectedIndices.length; j++) {
                subsampled.set(i, j, diskMatrix.get(i, selectedIndices[j]));
            }
        }
        return subsampled;
    }

    private FMatrixRMaj randomizedProjection(OptimizedDiskFloatMatrix diskMatrix, int n, int m, int l) {
        // Generate random matrix Ω (m × l)
        FMatrixRMaj omega = new FMatrixRMaj(m, l);
        Random rand = new Random(42);
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < l; j++) {
                omega.set(i, j, (float) rand.nextGaussian());
            }
        }
        
        // Compute Y = A × Ω without forming full A
        // A is n x m, Ω is m x l, so Y should be n x l
        FMatrixRMaj Y = new FMatrixRMaj(n, l);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < l; j++) {
                float sum = 0;
                for (int k = 0; k < m; k++) {
                    sum += diskMatrix.get(i, k) * omega.get(k, j);
                }
                Y.set(i, j, sum);
            }
        }
        
        return Y; // This returns Q = Y, which is n x l
    }

    private FMatrixRMaj projectToReducedSpace(OptimizedDiskFloatMatrix diskMatrix, FMatrixRMaj Q, 
                                            int n, int m, int l) {
        // Compute B = Q^T × A without forming full A
        // Q is n x l, A is n x m, so B = Q^T × A should be l x m
        FMatrixRMaj B = new FMatrixRMaj(l, m);
        
        // Process in blocks to control memory usage
        int blockSize = 1000;
        for (int colBlock = 0; colBlock < m; colBlock += blockSize) {
            int endCol = Math.min(colBlock + blockSize, m);
            int blockWidth = endCol - colBlock;
            
            // Extract block of A: n x blockWidth
            FMatrixRMaj A_block = new FMatrixRMaj(n, blockWidth);
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < blockWidth; j++) {
                    A_block.set(i, j, diskMatrix.get(i, colBlock + j));
                }
            }
            
            // Multiply Q^T × A_block: (l x n) × (n x blockWidth) = l x blockWidth
            FMatrixRMaj B_block = new FMatrixRMaj(l, blockWidth);
            CommonOps_FDRM.multTransA(Q, A_block, B_block);
            
            // Copy to B
            for (int i = 0; i < l; i++) {
                for (int j = 0; j < blockWidth; j++) {
                    B.set(i, colBlock + j, B_block.get(i, j));
                }
            }
        }
        
        return B; // This returns B which is l x m
    }
    
    private FloatPCAResult performSVDOnReducedMatrix(FMatrixRMaj B, FMatrixRMaj Q, int n, int k) {
        // Create SVD that computes U, S, and V
        SvdImplicitQrDecompose_FDRM svd = new SvdImplicitQrDecompose_FDRM(true, true, true, true);
        if (!svd.decompose(B)) {
            throw new RuntimeException("SVD failed on reduced matrix");
        }
        
        // Get singular values and V matrix
        float[] S_small = svd.getSingularValues();
        FMatrixRMaj V_small = svd.getV(null, false);
        
        // Debug dimensions - this will help us see what's wrong
        LOG.info("Matrix dimensions - Q: " + Q.numRows + "x" + Q.numCols + 
                 ", V_small: " + V_small.numRows + "x" + V_small.numCols +
                 ", B: " + B.numRows + "x" + B.numCols);
        
        // The issue: We're trying to multiply incompatible matrices
        // For PCA, we want the eigenvectors of A^T A, which are the right singular vectors of A
        // So we should return V_small directly or transform it properly
        
        // Since B = Q^T A, and we did SVD on B, the eigenvectors we want are actually
        // the right singular vectors of B, which approximate the right singular vectors of A
        
        // Extract only k components from V_small
        float[] singularValues = Arrays.copyOf(S_small, k);
        float[] eigenValues = new float[k];
        for (int i = 0; i < k; i++) {
            eigenValues[i] = (singularValues[i] * singularValues[i]) / (n - 1);
        }
        
        // V_small contains the approximate eigenvectors for the SNP space
        FMatrixRMaj eigenVectors = new FMatrixRMaj(V_small.numRows, k);
        CommonOps_FDRM.extract(V_small, 0, V_small.numRows, 0, k, eigenVectors, 0, 0);
        
        return new FloatPCAResult(eigenVectors, eigenValues, singularValues);
    }

    private FloatPCAResult performExactSVDWithMemoryCheck(OptimizedDiskFloatMatrix diskMatrix, int k) {
        int R = diskMatrix.rows();
        int C = diskMatrix.cols();

        FMatrixRMaj matrix = null;
        float[] rowBuf = new float[C];

        try {
            matrix = new FMatrixRMaj(R, C);
            float[] matrixData = matrix.data; // Direct access to underlying array

            // Optimized loading: read row-wise into matrixData to improve locality
            for (int i = 0; i < R; i++) {
                diskMatrix.getRow(i, rowBuf);
                int rowOffset = i * C;
                System.arraycopy(rowBuf, 0, matrixData, rowOffset, C);

                // Less frequent memory checks
                if (i % 1000 == 0 && Runtime.getRuntime().freeMemory() < 50 * 1024 * 1024)
                    throw new OutOfMemoryError("Running out of memory during matrix loading");
            }

            // Perform SVD decomposition
            SvdImplicitQrDecompose_FDRM svd = new SvdImplicitQrDecompose_FDRM(true, false, true, true);
            if (!svd.decompose(matrix))
                throw new RuntimeException("SVD failed to decompose the matrix.");

            FMatrixRMaj V = svd.getV(null, false);
            float[] singularValues = svd.getSingularValues();
            int numSamples = svd.numRows();

            // Optimized eigenvalue calculation
            float[] eigenValues = new float[singularValues.length];
            float scale = 1.0f / (numSamples - 1);
            for (int i = 0; i < singularValues.length; i++) {
                float s = singularValues[i];
                eigenValues[i] = (s * s) * scale;
            }

            return new FloatPCAResult(V, eigenValues, singularValues);

        } finally {
            // Help GC: null underlying array
            if (matrix != null) {
                matrix.data = null;
            }
        }
    }


    /**
     * Safe exact SVD implementation
     */
    private FloatPCAResult performExactSVD(OptimizedDiskFloatMatrix diskMatrix) {
        int R = diskMatrix.rows();
        int C = diskMatrix.cols();
        
        LOG.info("Loading matrix into memory for exact SVD: " + R + " × " + C);
        
        // Use direct memory allocation with cleanup guarantee
        FMatrixRMaj matrix = null;
        try {
            matrix = new FMatrixRMaj(R, C);
            
            // Load data with periodic memory checks
            for (int i = 0; i < R; i++) {
                for (int j = 0; j < C; j++) {
                    matrix.set(i, j, diskMatrix.get(i, j));
                }
                
                // Check memory every 100 rows
                if (i % 100 == 0) {
                    long freeMem = Runtime.getRuntime().freeMemory();
                    if (freeMem < 100 * 1024 * 1024) { // Less than 100MB free
                        throw new OutOfMemoryError("Running out of memory during matrix loading");
                    }
                }
            }
            
            LOG.info("Matrix loaded, performing SVD decomposition...");
            
            SvdImplicitQrDecompose_FDRM svd = new SvdImplicitQrDecompose_FDRM(true, false, true, true);
            if (!svd.decompose(matrix)) {
                throw new RuntimeException("SVD failed to decompose the matrix.");
            }

            FMatrixRMaj V = svd.getV(null, false);
            float[] singularValues = svd.getSingularValues();
            int numSamples = svd.numRows();
            
            float[] eigenValues = new float[singularValues.length];
            for (int i = 0; i < singularValues.length; i++) {
                eigenValues[i] = (singularValues[i] * singularValues[i]) / (numSamples - 1);
            }

            return new FloatPCAResult(V, eigenValues, singularValues);
            
        } finally {
            // Explicit cleanup to help GC
            if (matrix != null) {
                matrix.data = null; // Help GC by nulling the data array
            }
        }
    }

    /**
     * Transform data using streaming approach to avoid loading entire matrix
     */
    public float[][] transformDataStreaming(OptimizedDiskFloatMatrix diskMatrix, FMatrixRMaj pcaMatrix, Integer numEigenvectors) {
        int numSamples = diskMatrix.rows();
        int numVariants = diskMatrix.cols();
        
        int numComponents = numEigenvectors != null ? Math.min(numEigenvectors, pcaMatrix.numCols) : pcaMatrix.numCols;
        
        float[][] result = new float[numSamples][numComponents];
        
        // Use block processing for better cache performance
        int blockSize = 1000; // Adjust based on available memory
        int numBlocks = (numSamples + blockSize - 1) / blockSize;
        
        for (int block = 0; block < numBlocks; block++) {
            int startSample = block * blockSize;
            int endSample = Math.min(startSample + blockSize, numSamples);
            
            // Process block of samples
            for (int sample = startSample; sample < endSample; sample++) {
                for (int comp = 0; comp < numComponents; comp++) {
                    float sum = 0;
                    // Process variants in inner loop for better cache locality
                    for (int variant = 0; variant < numVariants; variant++) {
                        sum += diskMatrix.get(sample, variant) * pcaMatrix.get(variant, comp);
                    }
                    result[sample][comp] = sum;
                }
            }
        }
        
        return result;
    }

    public static class FloatPCAResult {
        public final FMatrixRMaj eigenVectors;
        public final float[] eigenValues;
        public final float[] singularValues;

        public FloatPCAResult(FMatrixRMaj ev, float[] eig, float[] sv) {
            this.eigenVectors = ev;
            this.eigenValues = eig;
            this.singularValues = sv;
        }
        
        public FMatrixRMaj getEigenVectors() {
            return eigenVectors;
        }
        
        public float[] getEigenValues() {
            return eigenValues;
        }
        
        public float[] getSingularValues() {
            return singularValues;
        }
    }
    
    /**
     * Optimized disk-backed matrix with memory-mapped storage
     */
    public static class OptimizedDiskFloatMatrix implements Closeable {
        private final File file;
        private final int rows, cols;
        private final MappedByteBuffer mmap;
        private final FloatBuffer buf;

        public OptimizedDiskFloatMatrix(int rows, int cols) throws IOException {
            this.rows = rows;
            this.cols = cols;

            this.file = Files.createTempFile("pca_opt_matrix_", ".bin").toFile();
            this.file.deleteOnExit();

            long size = (long) rows * cols * 4L;
            if (size > Integer.MAX_VALUE) {
                throw new IOException("Matrix too large for memory mapping: " + size + " bytes");
            }

            try (RandomAccessFile raf = new RandomAccessFile(file, "rw");
                 FileChannel fc = raf.getChannel()) {
                mmap = fc.map(FileChannel.MapMode.READ_WRITE, 0, size);
            }
            buf = mmap.asFloatBuffer();
        }

        public void set(int row, int col, float value) {
            buf.put(row * cols + col, value);  // ROW-MAJOR: [sample][variant]
        }

        public float get(int row, int col) {
            return buf.get(row * cols + col);  // ROW-MAJOR: [sample][variant]
        }
        
        /**
         * Fill dest with the row content. Assumes dest.length >= cols.
         * Reads the row in a contiguous loop for better locality.
         */
        public void getRow(int row, float[] dest) {
            int offset = row * cols;
            // Bulk copy via loop (FloatBuffer has no bulk into float[] with offset reliably)
            for (int j = 0; j < cols; j++) {
                dest[j] = buf.get(offset + j);
            }
        }

        /**
         * Write back the row from src (src.length >= cols).
         */
        public void setRow(int row, float[] src) {
            int offset = row * cols;
            for (int j = 0; j < cols; j++) {
                buf.put(offset + j, src[j]);
            }
        }


        public int rows() { return rows; }
        public int cols() { return cols; }

        @Override
        public void close() {
            // Clean up memory mapping
            try {
                java.lang.reflect.Method cleaner = mmap.getClass().getMethod("cleaner");
                cleaner.setAccessible(true);
                Object clean = cleaner.invoke(mmap);
                if (clean != null) {
                    clean.getClass().getMethod("clean").invoke(clean);
                }
            } catch (Exception e) {
                // Fallback to GC
                System.gc();
            }
            
            // Delete temporary file
            if (file.exists()) {
                file.delete();
            }
        }
    }
}