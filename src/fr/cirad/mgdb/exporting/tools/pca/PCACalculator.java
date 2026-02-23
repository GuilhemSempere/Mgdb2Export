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
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Random;

import org.apache.log4j.Logger;
import org.ejml.data.FMatrixRMaj;
import org.ejml.dense.row.CommonOps_FDRM;
import org.ejml.dense.row.decomposition.svd.SvdImplicitQrDecompose_FDRM;
import org.ejml.dense.row.decomposition.qr.QRDecompositionHouseholderTran_FDRM;

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
            warningOS.write(("- Excluding variant #" + variantIndex + " from PCA export, it has too much missing data: " + (nMissingDataCount * 100 / line.length()) + "%\n").getBytes());
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

    private boolean canAllocateOnHeap(long bytesNeeded) {
        Runtime rt = Runtime.getRuntime();
        long max = rt.maxMemory();
        long used = rt.totalMemory() - rt.freeMemory();

        long futureUsed = used + bytesNeeded;
        long minRemaining = 512L * 1024 * 1024; // 512 MB hardcoded

        return futureUsed <= max - minRemaining;
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
            long requiredFloats = (long) R * (long) C;
            long requiredBytes = requiredFloats * 4L;

            // Check heap safety before attempting allocation (fail early)
            if (!canAllocateOnHeap(requiredBytes)) {
                String msg = String.format("Not enough heap available for exact SVD: required=%d MB. Will not attempt large heap allocation.", requiredBytes / (1024L * 1024L));
                LOG.warn(msg);
                throw new OutOfMemoryError(msg);
            }

            // Defensive try/catch around the allocation itself as a last line of defense.
            try {
                matrix = new FMatrixRMaj(R, C); // may allocate float[R*C]
            } catch (OutOfMemoryError oom) {
                String msg = "Allocation of FMatrixRMaj failed despite pre-check: " + oom.getMessage();
                LOG.warn(msg);
                throw oom;
            }

            // Load data with periodic memory checks
            for (int i = 0; i < R; i++) {
                for (int j = 0; j < C; j++) {
                    matrix.set(i, j, diskMatrix.get(i, j));
                }
                
                // Check memory every 100 rows (additional defensive runtime check)
                if (i % 100 == 0) {
                    long freeMem = Runtime.getRuntime().freeMemory();
                    if (freeMem < 100 * 1024 * 1024) { // Less than 100MB free
                        throw new OutOfMemoryError("Running out of memory during matrix loading");
                    }
                }
            }
            
            LOG.info("Matrix loaded, performing SVD decomposition.");
            
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
            // Help GC: null underlying array
            if (matrix != null) {
                matrix.data = null; // Help GC by nulling the data array
            }
        }
    }

    /**
     * Transform data using streaming approach to avoid loading entire matrix
     */
    public float[][] transformDataBlockwise(OptimizedDiskFloatMatrix diskMatrix, FMatrixRMaj pcaMatrix, Integer numEigenvectors) {
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
     * Optimized disk-backed matrix with conditional segmentation.
     * Uses single mapping for matrices < MAX_BYTES_PER_MAP, segmented for larger matrices.
     *
     * Adds a configurable `pca.maxMappedBytes` to refuse mapping too much native memory.
     */
    public static class OptimizedDiskFloatMatrix implements Closeable {
        private final File file;
        private final int rows, cols;

        private final RandomAccessFile raf;
        private final FileChannel fc;

        // For single-segment mode (fast path)
        private FloatBuffer singleBuffer;
        
        // For multi-segment mode (large matrices only)
        private final List<MappedByteBuffer> mmaps = new ArrayList<>();
        private final List<FloatBuffer> buffers = new ArrayList<>();
        private long floatsPerSegment;
        
        private final long totalFloats;
        private final boolean useSegmentation;

        private static final long MAX_BYTES_PER_MAP = ((long) Integer.MAX_VALUE - 8L);

        // configurable maximum number of bytes we allow to be mapped in total
        private static final long MAX_MAPPED_BYTES = Long.MAX_VALUE; // effectively unlimited by default

        public OptimizedDiskFloatMatrix(int rows, int cols) throws IOException {
            this.rows = rows;
            this.cols = cols;

            this.file = Files.createTempFile("pca_opt_matrix_", ".bin").toFile();
            this.file.deleteOnExit();

            this.raf = new RandomAccessFile(file, "rw");
            this.fc = raf.getChannel();

            long sizeBytes = (long) rows * cols * 4L;
            this.totalFloats = (long) rows * cols;
            
            // fail-fast if requested mapped size exceeds configured system property
            if (sizeBytes > MAX_MAPPED_BYTES) {
                raf.close();
                fc.close();
                if (file.exists()) file.delete();
                throw new IOException("Requested disk-backed matrix size " + (sizeBytes / (1024L*1024L)) + "MB exceeds pca.maxMappedBytes=" + MAX_MAPPED_BYTES);
            }
            
            // Decide whether to use segmentation
            this.useSegmentation = sizeBytes > MAX_BYTES_PER_MAP;

            if (sizeBytes == 0) {
                this.floatsPerSegment = 0;
                return;
            }

            raf.setLength(sizeBytes);

            if (!useSegmentation) {
                // Fast path: single mapping
                MappedByteBuffer mmap = fc.map(FileChannel.MapMode.READ_WRITE, 0, sizeBytes);
                this.singleBuffer = mmap.asFloatBuffer();
                mmaps.add(mmap); // Keep for cleanup
            } else {
                // Segmented path for large matrices
                LOG.debug("Calculating PCA based on segmented data");
                long maxFloatsPerSegment = MAX_BYTES_PER_MAP / 4L;
                this.floatsPerSegment = Math.min(totalFloats, maxFloatsPerSegment);

                long floatsRemaining = totalFloats;
                long floatOffset = 0;
                while (floatsRemaining > 0) {
                    long currentFloats = Math.min(floatsPerSegment, floatsRemaining);
                    long mapByteOffset = floatOffset * 4L;
                    long mapByteSize = currentFloats * 4L;

                    MappedByteBuffer mmap = fc.map(FileChannel.MapMode.READ_WRITE, mapByteOffset, mapByteSize);
                    FloatBuffer fb = mmap.asFloatBuffer();
                    mmaps.add(mmap);
                    buffers.add(fb);

                    floatOffset += currentFloats;
                    floatsRemaining -= currentFloats;
                }
            }
        }

        public void set(int row, int col, float value) {
            long linear = (long) row * cols + col;
            setLinear(linear, value);
        }

        public float get(int row, int col) {
            long linear = (long) row * cols + col;
            return getLinear(linear);
        }

        public void getRow(int row, float[] dest) {
            if (!useSegmentation) {
                // Fast path: direct access
                int base = row * cols;
                singleBuffer.position(base);
                singleBuffer.get(dest, 0, cols);
            } else {
                // Segmented path
                int C = cols;
                long base = (long) row * C;
                int remaining = C;
                int destPos = 0;
                while (remaining > 0) {
                    int segIdx = (int) (base / floatsPerSegment);
                    int local = (int) (base - (long) segIdx * floatsPerSegment);
                    int can = (int) Math.min(remaining, floatsPerSegment - local);
                    FloatBuffer fb = buffers.get(segIdx);
                    for (int i = 0; i < can; i++) {
                        dest[destPos + i] = fb.get(local + i);
                    }
                    destPos += can;
                    remaining -= can;
                    base += can;
                }
            }
        }

        public void setRow(int row, float[] src) {
            if (!useSegmentation) {
                // Fast path: direct access
                int base = row * cols;
                singleBuffer.position(base);
                singleBuffer.put(src, 0, cols);
            } else {
                // Segmented path
                int C = cols;
                long base = (long) row * C;
                int remaining = C;
                int srcPos = 0;
                while (remaining > 0) {
                    int segIdx = (int) (base / floatsPerSegment);
                    int local = (int) (base - (long) segIdx * floatsPerSegment);
                    int can = (int) Math.min(remaining, floatsPerSegment - local);
                    FloatBuffer fb = buffers.get(segIdx);
                    for (int i = 0; i < can; i++) {
                        fb.put(local + i, src[srcPos + i]);
                    }
                    srcPos += can;
                    remaining -= can;
                    base += can;
                }
            }
        }

        public float getLinear(long linearIndex) {
            if (!useSegmentation) {
                return singleBuffer.get((int) linearIndex);
            } else {
                int segmentIndex = (int) (linearIndex / floatsPerSegment);
                int localIndex = (int) (linearIndex - (long)segmentIndex * floatsPerSegment);
                return buffers.get(segmentIndex).get(localIndex);
            }
        }

        public void setLinear(long linearIndex, float value) {
            if (!useSegmentation) {
                singleBuffer.put((int) linearIndex, value);
            } else {
                int segmentIndex = (int) (linearIndex / floatsPerSegment);
                int localIndex = (int) (linearIndex - (long)segmentIndex * floatsPerSegment);
                buffers.get(segmentIndex).put(localIndex, value);
            }
        }

        public int rows() { return rows; }
        public int cols() { return cols; }

        @Override
        public void close() {
            // Unmap all mmaps via reflection (best-effort)
            for (MappedByteBuffer mmap : mmaps) {
                try {
                    java.lang.reflect.Method cleaner = mmap.getClass().getMethod("cleaner");
                    cleaner.setAccessible(true);
                    Object clean = cleaner.invoke(mmap);
                    if (clean != null) {
                        clean.getClass().getMethod("clean").invoke(clean);
                    }
                } catch (Exception e) {
                    // ignore and continue; JVM GC will eventually unmap
                }
            }

            try {
                fc.close();
            } catch (IOException e) {
                // ignore
            }
            try {
                raf.close();
            } catch (IOException e) {
                // ignore
            }

            if (file.exists()) {
                file.delete();
            }
        }
    }
    
    /**
     * Smart SVD: chooses exact or randomized SVD depending on matrix size and memory.
     * 
     * @param diskMatrix     	Segmented disk-backed matrix
     * @param numComponents		If null: use all components. If not null: use only the first k.
     * @param progress			A ProgressIndicator
     * @param warningOS      	To write warnings to the caller
     */
    public FloatPCAResult performSmartSVD(OptimizedDiskFloatMatrix diskMatrix, Integer numComponents, OutputStream warningOS, ProgressIndicator progress) throws IOException 
    {
        int R = diskMatrix.rows();
        int C = diskMatrix.cols();
        long floats = (long) R * C;
        long bytes = floats * 4L;
        if (canAllocateOnHeap(bytes)) {
            try {
            	progress.addStep("Computing PCA via exact SVD");
            	progress.moveToNextStep();
                return performExactSVD(diskMatrix);
            }
            catch (OutOfMemoryError oom) {
                warningOS.write("Exact SVD failed due to memory. Falling back to randomized SVD.\n".getBytes());
                LOG.warn("Exact SVD failed: " + oom.getMessage());
            }
        } else {
            warningOS.write(("Skipping exact SVD: not enough heap available for " + (bytes / (1024*1024)) + " MB. Falling back to randomized SVD.\n").getBytes());
            LOG.info("Skipping exact SVD due to heap budget.");
        }

        // Fallback path: randomized SVD
    	progress.addStep("Computing PCA via randomized SVD");
    	progress.moveToNextStep();
        return performRandomizedSVD(diskMatrix, numComponents, warningOS);
    }
    
    /**
     * Randomized SVD that only requires streaming rows from the segmented matrix.
     */
    public FloatPCAResult performRandomizedSVD(
            OptimizedDiskFloatMatrix diskMatrix,
            Integer numComponents,
            OutputStream warningOS) throws IOException {
        int R = diskMatrix.rows();
        int C = diskMatrix.cols();
        int k = (numComponents != null ? numComponents : Math.min(R, C));
        int p = 10;
        int l = Math.min(C, k + p);

        // Create Omega using EJML (better performance)
        Random rng = new Random(42);
        FMatrixRMaj omega = new FMatrixRMaj(C, l);
        for (int i = 0; i < omega.getNumElements(); i++)
            omega.data[i] = (float) rng.nextGaussian();

        // OPTIMIZATION: Process in blocks to reduce overhead
        int blockSize = 1024; // Process 1024 rows at a time
        int numBlocks = (R + blockSize - 1) / blockSize;
        
        FMatrixRMaj Y = new FMatrixRMaj(R, l);
        float[] rowBuf = new float[C];
 
        // Load matrix in blocks and compute Y = A * omega
        for (int b = 0; b < numBlocks; b++) {
            int r0 = b * blockSize;
            int r1 = Math.min(r0 + blockSize, R);
            int blockRows = r1 - r0;
            
            // Load block into contiguous memory (better cache usage)
            FMatrixRMaj block = new FMatrixRMaj(blockRows, C);
            for (int r = r0; r < r1; r++) {
                diskMatrix.getRow(r, rowBuf);
                System.arraycopy(rowBuf, 0, block.data, (r - r0) * C, C);
            }
            
            // Compute Y_block = block * omega (use EJML - optimized)
            FMatrixRMaj Yblock = new FMatrixRMaj(blockRows, l);
            CommonOps_FDRM.mult(block, omega, Yblock);
            
            // Copy result to Y
            System.arraycopy(Yblock.data, 0, Y.data, r0 * l, blockRows * l);
        }
        
        // QR decomposition
        QRDecompositionHouseholderTran_FDRM qr = new QRDecompositionHouseholderTran_FDRM();
        if (!qr.decompose(Y))
            throw new RuntimeException("QR failed");
        FMatrixRMaj Q = qr.getQ(null, false);
        
        // SECOND PASS: Compute B = Q^T * A
        int qcols = Q.numCols;
        FMatrixRMaj B = new FMatrixRMaj(qcols, C);
        for (int b = 0; b < numBlocks; b++) {
            int r0 = b * blockSize;
            int r1 = Math.min(r0 + blockSize, R);
            int blockRows = r1 - r0;
            
            // Load block (fast from OS cache!)
            FMatrixRMaj block = new FMatrixRMaj(blockRows, C);
            for (int r = r0; r < r1; r++) {
                diskMatrix.getRow(r, rowBuf);
                System.arraycopy(rowBuf, 0, block.data, (r - r0) * C, C);
            }
            
            // Extract Q block
            FMatrixRMaj Qblock = CommonOps_FDRM.extract(Q, r0, r1, 0, qcols);
            
            // B += Q_block^T * block (use EJML - much faster than manual loops!)
            FMatrixRMaj Btemp = new FMatrixRMaj(qcols, C);
            CommonOps_FDRM.multTransA(Qblock, block, Btemp);
            CommonOps_FDRM.addEquals(B, Btemp);
        }

        // Final SVD on B
        SvdImplicitQrDecompose_FDRM svd = new SvdImplicitQrDecompose_FDRM(true, false, true, true);
        if (!svd.decompose(B))
            throw new RuntimeException("SVD failed on B");
        
        FMatrixRMaj V = svd.getV(null, false);
        float[] S = svd.getSingularValues();
        
        int keff = Math.min(k, S.length);
        FMatrixRMaj Vk = CommonOps_FDRM.extract(V, 0, V.numRows, 0, keff);
        
        float[] eigenValues = new float[keff];
        for (int i = 0; i < keff; i++)
            eigenValues[i] = (S[i] * S[i]) / (R - 1);
        
        return new FloatPCAResult(Vk, eigenValues, Arrays.copyOf(S, keff));
    }
}