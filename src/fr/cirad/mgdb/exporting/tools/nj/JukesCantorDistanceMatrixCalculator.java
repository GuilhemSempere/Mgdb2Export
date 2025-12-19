package fr.cirad.mgdb.exporting.tools.nj;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.*;
import java.util.concurrent.atomic.LongAdder;

import org.apache.log4j.Logger;

import fr.cirad.tools.ProgressIndicator;

/**
 * Safe, server-friendly Jukes-Cantor distance calculator using
 * UPPER-TRIANGLE storage to reduce memory usage by ~50%.
 *
 * Optimized calculator with:
 * - Dedicated bounded thread pool
 * - Lookup table for character validation (no branches in hot loop)
 * - char[] arrays for direct access
 * - Progress cancellation support at row level only
 * - Upper-triangle storage for 50% memory reduction
 * - Reduced contention on progress updates
 */
public class JukesCantorDistanceMatrixCalculator {
	
	private static final Logger LOG = Logger.getLogger(JukesCantorDistanceMatrixCalculator.class);

    private static final int BLOCK_SIZE = 32;          // rows per task
    private static final long PROGRESS_STEP = 50_000;  // distance updates per progress refresh (reduced frequency)

    // Lookup table for fast base validation
    private static final boolean[] IS_VALID_BASE = new boolean[256];
    static {
        IS_VALID_BASE['A'] = true;
        IS_VALID_BASE['T'] = true;
        IS_VALID_BASE['C'] = true;
        IS_VALID_BASE['G'] = true;
    }

    /**
     * Computes the Jukes-Cantor distance matrix using upper-triangle storage.
     *
     * Returned structure:
     *   dist[i][j - i - 1] == distance(i, j) for j > i
     */
    public static double[][] calculateDistanceMatrix(List<String> sequences, ProgressIndicator progress) throws InterruptedException {

        final int n = sequences.size();
        if (n == 0)
            return new double[0][];
        
        final int seqLen = sequences.get(0).length();
        
        // Pre-convert sequences to char arrays for faster access
        final char[][] seqArrays = new char[n][];
        for (int i = 0; i < n; i++)
            seqArrays[i] = sequences.get(i).toCharArray();

        // Upper-triangle allocation
        final double[][] dist = new double[n][];
        for (int i = 0; i < n; i++)
            dist[i] = new double[n - i - 1];

        final long totalPairs = (long) n * (n - 1) / 2;
        final LongAdder completed = new LongAdder();

        final int cpu = Runtime.getRuntime().availableProcessors();
        final int threads = Math.max(1, cpu / 4);
        
        LOG.debug("Launching distance calculation on " + threads + " threads for " +  n + " sequences of length " + seqLen + " (" + totalPairs + " pairs)");

        ExecutorService pool = Executors.newFixedThreadPool(threads);
        List<Future<Double>> futures = new ArrayList<>();

        // Pre-compute block boundaries to avoid Math.min() calls in loop
        List<int[]> blocks = new ArrayList<>();
        for (int start = 0; start < n; start += BLOCK_SIZE) {
            int from = start;
            int to = Math.min(n, start + BLOCK_SIZE);
            blocks.add(new int[]{from, to});
        }

        for (int[] block : blocks) {
            final int from = block[0];
            final int to = block[1];

            futures.add(pool.submit(() -> {
                double localMax = 0.0;
                long localCompleted = 0; // Local counter to reduce contention

                for (int i = from; i < to; i++) {
                    // Check for cancellation only once per row
                    if (progress != null && (progress.getError() != null || progress.isAborted()))
                        return localMax;
                    
                    final char[] seqI = seqArrays[i];

                    for (int j = i + 1; j < n; j++) {
                        double d = calculateJukesCantorDistanceOptimized(seqI, seqArrays[j], seqLen);
                        
                        if (Double.isNaN(d) || d < 0)
                            d = 3.0;

                        dist[i][j - i - 1] = d;

                        if (d > localMax)
                            localMax = d;

                        localCompleted++;
                    }
                    
                    // Update shared counter after completing each row to reduce contention
                    completed.add(localCompleted);
                    long done = completed.sum();
                    localCompleted = 0;
                    
                    // Throttle progress updates to avoid flooding UI
                    if (done % PROGRESS_STEP < (n - i - 1) && progress != null) {
                        int pct = (int) ((done / (double) totalPairs) * 100);
                        progress.setCurrentStepProgress(pct);
                    }
                }
                
                // Add any remaining local count
                if (localCompleted > 0) {
                    completed.add(localCompleted);
                }
                
                return localMax;
            }));
        }

        pool.shutdown();
        
        // Wait with regular cancellation checks
        while (!pool.awaitTermination(1, TimeUnit.SECONDS)) {
            if (progress != null && (progress.getError() != null || progress.isAborted())) {
                pool.shutdownNow();
                return null;
            }
        }

        if (progress != null && (progress.getError() != null || progress.isAborted()))
            return null;

        // Compute global maximum distance
        double maxDist = 0.0;
        for (Future<Double> f : futures)
            try {
                maxDist = Math.max(maxDist, f.get());
            } catch (ExecutionException e) {
                if (progress != null)
                    progress.setError("Distance computation failed: " + e.getCause().getMessage());
                throw new RuntimeException("Distance computation failed", e.getCause());
            }

        if (progress != null)
            progress.setCurrentStepProgress(100);
        
        LOG.debug("Distance matrix calculation completed. Max distance: " + maxDist);
        return dist;
    }

    /**
     * Optimized Jukes-Cantor distance calculation with:
     * - Lookup table for base validation (no branches)
     * - char[] direct access
     * - No progress checks (handled at row level in caller)
     */
    public static double calculateJukesCantorDistanceOptimized(char[] seq1, char[] seq2, final int seqLen) {
        
        int validSites = 0;
        int differences = 0;
        
        for (int i = 0; i < seqLen; i++) {
            char b1 = seq1[i];
            char b2 = seq2[i];
            
            // Fast lookup table validation - no redundant checks
            if (IS_VALID_BASE[b1] && IS_VALID_BASE[b2]) {
                validSites++;
                if (b1 != b2) {
                    differences++;
                }
            }
        }

        if (validSites == 0)
            return -1;

        double p = (double) differences / validSites;
        if (p == 0)
            return 0.0;

        return -0.75 * Math.log(1 - (4.0 / 3.0) * p);
    }

    /**
     * Accessor helper for upper-triangle matrix.
     */
    public static double getDistance(double[][] dist, int i, int j) {
        if (i == j) {
            return 0.0;
        }
        if (i < j) {
            return dist[i][j - i - 1];
        }
        return dist[j][i - j - 1];
    }

    /**
     * Check if base is valid (A, T, C, G)
     * Kept for external use if needed
     */
    public static boolean isValidBase(char b) {
        return b >= 0 && b < 256 && IS_VALID_BASE[b];
    }

    /**
     * Check if base is a gap (- or N)
     * Kept for external use if needed
     */
    public static boolean isGap(char b) {
        return b == '-' || b == 'N';
    }

    /**
     * Check if comparison between two bases is valid
     * Kept for external use if needed
     */
    public static boolean isValidComparison(char b1, char b2) {
        return isValidBase(b1) && isValidBase(b2) && !isGap(b1) && !isGap(b2);
    }
}