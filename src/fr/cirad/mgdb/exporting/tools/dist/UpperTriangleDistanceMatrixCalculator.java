package fr.cirad.mgdb.exporting.tools.dist;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.*;
import java.util.concurrent.atomic.LongAdder;

import fr.cirad.tools.ProgressIndicator;

/**
 * Abstract base class for parallel upper-triangle distance matrix computation.
 *
 * This class handles:
 *  - memory layout
 *  - block partitioning
 *  - threading
 *  - progress reporting
 *  - cancellation
 *
 * Subclasses implement ONLY the row-level distance logic.
 */
public abstract class UpperTriangleDistanceMatrixCalculator {

    protected static final int BLOCK_SIZE = 32;
    protected static final long PROGRESS_STEP = 50_000;

    /**
     * Allocate upper-triangle distance matrix.
     */
    protected final double[][] allocateMatrix(int n) {
        double[][] dist = new double[n][];
        for (int i = 0; i < n; i++)
            dist[i] = new double[n - i - 1];
        return dist;
    }

    /**
     * Compute the distance matrix.
     */
    public final double[][] compute(int n, ProgressIndicator progress) throws InterruptedException {
        final double[][] dist = allocateMatrix(n);
        final long totalPairs = (long) n * (n - 1) / 2;
        final LongAdder completed = new LongAdder();

        final int threads = getThreadCount();
        ExecutorService pool = Executors.newFixedThreadPool(threads);
        List<Future<?>> futures = new ArrayList<>();

        for (int start = 0; start < n; start += BLOCK_SIZE) {
            final int from = start;
            final int to = Math.min(n, start + BLOCK_SIZE);

            futures.add(pool.submit(() -> {
                long localCompleted = 0;

                for (int i = from; i < to; i++) {

                    if (progress != null &&
                        (progress.getError() != null || progress.isAborted()))
                        return null;

                    localCompleted += computeRow(i, dist);

                    completed.add(localCompleted);
                    long done = completed.sum();
                    localCompleted = 0;

                    if (progress != null &&
                        done % PROGRESS_STEP < (n - i - 1)) {
                        progress.setCurrentStepProgress((int) ((done / (double) totalPairs) * 100));
                    }
                }
                return null;
            }));
        }

        shutdownAndAwait(pool, futures, progress);
        return dist;
    }

    /**
     * Compute distances for row i (j > i).
     *
     * Must fill dist[i][j - i - 1].
     *
     * @return number of computed pairs for this row
     */
    protected abstract int computeRow(int i, double[][] dist);

    /**
     * Thread count heuristic (subclasses may override).
     */
    protected int getThreadCount() {
        return Math.max(1, Runtime.getRuntime().availableProcessors() / 2);
    }

    private void shutdownAndAwait(
            ExecutorService pool,
            List<Future<?>> futures,
            ProgressIndicator progress)
            throws InterruptedException {

        pool.shutdown();

        while (!pool.awaitTermination(1, TimeUnit.SECONDS)) {
            if (progress != null &&
                (progress.getError() != null || progress.isAborted())) {
                pool.shutdownNow();
                return;
            }
        }

        for (Future<?> f : futures) {
            try {
                f.get();
            } catch (ExecutionException e) {
                if (progress != null)
                    progress.setError(e.getCause().getMessage());
                throw new RuntimeException(e.getCause());
            }
        }

        if (progress != null)
            progress.setCurrentStepProgress(100);
    }
}
