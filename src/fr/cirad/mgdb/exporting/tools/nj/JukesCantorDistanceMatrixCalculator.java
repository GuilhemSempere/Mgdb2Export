package fr.cirad.mgdb.exporting.tools.nj;

import java.util.concurrent.RecursiveTask;
import java.util.concurrent.ForkJoinPool;
import java.util.concurrent.atomic.AtomicInteger;

import org.apache.log4j.Logger;

import fr.cirad.tools.ProgressIndicator;

public class JukesCantorDistanceMatrixCalculator {

    private static final Logger LOG = Logger.getLogger(JukesCantorDistanceMatrixCalculator.class);
    
    private static AtomicInteger completedTasks = new AtomicInteger(0);
    private static int totalTasks;

    static public double[][] calculateDistanceMatrix(String[] sequences, ProgressIndicator progress) {
        int n = sequences.length;
        double[][] distanceMatrix = new double[n][n];
        int nProcessCount = Math.max(1, (int) Math.ceil(Runtime.getRuntime().availableProcessors() / 3));
        LOG.info("Launching calculateDistanceMatrix on " + nProcessCount + " thread(s)");
        ForkJoinPool pool = new ForkJoinPool(nProcessCount);

        totalTasks = n * (n - 1) / 2; // Total number of distance calculations needed
        completedTasks.set(0); // Reset the counter

        DistanceCalculatorTask task = new DistanceCalculatorTask(sequences, distanceMatrix, 0, n, progress);
        Double maxDistance = pool.invoke(task);
    	if (progress.getError() != null || progress.isAborted())
    		return null;

        // fixed undefined distances
        for (int i = 0; i < distanceMatrix.length; i++) {
            for (int j = 0; j < distanceMatrix[i].length; j++) {
                if (distanceMatrix[i][j] == -1) {
                	distanceMatrix[i][j] = maxDistance * 2;
                }
            }
        }
        return distanceMatrix;
    }

    static public double calculateJukesCantorDistance(String seq1, String seq2) {
        int validSites = 0;
        int differences = 0;

        for (int i = 0; i < seq1.length(); i++) {
            char base1 = seq1.charAt(i);
            char base2 = seq2.charAt(i);

            if (isValidComparison(base1, base2)) {
                validSites++;
                if (base1 != base2) {
                    differences++;
                }
            }
        }

        if (validSites == 0) {
            return -1; // Indicate undefined distance
        }

        double p = (double) differences / validSites;

        if (p == 0) {
            return 0.0;
        }

        return -0.75 * Math.log(1 - (4.0 / 3.0) * p);
    }

    public static boolean isValidBase(char base) {
        return base == 'A' || base == 'T' || base == 'C' || base == 'G';
    }

    public static boolean isGap(char base) {
        return base == '-' || base == 'N';
    }

    public static boolean isValidComparison(char base1, char base2) {
        return (isValidBase(base1) && isValidBase(base2)) && !(isGap(base1) || isGap(base2));
    }

    static class DistanceCalculatorTask extends RecursiveTask<Double> {
        private static final int THRESHOLD = 10;
        private final String[] sequences;
        private final double[][] distanceMatrix;
        private final int start;
        private final int end;
        private final ProgressIndicator progress;

        DistanceCalculatorTask(String[] sequences, double[][] distanceMatrix, int start, int end, ProgressIndicator progress) {
            this.sequences = sequences;
            this.distanceMatrix = distanceMatrix;
            this.start = start;
            this.end = end;
            this.progress = progress;
        }

        @Override
        protected Double compute() {
            if (end - start <= THRESHOLD) {
                double maxDistance = 0.0;
                for (int i = start; i < end; i++) {
                    for (int j = i; j < sequences.length; j++) {
                        if (i != j) {
                            double distance = calculateJukesCantorDistance(sequences[i], sequences[j]);
                            if (Double.isNaN(distance))
                            	distance = 3;
                            else if (distance > maxDistance)
                                maxDistance = distance;
                            distanceMatrix[i][j] = distance;
                            distanceMatrix[j][i] = distance;
                            progress.setCurrentStepProgress((int) ((completedTasks.getAndIncrement() / (double) totalTasks) * 100));
                        }
                    	if (progress.getError() != null || progress.isAborted())
                    		return null;
                    }
                }
                return maxDistance;
            } else {
                int mid = (start + end) / 2;
                DistanceCalculatorTask task1 = new DistanceCalculatorTask(sequences, distanceMatrix, start, mid, progress);
                DistanceCalculatorTask task2 = new DistanceCalculatorTask(sequences, distanceMatrix, mid, end, progress);
                invokeAll(task1, task2);

            	if (progress.getError() != null || progress.isAborted())
            		return null;
                return Math.max(task1.join(), task2.join());
            }
        }
    }
}