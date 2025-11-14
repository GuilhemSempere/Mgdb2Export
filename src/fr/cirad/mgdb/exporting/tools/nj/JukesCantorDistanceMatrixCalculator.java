package fr.cirad.mgdb.exporting.tools.nj;

import java.util.List;
import java.util.concurrent.atomic.LongAdder;
import java.util.stream.IntStream;

import org.apache.log4j.Logger;

import fr.cirad.tools.ProgressIndicator;

public class JukesCantorDistanceMatrixCalculator {

    private static final Logger LOG = Logger.getLogger(JukesCantorDistanceMatrixCalculator.class);

    public static double[][] calculateDistanceMatrix(List<String> sequences, ProgressIndicator progress) {
        final int n = sequences.size();
        final int seqLen = sequences.get(0).length();   // required to be final

        double[][] distanceMatrix = new double[n][n];

        // Progress counters
        final LongAdder completed = new LongAdder();
        final long totalTasks = (long) n * (n - 1) / 2;

        LOG.info("Launching calculateDistanceMatrix on "
                 + Runtime.getRuntime().availableProcessors() + " thread(s)");

        // Parallel row computation
        double maxDist = IntStream.range(0, n).parallel().mapToDouble(i -> {

            double localMax = 0.0;
            String seqI = sequences.get(i);

            // diagonal
            distanceMatrix[i][i] = 0.0;

            for (int j = i + 1; j < n; j++) {
                if (progress.getError() != null || progress.isAborted())
                    return 0; // signal interruption

                double dist = calculateJukesCantorDistance(seqI, sequences.get(j), seqLen);

                if (Double.isNaN(dist))
                    dist = 3.0;

                distanceMatrix[i][j] = dist;
                distanceMatrix[j][i] = dist;

                if (dist > localMax)
                    localMax = dist;

                completed.increment();
                progress.setCurrentStepProgress(
                    (int) ((completed.sum() / (double) totalTasks) * 100)
                );
            }

            return localMax;

        }).max().orElse(0.0);

        if (progress.getError() != null || progress.isAborted())
            return null;

        // Fix undefined values (-1) if any remain
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (distanceMatrix[i][j] == -1)
                    distanceMatrix[i][j] = maxDist * 2;
            }
        }

        return distanceMatrix;
    }

    public static double calculateJukesCantorDistance(String seq1, String seq2, final int seqLen) {
        int validSites = 0;
        int differences = 0;

        for (int i = 0; i < seqLen; i++) {
            char b1 = seq1.charAt(i);
            char b2 = seq2.charAt(i);

            if (isValidComparison(b1, b2)) {
                validSites++;
                if (b1 != b2)
                    differences++;
            }
        }

        if (validSites == 0)
            return -1;

        double p = (double) differences / validSites;
        if (p == 0)
            return 0.0;

        return -0.75 * Math.log(1 - (4.0 / 3.0) * p);
    }

    public static boolean isValidBase(char b) {
        return b == 'A' || b == 'T' || b == 'C' || b == 'G';
    }

    public static boolean isGap(char b) {
        return b == '-' || b == 'N';
    }

    public static boolean isValidComparison(char b1, char b2) {
        return isValidBase(b1) && isValidBase(b2) && !(isGap(b1) || isGap(b2));
    }
}