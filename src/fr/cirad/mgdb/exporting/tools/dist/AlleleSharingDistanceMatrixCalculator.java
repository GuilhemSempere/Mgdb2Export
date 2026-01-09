package fr.cirad.mgdb.exporting.tools.dist;

import java.util.Arrays;

import fr.cirad.tools.ProgressIndicator;

/**
 * Allele Sharing Distance (ASD) matrix calculator for Eigenstrat numeric genotypes.
 *
 * Input:
 * - byte[][] genotypes, with values:
 *   0 = hom ref, 1 = het, 2 = hom alt, 9 = missing
 *
 * Optimized:
 * - Lookup table for shared alleles
 * - Upper-triangle storage
 * - Threaded execution
 */
public class AlleleSharingDistanceMatrixCalculator
        extends UpperTriangleDistanceMatrixCalculator {

    private final byte[][] genotypes;
    private final int nMarkers;
    private final int nIndividuals;

    // Shared allele lookup table
    private static final byte[] SHARED_ALLELES_LOOKUP = new byte[100];
    static {
        Arrays.fill(SHARED_ALLELES_LOOKUP, (byte) -1);

        // 0=hom ref, 1=het, 2=hom alt
        SHARED_ALLELES_LOOKUP[0 * 10 + 0] = 2;
        SHARED_ALLELES_LOOKUP[0 * 10 + 1] = 1;
        SHARED_ALLELES_LOOKUP[1 * 10 + 0] = 1;
        SHARED_ALLELES_LOOKUP[0 * 10 + 2] = 0;
        SHARED_ALLELES_LOOKUP[2 * 10 + 0] = 0;
        SHARED_ALLELES_LOOKUP[1 * 10 + 1] = 2;
        SHARED_ALLELES_LOOKUP[1 * 10 + 2] = 1;
        SHARED_ALLELES_LOOKUP[2 * 10 + 1] = 1;
        SHARED_ALLELES_LOOKUP[2 * 10 + 2] = 2;
        // Missing genotypes (9) stay as -1
    }

    public AlleleSharingDistanceMatrixCalculator(byte[][] genotypes, int nMarkers) {
        this.genotypes = genotypes;
        this.nMarkers = nMarkers;
        this.nIndividuals = genotypes.length;
    }

    /**
     * Calculate the ASD distance matrix.
     */
    public double[][] calculate(ProgressIndicator progress) throws InterruptedException {
        return compute(nIndividuals, progress);
    }

    /**
     * Compute a single row of the upper-triangle matrix.
     */
    @Override
    protected int computeRow(int i, double[][] dist) {
        byte[] g1 = genotypes[i];
        int pairs = 0;

        for (int j = i + 1; j < nIndividuals; j++) {
            dist[i][j - i - 1] = computeASD(g1, genotypes[j]);
            pairs++;
        }

        return pairs;
    }

    /**
     * Calculate ASD between two individuals.
     */
    private double computeASD(byte[] g1, byte[] g2) {
        int totalComparisons = 0;
        int sharedAlleles = 0;

        for (int m = 0; m < nMarkers; m++) {
            int a = g1[m];
            int b = g2[m];

            // Skip missing or invalid genotypes
            if (a < 0 || a > 2 || b < 0 || b > 2)
                continue;

            byte s = SHARED_ALLELES_LOOKUP[a * 10 + b];
            totalComparisons++;
            sharedAlleles += s;
        }

        if (totalComparisons == 0)
            return 1.0; // Maximum distance if no comparable markers

        return 1.0 - (sharedAlleles / (2.0 * totalComparisons));
    }

    @Override
    protected int getThreadCount() {
        // ASD calculation is mostly memory-bound; use more threads
        return Math.max(1, Runtime.getRuntime().availableProcessors() / 2);
    }
}
