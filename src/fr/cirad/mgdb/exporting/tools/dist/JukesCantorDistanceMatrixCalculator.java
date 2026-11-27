package fr.cirad.mgdb.exporting.tools.dist;

import java.util.Arrays;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import fr.cirad.tools.ProgressIndicator;

/**
 * Jukes–Cantor distance matrix calculator for Eigenstrat numeric genotypes.
 *
 * Input:
 * - byte[][] genotypes, with values:
 *   0 = hom ref, 1 = het, 2 = hom alt, 9 = missing
 *
 * Produces results identical to the previous (FASTA / nucleotide) version
 * when fed the equivalent "2 characters per SNP" alignment
 * (hom ref -> RR, hom alt -> AA, het -> RA or AR with a per-SNP-consistent
 * ordering, missing -> non-ACGT), because in that encoding the per-SNP
 * character-mismatch count is fully determined by the genotype pair and is
 * independent of the ref/alt nucleotide identities.
 *
 * Optimized:
 * - Lookup table for per-SNP character mismatches
 * - Upper-triangle storage
 * - Threaded execution
 */
public class JukesCantorDistanceMatrixCalculator extends UpperTriangleDistanceMatrixCalculator {

    private static final Logger LOG = LoggerFactory.getLogger(JukesCantorDistanceMatrixCalculator.class);

    private final byte[][] genotypes;
    private final int nMarkers;
    private final int nIndividuals;

    /**
     * Per-SNP character-mismatch count for a genotype pair.
     *
     * Derived from the 2-char FASTA encoding:
     *   0 -> RR, 2 -> AA, 1 -> RA (or AR, per-SNP-consistent)
     *
     *   gi \ gj    0    1    2
     *   -----------------------
     *     0        0    1    2
     *     1        1    0    1
     *     2        2    1    0
     *
     * Missing genotypes (9) stay at -1 and are skipped.
     */
    private static final byte[] CHAR_MISMATCH_LOOKUP = new byte[100];
    static {
        Arrays.fill(CHAR_MISMATCH_LOOKUP, (byte) -1);

        CHAR_MISMATCH_LOOKUP[0 * 10 + 0] = 0;
        CHAR_MISMATCH_LOOKUP[0 * 10 + 1] = 1;
        CHAR_MISMATCH_LOOKUP[0 * 10 + 2] = 2;

        CHAR_MISMATCH_LOOKUP[1 * 10 + 0] = 1;
        CHAR_MISMATCH_LOOKUP[1 * 10 + 1] = 0;
        CHAR_MISMATCH_LOOKUP[1 * 10 + 2] = 1;

        CHAR_MISMATCH_LOOKUP[2 * 10 + 0] = 2;
        CHAR_MISMATCH_LOOKUP[2 * 10 + 1] = 1;
        CHAR_MISMATCH_LOOKUP[2 * 10 + 2] = 0;
    }

    /**
     * @param genotypes one row per individual, one byte per marker
     *                  (values 0/1/2/9; used as-is, not copied)
     * @param nMarkers  number of markers (columns) in the genotypes array
     */
    public JukesCantorDistanceMatrixCalculator(byte[][] genotypes, int nMarkers) {
        this.genotypes = genotypes;
        this.nMarkers = nMarkers;
        this.nIndividuals = genotypes.length;
    }

    /**
     * Calculate the Jukes–Cantor distance matrix.
     */
    public double[][] calculate(ProgressIndicator progress) throws InterruptedException {
        long before = System.currentTimeMillis();
        double[][] result = compute(nIndividuals, progress);
        LOG.debug("Distance matrix for " + nMarkers + " markers and " + nIndividuals + " individuals computed in " + (System.currentTimeMillis() - before) / 1000d + "s");
        return result;
    }

    /**
     * Compute a single row of the upper-triangle matrix.
     */
    @Override
    protected int computeRow(int i, double[][] dist) {
        byte[] g1 = genotypes[i];
        int count = 0;

        for (int j = i + 1; j < nIndividuals; j++) {
            dist[i][j - i - 1] = jukesCantorDistance(g1, genotypes[j]);
            count++;
        }
        return count;
    }

    @Override
    protected int getThreadCount() {
        return Math.max(1, Runtime.getRuntime().availableProcessors() / 3);
    }

    /**
     * Jukes–Cantor distance over the equivalent 2-char-per-SNP encoding.
     *
     * 'valid' counts characters (2 per non-missing SNP), matching the
     * original FASTA-based implementation's semantics exactly.
     */
    private double jukesCantorDistance(byte[] g1, byte[] g2) {
        int valid = 0;   // characters: 2 per non-missing SNP pair
        int diff  = 0;

        for (int m = 0; m < nMarkers; m++) {
            int a = g1[m];
            int b = g2[m];

            // Skip missing or invalid genotypes (as before: any non-ACGT char)
            if (a < 0 || a > 2 || b < 0 || b > 2)
                continue;

            diff  += CHAR_MISMATCH_LOOKUP[a * 10 + b];
            valid += 2;
        }

        if (valid == 0)
            return 3.0;

        double p = (double) diff / valid;
        return (p == 0.0)
                ? 0.0
                : -0.75 * Math.log(1.0 - (4.0 / 3.0) * p);
    }
}