package fr.cirad.mgdb.exporting.tools.dist;

import java.util.Arrays;
import java.util.List;

import fr.cirad.tools.ProgressIndicator;

/**
 * Jukes–Cantor distance matrix calculator.
 *
 * Input: aligned nucleotide sequences.
 * Optimized for:
 *  - numeric base encoding
 *  - minimal branching
 *  - hot-loop efficiency
 */
public class JukesCantorDistanceMatrixCalculator
        extends UpperTriangleDistanceMatrixCalculator {

    private final byte[][] sequences;
    private final int seqLen;
    private final int n;

    private static final byte[] BASE_TO_CODE = new byte[256];
    static {
        Arrays.fill(BASE_TO_CODE, (byte) -1);
        BASE_TO_CODE['A'] = 0;
        BASE_TO_CODE['C'] = 1;
        BASE_TO_CODE['G'] = 2;
        BASE_TO_CODE['T'] = 3;
    }

    public JukesCantorDistanceMatrixCalculator(List<String> seqs) {
        this.n = seqs.size();
        this.seqLen = seqs.get(0).length();
        this.sequences = new byte[n][seqLen];

        for (int i = 0; i < n; i++) {
            char[] s = seqs.get(i).toCharArray();
            for (int j = 0; j < seqLen; j++)
                sequences[i][j] = (byte) s[j];
        }
    }

    public double[][] calculate(ProgressIndicator progress) throws InterruptedException {
        return compute(n, progress);
    }

    @Override
    protected int computeRow(int i, double[][] dist) {
        byte[] s1 = sequences[i];
        int count = 0;

        for (int j = i + 1; j < n; j++) {
            dist[i][j - i - 1] =
                    jukesCantorDistance(s1, sequences[j]);
            count++;
        }
        return count;
    }

    @Override
    protected int getThreadCount() {
        return Math.max(1, Runtime.getRuntime().availableProcessors() / 4);
    }

    private double jukesCantorDistance(byte[] s1, byte[] s2) {
        int valid = 0;
        int diff = 0;

        for (int k = 0; k < seqLen; k++) {
            int a = BASE_TO_CODE[s1[k]];
            int b = BASE_TO_CODE[s2[k]];
            if (a >= 0 && b >= 0) {
                valid++;
                diff += (a ^ b) != 0 ? 1 : 0;
            }
        }

        if (valid == 0)
            return 3.0;

        double p = (double) diff / valid;
        return (p == 0.0)
                ? 0.0
                : -0.75 * Math.log(1.0 - (4.0 / 3.0) * p);
    }
}
