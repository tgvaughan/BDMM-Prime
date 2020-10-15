package bdmmprime.util;

import beast.math.Binomial;
import beast.util.Randomizer;

/**
 * @author dkuh004
 *         Date: May 28, 2012
 *         Time: 11:25:01 AM
 */
public class Utils {


    public final static double globalPrecisionThreshold = 1e-10;

    /**
     * Determine whether a and b are equal to within a precision threshold.
     *
     * @param a first number
     * @param b second number
     * @return true if a and b can be considered equal.
     */
    public static boolean equalWithPrecision(double a, double b) {
        return Math.abs(a-b) < globalPrecisionThreshold;
    }

    /**
     * Returns true when a-b exceeds epsilon.
     *
     * @param a first number
     * @param b second number
     * @return true if a > b + epsilon
     */
    public static boolean greaterThanWithPrecision(double a, double b) {
        return a > b + globalPrecisionThreshold;
    }

    /**
     * Returns true when b-a exceeds epsilon.
     *
     * @param a first number
     * @param b second number
     * @return true if a < b - epsilon
     */
    public static boolean lessThanWithPrecision(double a, double b) {
        return greaterThanWithPrecision(b, a);
    }

    /**
     * In-place reversal of array of (little d) doubles.
     * @param array array to reverse
     */
    public static void reverseDoubleArray(double[] array) {
        double tmp;
        for (int i=0; i<array.length/2; i++) {
            tmp = array[i];
            array[i] = array[array.length - 1 - i];
            array[array.length - 1 - i] = tmp;
        }
    }

    /**
     * In-place reversal of an array.
     * @param array array to reverse.
     * @param <T> type of elements in array
     */
    public static <T> void reverseArray(T[] array) {
        T tmp;
        for (int i=0; i<array.length/2; i++) {
            tmp = array[i];
            array[i] = array[array.length - 1 - i];
            array[array.length - 1 - i] = tmp;
        }
    }

    /**
     * Inefficient (expected time complexity O(np)) ICDF-based binomial sampler.
     *
     * @param n number of trials
     * @param p success probability
     * @return sampled number of successes
     */
    public static int nextBinomial(int n, double p) {
        double u = Randomizer.nextDouble();

        double acc = Math.pow(1-p, n);
        int m = 0;

        while (u > acc && m < n) {
            m += 1;
            acc += Math.exp(Binomial.logChoose(n, m) + (n - m) * Math.log(1 - p) + m * Math.log(p));
        }

        return m;
    }
}
