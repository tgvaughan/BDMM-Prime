package bdmmprime.util;

import beast.base.util.GammaFunction;
import beast.base.util.Randomizer;

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
     * This is a comparison method that can be used (for example) to construct
     * ordered sets where elements are guaranteed to differ by at least
     * Utils.globalPrecisionThreshold.
     *
     * Note: this comparitor imposes orderings that are inconsistent with equals.
     *
     * @param a first value to compare
     * @param b second value to compare
     * @return 0 a==b to the global precision, or sign(b-a) otherwise
     */
    public static int precisionLimitedComparator(Double a, Double b) {
        if (equalWithPrecision(a,b))
            return 0;
        if (a<b)
            return -1;
        else
            return 1;
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
            acc += Math.exp(logChoose(n, m) + (n - m) * Math.log(1 - p) + m * Math.log(p));
        }

        return m;
    }

    /**
     * Compute (log) Binomial coefficients n choose k in the case that
     * n and k are integers, but generalizes to real valued n and k.
     *
     * I'm using this in place of Binomial.logChoose(int n, int k) which
     * has exactly the same definition but takes integer arguments for some
     * reason.
     *
     * @param n coefficient parameter
     * @param k coefficient parameter
     * @return  log n choose k
     */
    public static double logChoose(double n, double k) {
        return GammaFunction.lnGamma(n + 1.0)
                - GammaFunction.lnGamma(k + 1.0)
                - GammaFunction.lnGamma(n - k + 1.0);
    }
}
