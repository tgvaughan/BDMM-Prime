package bdmm.util;

import java.util.Arrays;

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
}
