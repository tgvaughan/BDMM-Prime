package bdmm.util;

import java.util.Arrays;

/**
 * @author dkuh004
 *         Date: May 28, 2012
 *         Time: 11:25:01 AM
 */
public class Utils {

    /**
     * Finds the index of the time interval t lies in.  Note that if t
     * lies on a boundary between intervals, the interval returned will be
     * the _earlier_ of these two intervals.
     *
     * @param t time for which to identify interval
     * @param endTimes interval end times
     * @return index identifying interval.
     */
    public static int getIntervalIndex(double t, double[] endTimes) {

        int index = Arrays.binarySearch(endTimes, t);

        if (index < 0)
            index = -index - 1;

        // return at most the index of the last interval (m-1)
        return Math.max(0, Math.min(index, endTimes.length-1));
    }

    private final static double globalPrecisionThreshold = 1e-10;

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
