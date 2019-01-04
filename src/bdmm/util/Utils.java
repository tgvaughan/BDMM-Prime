package bdmm.util;

import java.util.Arrays;

/**
 * @author dkuh004
 *         Date: May 28, 2012
 *         Time: 11:25:01 AM
 */
public class Utils {

    /**
     * Finds the index of the time interval t lies in
     * @param t time for which to identify interval
     * @param times interval start times
     * @return
     */
    public static int getIntervalIndex(double t, double[] times) {

        if (t<0.0)
            return 0;

        int index = Arrays.binarySearch(times, t);

        // If t was not found in array times by binarySearch,
        // then epoch is negative and binarySearch returns
        // (-(insertion point) - 1).
        if (index < 0)
            index = -index - 2;

        // return at most the index of the last interval (m-1)
        return Math.min(index, times.length-1);
    }

    public static void reverseDoubleArray(double[] array) {
        double tmp;
        for (int i=0; i<array.length/2; i++) {
            tmp = array[i];
            array[i] = array[array.length - 1 - i];
            array[array.length - 1 - i] = tmp;
        }
    }

    public static <T> void reverseArray(T[] array) {
        T tmp;
        for (int i=0; i<array.length/2; i++) {
            tmp = array[i];
            array[i] = array[array.length - 1 - i];
            array[array.length - 1 - i] = tmp;
        }
    }
}
