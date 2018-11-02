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
     * @param t
     * @param times
     * @param m the total number of time intervals + 1 (the total number of time change events)
     * @return
     */
    public static int index(Double t, Double[] times, int m) {

        int epoch = Arrays.binarySearch(times, t);

        // If t was not found in array times by binarySearch, then epoch is negative and binarySearch returns (-(insertion point) - 1).
        // The insertion point is the point at which t would be inserted into the array times
        if (epoch < 0) {
            epoch = -epoch - 1;
        }

        // return at most the index of the last interval (m-1)
        return Math.min(epoch, m-1);
    }


    /**
     * Finds the index of the time interval t lies in
     * @param t
     * @param times
     * @return
     */
    public static int indexTimeIntervalBelow(Double t, Double[] times) {

        // Sort the array times
        Arrays.sort(times);

        if(t > times[times.length - 1]){
            throw new RuntimeException("t is bigger than the biggest value in times array. Rework on times array so that it includes the upper bound for t.");
        }

        // Perform binary search on array times
        int epoch = Arrays.binarySearch(times, t);

        // If t was not found in array times by binarySearch, then epoch is negative and binarySearch returns (-(insertion point) - 1).
        // The insertion point is the point at which t would be inserted into the array times, the index of the first element greater than the key
        // Therefore, change epoch in the corresponding interval number.
        if (epoch < -1) {
            epoch = -epoch - 2;
        }

        return epoch;
    }


    /**
     * Finds the index of the time interval t lies in
     * @param t
     * @param times
     * @return
     */
    public static int indexTimeIntervalAbove(Double t, Double[] times) {

        // Sort the array times
        Arrays.sort(times);

        if(t > times[times.length - 1]){
            throw new RuntimeException("t is bigger than the biggest value in times array. Rework on times array so that it includes the upper bound for t.");
        }

        // Perform binary search on array times
        int epoch = Arrays.binarySearch(times, t);

        // If t was not found in array times by binarySearch, then epoch is negative and binarySearch returns (-(insertion point) - 1).
        // The insertion point is the point at which t would be inserted into the array times, the index of the first element greater than the key
        // Therefore, change epoch in the corresponding interval number.
        if (epoch < -1) {
            epoch = -epoch - 1;
        }

        return epoch;
    }


}
