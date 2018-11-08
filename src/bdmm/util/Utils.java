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
    public static int index(double t, double[] times, int m) {

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

    public static void reverseDoubleArray(double[] array) {
        double tmp;
        for (int i=0; i<array.length/2; i++) {
            tmp = array[i];
            array[i] = array[array.length - 1 - i];
            array[array.length - 1 - i] = tmp;
        }
    }



	/**
	 * Obtain offset into "rate matrix" and associated flag arrays.
	 *
	 * @param i row index
	 * @param j column index
     * @param matrixSide side of square matrix
	 * @return Offset
	 */
	public static int getArrayOffset(int i, int j, int matrixSide) {

		if (i==j)
			throw new RuntimeException("Programmer error: requested migration "
					+ "rate array offset for diagonal element of "
					+ "migration rate matrix.");


		if (j>i)
			j -= 1;
		return i*(matrixSide-1)+j;   // todo: check if this is correct!!!
	}

    public static void printArray(double[] array) {
        System.out.print("{");
        for (int i=0; i<array.length; i++) {
            System.out.print(array[i]);
            if (i<array.length-1)
                System.out.print(", ");
        }
        System.out.println("}");
    }

    public static void main(String[] args) {
        double[] testArray = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};

        printArray(testArray);
        reverseDoubleArray(testArray);
        printArray(testArray);
    }

}
