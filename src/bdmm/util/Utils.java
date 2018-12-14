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
    public static int index(double t, double[] times) {

        int epoch = Arrays.binarySearch(times, t);

        // If t was not found in array times by binarySearch, then epoch is negative and binarySearch returns (-(insertion point) - 1).
        // The insertion point is the point at which t would be inserted into the array times
        if (epoch < 0) {
            epoch = -epoch - 1;
        }

        // return at most the index of the last interval (m-1)
        return Math.min(epoch, m-1);
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
}
