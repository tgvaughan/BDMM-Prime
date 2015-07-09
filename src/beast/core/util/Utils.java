package beast.core.util;

import java.util.Arrays;

/**
 * @author dkuh004
 *         Date: May 28, 2012
 *         Time: 11:25:01 AM
 */
public class Utils {

    public static int index(Double t, Double[] times, int m) {

        int epoch = Arrays.binarySearch(times, t);

        if (epoch < 0) {
            epoch = -epoch - 1;
        }

        return Math.min(epoch, m-1); //Math.max((epoch - 1), 0);
    }


}
