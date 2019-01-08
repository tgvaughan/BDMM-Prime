package bdmm.util;

import org.junit.Assert;
import org.junit.Test;

public class UtilsTest {

    @Test
    public void testGetIntervalIndex() {

        double[] startTimes = {0.0, 0.5, 0.9};

        Assert.assertEquals(0, Utils.getIntervalIndex(-3.0, startTimes));
        Assert.assertEquals(0, Utils.getIntervalIndex(0.0, startTimes));
        Assert.assertEquals(0, Utils.getIntervalIndex(0.1, startTimes));
        Assert.assertEquals(0, Utils.getIntervalIndex(0.5, startTimes));
        Assert.assertEquals(1, Utils.getIntervalIndex(0.7, startTimes));
        Assert.assertEquals(1, Utils.getIntervalIndex(0.9, startTimes));
        Assert.assertEquals(2, Utils.getIntervalIndex(1.0, startTimes));
    }
}
