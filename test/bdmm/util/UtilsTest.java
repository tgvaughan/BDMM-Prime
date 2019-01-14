package bdmm.util;

import org.junit.Assert;
import org.junit.Test;

public class UtilsTest {

    @Test
    public void testGetIntervalIndex() {

        double[] endTimes = {0.5, 0.9};

        Assert.assertEquals(0, Utils.getIntervalIndex(-3.0, endTimes));
        Assert.assertEquals(0, Utils.getIntervalIndex(0.0, endTimes));
        Assert.assertEquals(0, Utils.getIntervalIndex(0.1, endTimes));
        Assert.assertEquals(0, Utils.getIntervalIndex(0.5, endTimes));
        Assert.assertEquals(1, Utils.getIntervalIndex(0.7, endTimes));
        Assert.assertEquals(1, Utils.getIntervalIndex(0.9, endTimes));
        Assert.assertEquals(1, Utils.getIntervalIndex(1.0, endTimes));
    }
}
