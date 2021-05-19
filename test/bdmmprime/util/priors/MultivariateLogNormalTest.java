package bdmmprime.util.priors;

import beast.core.Distribution;
import beast.core.parameter.RealParameter;
import org.junit.Test;

import static org.junit.Assert.assertEquals;

public class MultivariateLogNormalTest {

    @Test
    public void testValue() {

        Distribution dist = new MultivariateLogNormal();
        dist.initByName("x", new RealParameter("1"),
                "x", new RealParameter("2"),
                "M", new RealParameter("0 0"),
                "S", new RealParameter("1 0.1 0.1 1"));

        assertEquals(-2.8181831410741953, dist.calculateLogP(), 1e-10);
    }

    @Test(expected = IllegalArgumentException.class)
    public void testNonSymmetric() {

        Distribution dist = new MultivariateLogNormal();
        dist.initByName("x", new RealParameter("1"),
                "x", new RealParameter("2"),
                "M", new RealParameter("0 0"),
                "S", new RealParameter("1 1.1 1.2 1"));

        System.out.println(dist.calculateLogP());
    }

    @Test
    public void testNonPositiveDefinite() {

        Distribution dist = new MultivariateLogNormal();
        dist.initByName("x", new RealParameter("1"),
                "x", new RealParameter("2"),
                "M", new RealParameter("0 0"),
                "S", new RealParameter("1 1.1 1.1 1"));

        assertEquals(Double.NEGATIVE_INFINITY, dist.calculateLogP(), 0.0);
    }
}
