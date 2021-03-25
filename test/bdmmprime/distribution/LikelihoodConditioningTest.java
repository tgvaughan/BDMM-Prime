package bdmmprime.distribution;

import bdmmprime.parameterization.EpiParameterization;
import bdmmprime.parameterization.Parameterization;
import bdmmprime.parameterization.SkylineVectorParameter;
import beast.core.parameter.RealParameter;
import beast.util.TreeParser;
import org.junit.Test;

import static org.junit.Assert.assertEquals;

public class LikelihoodConditioningTest {

    @Test
    public void testSampleCountProbabilities() {

        Parameterization param = new EpiParameterization();
        param.initByName(
                "R0",
                new SkylineVectorParameter(
                        new RealParameter("2.0"),
                        new RealParameter("2.0 1.8"), 1),
                "becomeUninfectiousRate",
                new SkylineVectorParameter(null,
                        new RealParameter("1.0"), 1),
                "samplingProportion",
                new SkylineVectorParameter(
                        new RealParameter("2.0"),
                        new RealParameter("0.1 0.2"), 1),
                "removalProb",
                new SkylineVectorParameter(null,
                        new RealParameter("1.0"), 1),
                "origin", new RealParameter("5.0"));

        BirthDeathMigrationDistribution distr = new BirthDeathMigrationDistribution();
        distr.initByName(
                "tree", new TreeParser(
                        "((A:1,B:1):1,C:1):0;",
                        false,
                        false,
                        true,
                        1),
                "parameterization", param);

        // Test against values computed independently in R:

        assertEquals(0.414446680,
                Math.exp(distr.calculateLogSampleProb(0, 5.0, 5)),
                1e-5);

        assertEquals( 0.014046736,
                Math.exp(distr.calculateLogSampleProb(3, 5.0, 5)),
                1e-5);

        assertEquals(0.011529376,
                Math.exp(distr.calculateLogSampleProb(7, 5.0, 5)),
                1e-5);

        assertEquals(0.3042198,
                Math.exp(distr.calculateLogMinSampleProb(20, 5.0)),
                1e-5);
    }
}
