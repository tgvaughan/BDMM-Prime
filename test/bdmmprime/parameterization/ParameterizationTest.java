package bdmmprime.parameterization;

import beast.core.parameter.RealParameter;
import org.junit.Assert;
import org.junit.Test;

public class ParameterizationTest {

    public double TOLERANCE = 1e-20;

    @Test
    public void basicTest() {

		RealParameter originParam = new RealParameter("2.0");

		Parameterization parameterization = new CanonicalParameterization();
		parameterization.initByName(
		        "typeSet", new TypeSet(2),
                "processLength", originParam,
                "birthRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("4.0"), 2),
                "deathRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("3.0"), 2),
                "birthRateAmongDemes", new SkylineMatrixParameter(
                        null,
                        new RealParameter("0.0"), 2),
                "migrationRate", new SkylineMatrixParameter(
                        new RealParameter("1.0"),
                        new RealParameter("0.1 0.2"), 2),
                "samplingRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("1.5"), 2),
                "removalProb", new SkylineVectorParameter(
                        null,
                        new RealParameter("1.0"), 2),
                "rhoSampling", new TimedParameter(
                        originParam,
                        new RealParameter("0.0 0.0")));

		Assert.assertEquals(2, parameterization.getTotalIntervalCount());

        Assert.assertEquals(2, parameterization.getBirthRates().length);
        Assert.assertEquals(2, parameterization.getDeathRates().length);
        Assert.assertEquals(2, parameterization.getSamplingRates().length);
        Assert.assertEquals(2, parameterization.getRemovalProbs().length);
        Assert.assertEquals(2, parameterization.getRhoValues().length);

        for (int interval=0; interval<2; interval++) {
            double migRate = interval < 1 ? 0.1 : 0.2;

            for (int state1 = 0; state1 < 2; state1++) {
                Assert.assertEquals(4.0, parameterization.getBirthRates()[interval][state1], TOLERANCE);
                Assert.assertEquals(3.0, parameterization.getDeathRates()[interval][state1], TOLERANCE);
                Assert.assertEquals(1.5, parameterization.getSamplingRates()[interval][state1], TOLERANCE);
                Assert.assertEquals(1.0, parameterization.getRemovalProbs()[interval][state1], TOLERANCE);
                Assert.assertEquals(0.0, parameterization.getRhoValues()[interval][state1], TOLERANCE);

                for (int state2 = 0; state2 < 2; state2++) {
                    if (state2 == state1)
                        continue;

                    Assert.assertEquals(migRate, parameterization.getMigRates()[interval][state1][state2], TOLERANCE);
                    Assert.assertEquals(0.0, parameterization.getCrossBirthRates()[interval][state1][state2], TOLERANCE);
                }
            }
        }
    }


    @Test
    public void testGetIntervalIndex() {
        RealParameter originParam = new RealParameter("2.0");

        Parameterization parameterization = new CanonicalParameterization();
		parameterization.initByName(
		        "typeSet", new TypeSet(2),
                "processLength", originParam,
                "birthRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("4.0"), 2),
                "deathRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("3.0"), 2),
                "birthRateAmongDemes", new SkylineMatrixParameter(
                        null,
                        new RealParameter("0.0"), 2),
                "migrationRate", new SkylineMatrixParameter(
                        new RealParameter("1.0"),
                        new RealParameter("0.1 0.2"), 2),
                "samplingRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("1.5"), 2),
                "removalProb", new SkylineVectorParameter(
                        null,
                        new RealParameter("1.0"), 2),
                "rhoSampling", new TimedParameter(
                        originParam,
                        new RealParameter("0.0 0.0")));

        Assert.assertEquals(0, parameterization.getIntervalIndex(-3.0));
        Assert.assertEquals(0, parameterization.getIntervalIndex(0.0));
        Assert.assertEquals(0, parameterization.getIntervalIndex(0.1));
        Assert.assertEquals(0, parameterization.getIntervalIndex(1.0));
        Assert.assertEquals(1, parameterization.getIntervalIndex(1.1));
        Assert.assertEquals(1, parameterization.getIntervalIndex(1.9));
        Assert.assertEquals(1, parameterization.getIntervalIndex(2.0));
    }
}
