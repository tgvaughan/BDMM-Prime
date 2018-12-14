package bdmm.distributions;

import beast.core.parameter.RealParameter;
import org.junit.Test;

public class ParameterizationTest {

    @Test
    public void basicTest() {

		RealParameter originParam = new RealParameter("2.0");

		Parameterization parameterization = new CanonicalParameterization();
		parameterization.initByName(
		        "nTypes", 2,
                "origin", originParam,
                "birthRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("4.0 4.0")),
                "deathRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("3.0 3.0")),
                "birthRateAmongDemes", new SkylineMatrixParameter(
                        null,
                        new RealParameter("0.0 0.0")),
                "migrationRate", new SkylineMatrixParameter(
                        new RealParameter("1.0"),
                        new RealParameter("0.1 0.1 0.2 0.2")),
                "samplingRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("1.5 1.5")),
                "removalProb", new SkylineVectorParameter(
                        null,
                        new RealParameter("1.0 1.0")),
                "rhoSampling", new TimedParameter(
                        originParam,
                        new RealParameter("0.0 0.0")));

		System.out.println("Number of intervals: " + parameterization.getTotalIntervalCount());
		System.out.println("Done.");

    }
}
