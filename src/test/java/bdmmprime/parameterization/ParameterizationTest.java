/*
 * Copyright (C) 2019-2025 ETH Zurich
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package bdmmprime.parameterization;

import beast.base.inference.parameter.RealParameter;
import beast.base.spec.domain.NonNegativeReal;
import beast.base.spec.domain.Real;
import beast.base.spec.domain.UnitInterval;
import beast.base.spec.inference.parameter.RealScalarParam;
import beast.base.spec.inference.parameter.RealVectorParam;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;

public class ParameterizationTest {

    public double TOLERANCE = 1e-20;

    @Test
    public void basicTest() {

		RealScalarParam<NonNegativeReal> originParam = new RealScalarParam<>(2.0, NonNegativeReal.INSTANCE);

		Parameterization parameterization = new CanonicalParameterization();
		parameterization.initByName(
		        "typeSet", new TypeSet(2),
                "processLength", originParam,
                "birthRate", new SkylineVectorParameter(
                        null,
                        new RealVectorParam<>(new double[] {4.0}, NonNegativeReal.INSTANCE), 2),
                "deathRate", new SkylineVectorParameter(
                        null,
                        new RealVectorParam<>(new double[] {3.0}, NonNegativeReal.INSTANCE), 2),
                "birthRateAmongDemes", new SkylineMatrixParameter(
                        null,
                        new RealVectorParam<>(new double[] {0.0}, NonNegativeReal.INSTANCE), 2),
                "migrationRate", new SkylineMatrixParameter(
                        new RealVectorParam<>(new double[] {1.0}, Real.INSTANCE),
                        new RealVectorParam<>(new double[] {0.1, 0.2}, NonNegativeReal.INSTANCE), 2),
                "samplingRate", new SkylineVectorParameter(
                        null,
                        new RealVectorParam<>(new double[] {1.5}, NonNegativeReal.INSTANCE), 2),
                "removalProb", new SkylineVectorParameter(
                        null,
                        new RealVectorParam<>(new double[] {1.0}, UnitInterval.INSTANCE), 2),
                "rhoSampling", new TimedParameter(
                        new RealVectorParam<>(new double[] {originParam.get()}, Real.INSTANCE),
                        new RealVectorParam<>(new double[] {0.0, 0.0}, UnitInterval.INSTANCE)));

		assertEquals(2, parameterization.getTotalIntervalCount());

        assertEquals(2, parameterization.getBirthRates().length);
        assertEquals(2, parameterization.getDeathRates().length);
        assertEquals(2, parameterization.getSamplingRates().length);
        assertEquals(2, parameterization.getRemovalProbs().length);
        assertEquals(2, parameterization.getRhoValues().length);

        for (int interval=0; interval<2; interval++) {
            double migRate = interval < 1 ? 0.1 : 0.2;

            for (int state1 = 0; state1 < 2; state1++) {
                assertEquals(4.0, parameterization.getBirthRates()[interval][state1], TOLERANCE);
                assertEquals(3.0, parameterization.getDeathRates()[interval][state1], TOLERANCE);
                assertEquals(1.5, parameterization.getSamplingRates()[interval][state1], TOLERANCE);
                assertEquals(1.0, parameterization.getRemovalProbs()[interval][state1], TOLERANCE);
                assertEquals(0.0, parameterization.getRhoValues()[interval][state1], TOLERANCE);

                for (int state2 = 0; state2 < 2; state2++) {
                    if (state2 == state1)
                        continue;

                    assertEquals(migRate, parameterization.getMigRates()[interval][state1][state2], TOLERANCE);
                    assertEquals(0.0, parameterization.getCrossBirthRates()[interval][state1][state2], TOLERANCE);
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

        assertEquals(0, parameterization.getIntervalIndex(-3.0));
        assertEquals(0, parameterization.getIntervalIndex(0.0));
        assertEquals(0, parameterization.getIntervalIndex(0.1));
        assertEquals(0, parameterization.getIntervalIndex(1.0));
        assertEquals(1, parameterization.getIntervalIndex(1.1));
        assertEquals(1, parameterization.getIntervalIndex(1.9));
        assertEquals(1, parameterization.getIntervalIndex(2.0));
    }
}
