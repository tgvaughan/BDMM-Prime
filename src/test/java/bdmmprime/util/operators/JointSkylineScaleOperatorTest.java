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

package bdmmprime.util.operators;

import bdmmprime.parameterization.SkylineMatrixParameter;
import bdmmprime.parameterization.SkylineVectorParameter;
import bdmmprime.testclasses.RealVectorParamFromString;
import bdmmprime.util.priors.SmartZeroExcludingRealIID;
import beast.base.inference.*;
import beast.base.spec.domain.NonNegativeReal;
import beast.base.spec.inference.distribution.Uniform;
import beast.base.spec.inference.operator.ScaleOperator;
import beast.base.spec.inference.parameter.RealVectorParam;
import beast.base.util.Randomizer;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotEquals;

public class JointSkylineScaleOperatorTest extends OperatorTestParent {

    @Test
    public void test() throws Exception {
        Randomizer.setSeed(42);

        RealVectorParam<NonNegativeReal> sv1vals = new RealVectorParamFromString<>("1 1 2 3 3", NonNegativeReal.INSTANCE);
        sv1vals.setID("sv1vals");
        RealVectorParam<NonNegativeReal> sv2vals = new RealVectorParamFromString<>("1 2 2 2 2", NonNegativeReal.INSTANCE);
        sv2vals.setID("sv2vals");
        RealVectorParam<NonNegativeReal> smvals = new RealVectorParamFromString<>("1 2 2 2 2 3", NonNegativeReal.INSTANCE);
        smvals.setID("smvals");

        SkylineVectorParameter sv1 = new SkylineVectorParameter(
                null, sv1vals, 5);
        sv1.linkIdenticalValuesInput.setValue(true, sv1);

        SkylineVectorParameter sv2 = new SkylineVectorParameter(
                null, sv2vals, 5);
        sv2.linkIdenticalValuesInput.setValue(false, sv2);

        SkylineMatrixParameter sm = new SkylineMatrixParameter(
                null, smvals, 3);
        sm.linkIdenticalValuesInput.setValue(true, sm);

        Uniform unif = new Uniform();
        unif.initByName("lower", "0.0", "upper", "10.0");

        SmartZeroExcludingRealIID sv1valsPrior = new SmartZeroExcludingRealIID(sv1vals, unif);
        SmartZeroExcludingRealIID sv2valsPrior = new SmartZeroExcludingRealIID(sv2vals, unif);
        SmartZeroExcludingRealIID smvalsPrior = new SmartZeroExcludingRealIID(smvals, unif);

        Distribution target = new CompoundDistribution();
        target.initByName("distribution", sv1valsPrior,
                "distribution", sv2valsPrior,
                "distribution", smvalsPrior);

        Operator sv1Op = new SmartScaleOperator();
        sv1Op.initByName("weight", 1.0,
                "parameter", sv1vals);

        Operator sv2Op = new ScaleOperator();
        sv2Op.initByName("weight", 2.0,
                "parameter", sv2vals);

        Operator smOp = new SmartScaleOperator();
        smOp.initByName("weight", 1.0,
                "parameter", smvals);

        Operator opJoint = new JointSkylineScaleOperator();
        opJoint.initByName("weight", 1.0,
                "skylineParameter", sv1,
                "skylineParameter", sv2,
                "skylineParameter", sm);


        State state = new State();
        state.initByName("stateNode", sv1vals,
                "stateNode", sv2vals,
                "stateNode", smvals);

        TestLogger testLogger = new TestLogger();
        testLogger.initByName("log", sv1vals,
                "log", sv2vals,
                "log", smvals,
                "logEvery", 1);

        MCMC mcmc = new MCMC();
        mcmc.initByName("chainLength", (long)10000000,
                "preBurnin", 10000,
                "state", state,
                "distribution", target,
                "operator", sv1Op,
                "operator", sv2Op,
                "operator", smOp,
                "operator", opJoint,
                "logger", testLogger);

        mcmc.run();

        for (int idx=0; idx<5; idx++) {
            assertEquals(5.0, testLogger.getMeans(sv1vals)[idx], 0.1);
            assertEquals(25.0/3.0, testLogger.getVariances(sv1vals)[idx], 0.2);
        }
        for (int idx=0; idx<5; idx++) {
            assertEquals(5.0, testLogger.getMeans(sv2vals)[idx], 0.1);
            assertEquals(25.0/3.0, testLogger.getVariances(sv2vals)[idx], 0.2);
        }
        for (int idx=0; idx<6; idx++) {
            assertEquals(5.0, testLogger.getMeans(smvals)[idx], 0.1);
            assertEquals(25.0/3.0, testLogger.getVariances(smvals)[idx], 0.2);
        }

        assertEquals(sv1vals.get(0),sv1vals.get(1), 1e-10);
        assertNotEquals(sv1vals.get(0),sv1vals.get(2), 1e-10);
        assertEquals(sv1vals.get(3),sv1vals.get(4), 1e-10);

        assertNotEquals(sv2vals.get(0),sv2vals.get(1), 1e-10);
        assertNotEquals(sv2vals.get(1),sv2vals.get(2), 1e-10);
        assertNotEquals(sv2vals.get(1),sv2vals.get(3), 1e-10);
        assertNotEquals(sv2vals.get(1),sv2vals.get(4), 1e-10);

        assertNotEquals(smvals.get(0),smvals.get(1), 1e-10);
        assertEquals(smvals.get(1),smvals.get(2), 1e-10);
        assertEquals(smvals.get(1),smvals.get(3), 1e-10);
        assertEquals(smvals.get(1),smvals.get(4), 1e-10);
        assertNotEquals(smvals.get(1),smvals.get(5), 1e-10);
    }
}
