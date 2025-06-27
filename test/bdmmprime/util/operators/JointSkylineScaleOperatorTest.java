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
import bdmmprime.util.priors.SmartZeroExcludingPrior;
import bdmmprime.util.priors.ZeroExcludingPrior;
import beast.base.evolution.operator.ScaleOperator;
import beast.base.inference.*;
import beast.base.inference.distribution.Prior;
import beast.base.inference.distribution.Uniform;
import beast.base.inference.parameter.RealParameter;
import beast.base.util.Randomizer;
import org.junit.Assert;
import org.junit.Test;
import org.xml.sax.SAXException;

import javax.xml.parsers.ParserConfigurationException;
import java.io.IOException;

public class JointSkylineScaleOperatorTest extends OperatorTest {

    @Test
    public void test() throws IOException, ParserConfigurationException, SAXException {
        Randomizer.setSeed(42);

        RealParameter sv1vals = new RealParameter("1 1 2 3 3");
        sv1vals.setLower(0.0);
        sv1vals.setUpper(10.0);
        sv1vals.setID("sv1vals");
        RealParameter sv2vals = new RealParameter("1 2 2 2 2");
        sv2vals.setLower(0.0);
        sv2vals.setUpper(10.0);
        sv2vals.setID("sv2vals");
        RealParameter smvals = new RealParameter("1 2 2 2 2 3");
        smvals.setLower(0.0);
        smvals.setUpper(10.0);
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
        unif.initByName("lower", 0.0,
                "upper", 10.0);

        Prior sv1valsPrior = new SmartZeroExcludingPrior();
        sv1valsPrior.initByName("x", sv1vals, "distr", unif);
        Prior sv2valsPrior = new ZeroExcludingPrior();
        sv2valsPrior.initByName("x", sv2vals, "distr", unif);
        Prior smvalsPrior = new SmartZeroExcludingPrior();
        smvalsPrior.initByName("x", smvals, "distr", unif);

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
            Assert.assertEquals(5.0, testLogger.getMeans(sv1vals)[idx], 0.1);
            Assert.assertEquals(25.0/3.0, testLogger.getVariances(sv1vals)[idx], 0.2);
        }
        for (int idx=0; idx<5; idx++) {
            Assert.assertEquals(5.0, testLogger.getMeans(sv2vals)[idx], 0.1);
            Assert.assertEquals(25.0/3.0, testLogger.getVariances(sv2vals)[idx], 0.2);
        }
        for (int idx=0; idx<6; idx++) {
            Assert.assertEquals(5.0, testLogger.getMeans(smvals)[idx], 0.1);
            Assert.assertEquals(25.0/3.0, testLogger.getVariances(smvals)[idx], 0.2);
        }

        Assert.assertEquals(sv1vals.getArrayValue(0),sv1vals.getArrayValue(1), 1e-10);
        Assert.assertNotEquals(sv1vals.getArrayValue(0),sv1vals.getArrayValue(2), 1e-10);
        Assert.assertEquals(sv1vals.getArrayValue(3),sv1vals.getArrayValue(4), 1e-10);

        Assert.assertNotEquals(sv2vals.getArrayValue(0),sv2vals.getArrayValue(1), 1e-10);
        Assert.assertNotEquals(sv2vals.getArrayValue(1),sv2vals.getArrayValue(2), 1e-10);
        Assert.assertNotEquals(sv2vals.getArrayValue(1),sv2vals.getArrayValue(3), 1e-10);
        Assert.assertNotEquals(sv2vals.getArrayValue(1),sv2vals.getArrayValue(4), 1e-10);

        Assert.assertNotEquals(smvals.getArrayValue(0),smvals.getArrayValue(1), 1e-10);
        Assert.assertEquals(smvals.getArrayValue(1),smvals.getArrayValue(2), 1e-10);
        Assert.assertEquals(smvals.getArrayValue(1),smvals.getArrayValue(3), 1e-10);
        Assert.assertEquals(smvals.getArrayValue(1),smvals.getArrayValue(4), 1e-10);
        Assert.assertNotEquals(smvals.getArrayValue(1),smvals.getArrayValue(5), 1e-10);
    }
}
