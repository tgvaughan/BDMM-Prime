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

import beast.base.inference.MCMC;
import beast.base.inference.Operator;
import beast.base.inference.State;
import beast.base.inference.distribution.Prior;
import beast.base.inference.distribution.Uniform;
import beast.base.inference.parameter.RealParameter;
import beast.base.util.Randomizer;
import org.junit.Assert;
import org.junit.Test;
import org.xml.sax.SAXException;

import javax.xml.parsers.ParserConfigurationException;
import java.io.IOException;

public class ChangeTimeOperatorTest extends OperatorTest {

    @Test
    public void test() throws IOException, ParserConfigurationException, SAXException {
        Randomizer.setSeed(26);

        RealParameter changeTimes= new RealParameter("0.1 0.2 0.3 0.4 0.5");
        changeTimes.setLower(-10.0);
        changeTimes.setUpper(10.0);

        Uniform unif = new Uniform();
        unif.initByName("lower", -10.0,
                "upper", 10.0);
        Prior prior = new Prior();
        prior.initByName("x", changeTimes,
                "distr", unif);

        Operator op = new ChangeTimeOperator();
        op.initByName("changeTimes", changeTimes,
                "weight", 1.0);

        State state = new State();
        state.initByName("stateNode", changeTimes);

        TestLogger testLogger = new TestLogger();
        testLogger.initByName("log", changeTimes, "logEvery", 1);

        MCMC mcmc = new MCMC();
        mcmc.initByName("chainLength", (long)10000000,
                "preBurnin", 100000,
                "state", state,
                "distribution", prior,
                "operator", op,
                "logger", testLogger);

        mcmc.run();

        // Change times should be distributed as per order statistics of a
        // uniform distribution on [-10,10], i.e. beta-distributed:
        Double[] means = testLogger.getMeans(changeTimes);
        for (int idx=0; idx<5; idx++) {
            Assert.assertEquals(((idx+1)/6.0)*20 - 10, means[idx], 0.05);
        }
    }
}
