/*
 * Copyright (c) 2017-2026 ETH Zürich
 *
 * This file is part of bdmm-prime.
 *
 * bdmm-prime is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or (at
 * your option) any later version.
 *
 * bdmm-prime is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with bdmm-prime. If not, see <https://www.gnu.org/licenses/>.
 */

package bdmmprime.util.operators;

import bdmmprime.testclasses.RealVectorParamFromString;
import beast.base.inference.MCMC;
import beast.base.inference.Operator;
import beast.base.inference.State;
import beast.base.spec.domain.Real;
import beast.base.spec.inference.distribution.IID;
import beast.base.spec.inference.distribution.Uniform;
import beast.base.spec.inference.parameter.RealVectorParam;
import beast.base.util.Randomizer;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;

public class ChangeTimeOperatorTest extends OperatorTestParent {

    @Test
    public void test() throws Exception {
        Randomizer.setSeed(26);

        RealVectorParam<Real> changeTimes= new RealVectorParamFromString<>("0.1 0.2 0.3 0.4 0.5", Real.INSTANCE);

        Uniform unif = new Uniform();
        unif.initByName("lower", "-10.0",
                "upper", "10.0");
        IID<?,?,?> prior = new IID<>(changeTimes, unif);

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
            assertEquals(((idx+1)/6.0)*20 - 10, means[idx], 0.05);
        }
    }
}
