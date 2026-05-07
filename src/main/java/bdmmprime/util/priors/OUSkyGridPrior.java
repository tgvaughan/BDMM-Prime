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

package bdmmprime.util.priors;

import beast.base.inference.Distribution;
import beast.base.core.Input;
import beast.base.inference.State;
import beast.base.spec.domain.NonNegativeReal;
import beast.base.spec.domain.PositiveReal;
import beast.base.spec.type.RealScalar;
import beast.base.spec.type.RealVector;

import java.util.List;
import java.util.Random;

public class OUSkyGridPrior extends Distribution {

    public Input<RealVector<?>> xInput = new Input<>("x",
            "Parameter to place prior on.", Input.Validate.REQUIRED);

    public Input<RealScalar<?>> MInput = new Input<>("M",
            "M parameter for log-normal distribution.", Input.Validate.REQUIRED);

    public Input<RealScalar<? extends PositiveReal>> SInput = new Input<>("S",
            "S parameter for log-normal distribution.", Input.Validate.REQUIRED);

    public Input<Boolean> meanInRealSpaceInput = new Input<>("meanInRealSpace",
            "Same as meanInRealSpace input for log normal distribution.",
            false);

    public Input<RealScalar<? extends NonNegativeReal>> thetaInput = new Input<>("theta",
            "Relaxation parameter for O-U process.", Input.Validate.REQUIRED);

    public OUSkyGridPrior() { }

    RealVector<?> x;
    int n;

    final double logOneOnSqrt2Pi = -0.5*Math.log(2*Math.PI);

    @Override
    public void initAndValidate() {
        x = xInput.get();
        n = x.size();
    }

    @Override
    public double calculateLogP() {
        logP = 0.0;

        // Parameters for O-U process:
        double M = MInput.get().get();
        double S = SInput.get().get();
        double theta = thetaInput.get().get();

        if (meanInRealSpaceInput.get())
            M = Math.log(M);

        double expNegTheta = Math.exp(-theta);
        double S2 = S*S;

        double var = S2*(1.0 - expNegTheta*expNegTheta);
        double logGausNorm = logOneOnSqrt2Pi - 0.5*Math.log(var);

        // Keep track of previous value:
        double prevEl = Math.log(x.get(0));

        // Log normal distribution for initial value
        logP +=  logOneOnSqrt2Pi - Math.log(S) - 0.5*(prevEl-M)*(prevEl-M)/S2 - prevEl;

        for (int i=1; i<n; i++) {

            double el = Math.log(x.get(i));
            double mean = prevEl*expNegTheta + M*(1.0 - expNegTheta);
            double delta = el - mean;

            logP += logGausNorm - 0.5*delta*delta/var - el;

            prevEl = el;
        }

        return logP;
    }

    @Override
    public List<String> getArguments() {
        return null;
    }

    @Override
    public List<String> getConditions() {
        return null;
    }

    @Override
    public void sample(State state, Random random) {
        throw new UnsupportedOperationException("Sampling from OUSkygridPrior not supported.");
    }
}
