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
import beast.base.spec.type.RealScalar;
import beast.base.spec.type.RealVector;

import java.util.List;
import java.util.Random;

public class SkyGridPrior extends Distribution {

    public Input<RealVector<?>> xInput = new Input<>("x",
            "Parameter to place prior on.", Input.Validate.REQUIRED);

    public Input<RealScalar<?>> MInput = new Input<>("M",
            "M parameter for log normal distribution of first element.",
            Input.Validate.REQUIRED);

    public Input<RealScalar<?>> SInput = new Input<>("S",
            "S parameter for log normal distribution of first element.",
            Input.Validate.REQUIRED);

    public Input<RealScalar<?>> sigmaInput = new Input<>("sigma",
            "Standard deviation of increment priors.", Input.Validate.REQUIRED);

    public SkyGridPrior() { }

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

        double sigma = sigmaInput.get().get();
        double M = MInput.get().get();
        double S = SInput.get().get();

        double prevEl = Math.log(x.get(0));

        // Log normal distribution for initial element:
        logP += logOneOnSqrt2Pi - Math.log(S) - prevEl - 0.5*(prevEl - M)*(prevEl - M)/S/S;

        double logGausNorm = logOneOnSqrt2Pi - Math.log(sigma);
        double sigma2 = sigma*sigma;

        for (int i=1; i<n; i++) {
            double el = Math.log(x.get(i));
            double delta = el - prevEl;

            logP += logGausNorm - el - 0.5*delta*delta/sigma2;

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
        throw new UnsupportedOperationException("Sampling from SkyGridPrior not supported.");
    }
}
