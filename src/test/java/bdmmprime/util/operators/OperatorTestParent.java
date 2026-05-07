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

import beast.base.core.BEASTObject;
import beast.base.inference.Logger;
import beast.base.spec.domain.Real;
import beast.base.spec.type.RealVector;

import java.util.*;

public abstract class OperatorTestParent {

    public static class TestLogger extends Logger {


        Map<RealVector<? extends Real>,Double[]> means, variances;
        int N;

        @Override
        public void initAndValidate() {
            means = new HashMap<>();
            variances = new HashMap<>();
        }


        @Override
        public void init() {
            N = 0;
            for (BEASTObject loggable : loggersInput.get()) {
                if (!(loggable instanceof RealVector<?> loggableRealVector))
                    throw new IllegalArgumentException("TestLogger only " +
                            "compatible with objects of type RealVector");

                Double[] theseMeans = new Double[loggableRealVector.size()];
                Arrays.fill(theseMeans, 0.0);
                means.put(loggableRealVector, theseMeans);

                Double[] theseVars = new Double[loggableRealVector.size()];
                Arrays.fill(theseVars, 0.0);
                variances.put(loggableRealVector, theseVars);
            }

        }

        @Override
        public void log(long sampleNr) {
            N += 1;

            for (BEASTObject loggable : loggersInput.get()) {
                RealVector<? extends Real> loggableFunction = (RealVector<? extends Real>) loggable;

                Double[] theseMeans = means.get(loggableFunction);
                Double[] theseVars = variances.get(loggableFunction);

                for (int idx=0; idx<loggableFunction.size(); idx++) {
                    double thisVal = loggableFunction.get(idx);
                    theseMeans[idx] += thisVal;
                    theseVars[idx] += thisVal*thisVal;
                }
            }
        }

        @Override
        public void close() {
            for (BEASTObject loggable : loggersInput.get()) {
                RealVector<? extends Real> loggableFunction = (RealVector<? extends Real>) loggable;

                Double[] theseMeans = means.get(loggableFunction);
                Double[] theseVars = variances.get(loggableFunction);

                for (int idx=0; idx<loggableFunction.size(); idx++) {
                    theseMeans[idx] /= N;
                    theseVars[idx] = theseVars[idx]/N - theseMeans[idx]*theseMeans[idx];
                }
            }
        }

        public Double[] getMeans(RealVector<? extends Real> realVector) {
            return means.get(realVector);
        }

        public Double[] getVariances(RealVector<? extends Real> realVector) {
            return variances.get(realVector);
        }
    }

}
