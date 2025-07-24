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

import beast.base.core.BEASTObject;
import beast.base.core.Function;
import beast.base.inference.Logger;

import java.util.*;

public abstract class OperatorTestParent {

    public class TestLogger extends Logger {


        Map<Function,Double[]> means, variances;
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
                if (!(loggable instanceof Function loggableFunction))
                    throw new IllegalArgumentException("TestLogger only " +
                            "compatible with objects of type Function");

                Double[] theseMeans = new Double[loggableFunction.getDimension()];
                Arrays.fill(theseMeans, 0.0);
                means.put(loggableFunction, theseMeans);

                Double[] theseVars = new Double[loggableFunction.getDimension()];
                Arrays.fill(theseVars, 0.0);
                variances.put(loggableFunction, theseVars);
            }

        }

        @Override
        public void log(long sampleNr) {
            N += 1;

            for (BEASTObject loggable : loggersInput.get()) {
                Function loggableFunction = (Function)loggable;

                Double[] theseMeans = means.get(loggableFunction);
                Double[] theseVars = variances.get(loggableFunction);

                for (int idx=0; idx<loggableFunction.getDimension(); idx++) {
                    double thisVal = loggableFunction.getArrayValue(idx);
                    theseMeans[idx] += thisVal;
                    theseVars[idx] += thisVal*thisVal;
                }
            }
        }

        @Override
        public void close() {
            for (BEASTObject loggable : loggersInput.get()) {
                Function loggableFunction = (Function)loggable;

                Double[] theseMeans = means.get(loggableFunction);
                Double[] theseVars = variances.get(loggableFunction);

                for (int idx=0; idx<loggableFunction.getDimension(); idx++) {
                    theseMeans[idx] /= N;
                    theseVars[idx] = theseVars[idx]/N - theseMeans[idx]*theseMeans[idx];
                }
            }
        }

        public Double[] getMeans(Function function) {
            return means.get(function);
        }

        public Double[] getVariances(Function function) {
            return variances.get(function);
        }
    }

}
