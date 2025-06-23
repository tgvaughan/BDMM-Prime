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

import bdmmprime.parameterization.SkylineParameter;
import beast.base.core.BEASTInterface;
import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.inference.Operator;
import beast.base.inference.StateNode;
import beast.base.inference.parameter.RealParameter;
import beast.base.util.Randomizer;

import java.util.*;

@Description("Jointly scales all non-zero values within the provided " +
        "SkylineParameters, leaving their change times untouched.  This " +
        "operator assumes that each skyline vector is independent, although " +
        "it pays attention to the \"linkIdenticalValues\" input when counting " +
        "the degrees of freedom for individual vectors.  This operator is " +
        "primarily used in BEAUti-generated BDMM-Prime XMLs to improve mixing. " +
        "Beware when using it elsewhere that its assumptions about independence " +
        "may not hold.")
public class JointSkylineScaleOperator extends Operator {

    public Input<List<SkylineParameter>> skylineParameterInput = new Input<>(
            "skylineParameter",
            "Skyline parameters on which to operate",
            new ArrayList<>());

    public Input<Double> scaleFactorInput = new Input<>("scaleFactor",
            "Scale factor is chosen uniformly between scaleFactor " +
                    "and 1/scaleFactor.",
            0.75);

    List<SkylineParameter> skylineParameters;
    int nClasses;

    @Override
    public void initAndValidate() {

        skylineParameters = skylineParameterInput.get();
        Set<Double> seenValuesSet = new TreeSet<>();
        nClasses = 0;

        for (SkylineParameter param : skylineParameters) {
            Function values = param.skylineValuesInput.get();

            if (!(values instanceof RealParameter)) {
                throw new IllegalArgumentException("JointSkylineOperator may " +
                        "only be applied to SkylineParameters with " +
                        "true RealParameter values. In particular, values with " +
                        "Non-RealParameter Functions are not compatible.");
            }

            if (param.linkIdenticalValuesInput.get()) {
                seenValuesSet.clear();
                for (int i = 0; i < values.getDimension(); i++) {
                    if (values.getArrayValue(i) != 0.0)
                        seenValuesSet.add(values.getArrayValue(i));
                }
                nClasses += seenValuesSet.size();
            } else {
                for (int i=0; i<values.getDimension(); i++) {
                    if (values.getArrayValue(i) != 0.0)
                        nClasses += 1;
                }
            }
        }
    }

    @Override
    public double proposal() {

        double minf = Math.min(scaleFactorInput.get(), 1.0/scaleFactorInput.get());
        double f = Randomizer.nextDouble()*(1.0/minf - minf) + minf;

        for (SkylineParameter param : skylineParameters) {
            RealParameter values = (RealParameter) param.skylineValuesInput.get();

            for (int i=0; i<values.getDimension(); i++) {
                double newVal = values.getValue(i)*f;
                if (newVal > values.getUpper() || newVal < values.getLower())
                    return Double.NEGATIVE_INFINITY;
                values.setValue(i, newVal);
            }
        }
        return (nClasses - 2)*Math.log(f);
    }

    @Override
    public List<StateNode> listStateNodes() {
        final List<StateNode> stateNodes = new ArrayList<>();
        for (SkylineParameter parameter : skylineParameters)
            stateNodes.add((RealParameter) parameter.skylineValuesInput.get());

        return stateNodes;
    }
}
