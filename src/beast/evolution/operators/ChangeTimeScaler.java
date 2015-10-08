/*
 * Copyright (C) 2015 Tim Vaughan <tgvaughan@gmail.com>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package beast.evolution.operators;

import beast.core.Description;
import beast.core.Input;
import beast.core.Operator;
import beast.core.parameter.RealParameter;
import beast.util.Randomizer;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
@Description("Scale elements of a change time vector, maintaining order " +
        "of elements.")
public class ChangeTimeScaler extends Operator {

    public Input<RealParameter> parameterInput = new Input<>("parameter",
            "Change time parameter to scale.",
            Input.Validate.REQUIRED);

    public Input<Double> scaleFactorInput = new Input<>("scaleFactor",
            "Maximum scale factor used to scale an element.", 0.8);

    RealParameter param;
    double alphaMin, alphaMax;

    @Override
    public void initAndValidate() throws Exception {
        param = parameterInput.get();
        alphaMin = Math.min(scaleFactorInput.get(), 1.0 / scaleFactorInput.get());
        alphaMax = 1.0/alphaMin;
    }

    @Override
    public double proposal() {

        double f = alphaMin + Randomizer.nextDouble()*(alphaMax - alphaMin);

        int idx = Randomizer.nextInt(param.getDimension());

        double lower = idx>0
                ? param.getValue(idx-1)
                : param.getLower();
        double upper = idx<param.getDimension() - 1
                ? param.getValue(idx+1)
                : param.getUpper();

        double x = param.getValue(idx)*f;

        if (param.getValue(idx)>0.0
                && x >= lower
                && x <= upper) {

            param.setValue(idx, x);
            return -Math.log(f);
        } else
            return Double.NEGATIVE_INFINITY;
    }

}
