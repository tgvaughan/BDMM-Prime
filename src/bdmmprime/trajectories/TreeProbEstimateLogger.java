/*
 * Copyright (C) 2019-2024 Tim Vaughan
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

package bdmmprime.trajectories;

import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.inference.CalculationNode;

import java.io.PrintStream;

public class TreeProbEstimateLogger extends CalculationNode implements Loggable {

    public Input<SampledTrajectory> sampledTrajectoryInput = new Input<>("sampledTrajectory",
            "Sampled trajectory object.", Input.Validate.REQUIRED);

    @Override
    public void initAndValidate() { }

    @Override
    public void init(PrintStream out) {
        if (getID() != null)
            out.print(getID() + "\t");
        else
            out.print("logProbEst\t");
    }

    @Override
    public void log(long sample, PrintStream out) {
        out.print(sampledTrajectoryInput.get().getLogTreeProbEstimate() + "\t");
    }

    @Override
    public void close(PrintStream out) { }
}
