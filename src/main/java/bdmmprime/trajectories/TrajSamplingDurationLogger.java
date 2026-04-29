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

package bdmmprime.trajectories;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.inference.CalculationNode;

import java.io.PrintStream;

@Description("Logger to record the duration of trajectory sampling calculations. " +
        "(Beware - this logger is sensitive to where it's placed in the XML!)")
public class TrajSamplingDurationLogger extends CalculationNode implements Loggable {

    public Input<SampledTrajectory> sampledTrajectoryInput = new Input<>(
            "sampledTrajectory",
            "SampledTrajectory object used to sample trajectories " +
                    "consistent with a type-mapped tree.",
            Input.Validate.REQUIRED);

    SampledTrajectory traj;

    @Override
    public void initAndValidate() {
        traj = sampledTrajectoryInput.get();
    }

    @Override
    public void init(PrintStream out) {
        String prefix = traj.getID() != null
                ? traj.getID() + "."
                : "";

        out.print(prefix + "trajSamplingCalcDuration" + "\t");
    }

    @Override
    public void log(long sample, PrintStream out) {
        out.print(traj.getPrevCalculationTime() + "\t");
    }

    @Override
    public void close(PrintStream out) {
    }
}
