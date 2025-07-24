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

package bdmmprime.trajectories.simulation;

import bdmmprime.trajectories.Trajectory;
import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.inference.CalculationNode;

import java.io.PrintStream;

public class SimulatedTrajectoryLogger extends CalculationNode implements Loggable {

    public Input<SimulatedTree> simulatedTreeInput = new Input<>("simulatedTree",
            "Simulated tree whose trajectory you want to log.",
            Input.Validate.REQUIRED);

    public Input<Double> discretizationTimeStepInput = new Input<>(
            "discretizationTimeStep",
            "If provided, trajectory event sequence will be coarse-grained " +
                    "using a time grid with this spacing. Useful for reducing size of " +
                    "trajectory log files.");

    SimulatedTree simulatedTree;

    @Override
    public void initAndValidate() {
        simulatedTree = simulatedTreeInput.get();
    }

    @Override
    public void init(PrintStream out) {
        Trajectory.init(out);
    }

    @Override
    public void log(long sample, PrintStream out) {
        if (simulatedTree.traj == null)
            Trajectory.logEmpty(out);
        else {
            Trajectory trajTolog;
            if (discretizationTimeStepInput.get() != null)
                trajTolog = simulatedTree.traj.getDiscretized(discretizationTimeStepInput.get());
            else
                trajTolog = simulatedTree.traj;

            Trajectory.log(sample, trajTolog.getStateList(),
                    trajTolog.events,
                    simulatedTree.parameterizationInput.get().getTotalProcessLength(),
                    out);
        }
    }

    @Override
    public void close(PrintStream out) {

    }
}
