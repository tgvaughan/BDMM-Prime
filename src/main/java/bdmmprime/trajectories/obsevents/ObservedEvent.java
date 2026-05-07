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

package bdmmprime.trajectories.obsevents;

import bdmmprime.parameterization.Parameterization;
import bdmmprime.trajectories.Trajectory;

public abstract class ObservedEvent {

    public double time;
    public int type, multiplicity;
    public int[] lineages; // before the event

    public boolean isFinalEvent() {
        return false;
    }

    public abstract int[] getNextLineageCounts();

    /**
     * Apply the observed event to the chosen trajectory, returning the probability of the observed
     * event given the trajectory.
     *
     * @param param parameterization object
     * @param interval skyline interval
     * @param trajectory trajectory to update
     * @return log probability of observed event given trajectory.
     */
    public abstract double applyToTrajectory(Parameterization param, int interval, Trajectory trajectory);
}
