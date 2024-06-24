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

package bdmmprime.trajectories.obsevents;

import bdmmprime.parameterization.Parameterization;
import bdmmprime.trajectories.Trajectory;

public class ObservationEndEvent extends ObservedEvent {

    public ObservationEndEvent(double time) {
        this.time = time;
    }

    @Override
    public int[] getNextLineageCounts() {
        throw new RuntimeException("getNextLineageCounts() not implemented for ObservationEndEvent!");
    }

    @Override
    public double applyToTrajectory(Parameterization param, int interval, Trajectory trajectory) {
        throw new RuntimeException("applyToTrajectory() not implemented for ObservationEndEvent!");
    }

    @Override
    public boolean isFinalEvent() {
        return true;
    }
}
