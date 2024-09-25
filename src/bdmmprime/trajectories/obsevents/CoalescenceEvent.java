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
import bdmmprime.trajectories.trajevents.BirthEvent;
import bdmmprime.trajectories.trajevents.CrossBirthEvent;

public class CoalescenceEvent extends ObservedEvent {

    public int childType1, childType2;

    public CoalescenceEvent(double time, int parentType, int childType1, int childType2, int multiplicity) {
        this.time = time;
        this.type = parentType;
        this.childType1 = childType1;
        this.childType2 = childType2;
        this.multiplicity = multiplicity;
    }

    @Override
    public int[] getNextLineageCounts() {
        int[] nextLineageCounts = lineages.clone();

        nextLineageCounts[type] -= multiplicity;
        nextLineageCounts[childType1] += multiplicity;
        nextLineageCounts[childType2] += multiplicity;

        return nextLineageCounts;
    }

    @Override
    public double applyToTrajectory(Parameterization param, int interval, Trajectory trajectory) {
        double logWeightContrib = 0.0;

        int s = type;
        int sc1 = childType1;
        int sc2 = childType2;

        int sp = sc1 != s ? sc1 : sc2;

        if (sp == s) {
            // Birth

            double birth_prop = trajectory.currentState[s]*param.getBirthRates()[interval][s];
            logWeightContrib += Math.log(birth_prop)
                    - Math.log(0.5*(trajectory.currentState[s]*(trajectory.currentState[s] + 1.0)));

            trajectory.addEvent(new BirthEvent(time, s));

        } else {
            // Cross-birth

            double crossbirth_prop = trajectory.currentState[s]*param.getCrossBirthRates2()[interval][s][sp];
            logWeightContrib += Math.log(crossbirth_prop)
                    - Math.log(trajectory.currentState[s]*(trajectory.currentState[sp] + 1.0));

            trajectory.addEvent(new CrossBirthEvent(time, s, sp));
        }

        return logWeightContrib;
    }
}
