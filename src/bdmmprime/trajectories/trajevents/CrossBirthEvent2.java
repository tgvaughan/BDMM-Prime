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

package bdmmprime.trajectories.trajevents;

import bdmmprime.trajectories.simulation.NodeFactory;
import beast.base.evolution.tree.Node;
import beast.base.util.Randomizer;

import java.util.List;

public class CrossBirthEvent2 extends TrajectoryEvent {

    int srcType, destType;

    public CrossBirthEvent2(double time, int srcType, int destType, int multiplicity) {
        this.time = time;
        this.srcType = srcType;
        this.destType = destType;
        this.multiplicity = multiplicity;
    }

    public CrossBirthEvent2(double time, int srcType, int destType) {
        this(time, srcType, destType, 1);
    }

    @Override
    public void updateState(double[] state) {
        state[destType] += multiplicity;
    }

    @Override
    public void reverseUpdateState(double[] state) {
        state[destType] -= multiplicity;
    }

    @Override
    public void simulateSingleTreeEvent(double[] state, List<List<Node>> activeLineages, NodeFactory nodeFactory,
                                        Boolean untypedTree) {
        if (activeLineages.get(destType).isEmpty())
            return;

        double pObsStateChange = activeLineages.get(destType).size()/state[destType];
        double pCoal = pObsStateChange*activeLineages.get(srcType).size()/state[srcType];

        double u = Randomizer.nextDouble();

        if (u < pCoal) {
            // Coalescence

            Node child1 = activeLineages.get(srcType).remove(Randomizer.nextInt(activeLineages.get(srcType).size()));
            Node child2 = activeLineages.get(destType).remove(Randomizer.nextInt(activeLineages.get(destType).size()));

            Node parent = nodeFactory.newIntNode(untypedTree ? -1 : srcType, time);
            parent.addChild(child1);
            parent.addChild(child2);

            activeLineages.get(srcType).add(parent);

        } else if (u < pObsStateChange) {
            // Lineage state change

            Node child = activeLineages.get(destType).remove(Randomizer.nextInt(activeLineages.get(destType).size()));
            if (untypedTree) {
                activeLineages.get(srcType).add(child);
            } else {
                Node parent = nodeFactory.newIntNode(srcType, time);
                parent.addChild(child);
                activeLineages.get(srcType).add(parent);
            }
        }
    }

    @Override
    public String toString() {
        return "CrossBirthEvent{" +
                "srcType=" + srcType +
                ", destType=" + destType +
                ", multiplicity=" + multiplicity +
                "}";
    }

    @Override
    public String getEventCode(int nTypes) {
        return "B\t" + srcType + "\t" + destType + "\t" + multiplicity;
    }
}
