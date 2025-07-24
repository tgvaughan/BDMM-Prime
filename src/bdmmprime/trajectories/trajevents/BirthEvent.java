/*
 * Copyright (C) 2019-2024 ETH Zurich
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

public class BirthEvent extends TrajectoryEvent {

    int type;

    public BirthEvent(double time, int type, int multiplicity) {
        this.time = time;
        this.type = type;
        this.multiplicity = multiplicity;
    }

    public BirthEvent(double time, int type) {
        this(time, type, 1);
    }

    @Override
    public void updateState(double[] state) {
        state[type] += multiplicity;
    }

    @Override
    public void reverseUpdateState(double[] state) {
        state[type] -= multiplicity;
    }

    @Override
    public void simulateSingleTreeEvent(double[] state, List<List<Node>> activeLineages, NodeFactory nodeFactory,
                                        Boolean untypedTree) {
        double probCoal = activeLineages.get(type).size()*(activeLineages.get(type).size()-1)
                /(state[type]*(state[type]-1));

        if (Randomizer.nextDouble() >= probCoal)
            return;

        Node child1 = activeLineages.get(type).remove(Randomizer.nextInt(activeLineages.get(type).size()));
        Node child2 = activeLineages.get(type).remove(Randomizer.nextInt(activeLineages.get(type).size()));

        Node parent = nodeFactory.newIntNode(untypedTree ? -1 : type, time);
        parent.addChild(child1);
        parent.addChild(child2);

        activeLineages.get(type).add(parent);
    }

    @Override
    public String toString() {
        return "BirthEvent{" +
                "type=" + type +
                ", multiplicity=" + multiplicity +
                '}';
    }

    @Override
    public String getEventCode(int nTypes) {
        return "B\t" + type + "\t" + type + "\t" + multiplicity;
    }

    @Override
    public String getEventFingerprint() {
        return "B\t" + type + "\t" + type;
    }

    @Override
    public BirthEvent copy() {
        return new BirthEvent(time, type, multiplicity);
    }
}
