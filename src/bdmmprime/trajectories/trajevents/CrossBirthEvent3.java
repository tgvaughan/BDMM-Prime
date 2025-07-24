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

public class CrossBirthEvent3 extends TrajectoryEvent {

    int parentType, child1Type, child2Type;

    public CrossBirthEvent3(double time, int parentType,
                            int child1Type, int child2Type, int multiplicity) {
        this.time = time;
        this.parentType = parentType;
        this.child1Type = child1Type;
        this.child2Type = child2Type;
        this.multiplicity = multiplicity;
    }

    public CrossBirthEvent3(double time, int parentType,
                            int child1Type, int child2Type) {
        this(time, parentType, child1Type, child2Type, 1);
    }

    @Override
    public void updateState(double[] state) {
        state[parentType] -= multiplicity;
        state[child1Type] += multiplicity;
        state[child2Type] += multiplicity;
    }

    @Override
    public void reverseUpdateState(double[] state) {
        state[child1Type] -= multiplicity;
        state[child2Type] -= multiplicity;
        state[parentType] += multiplicity;
    }

    @Override
    public void simulateSingleTreeEvent(double[] state, List<List<Node>> activeLineages,
                                        NodeFactory nodeFactory, Boolean untypedTree) {

        int n1 = activeLineages.get(child1Type).size();
        double N1 = state[child1Type];

        int n2 = activeLineages.get(child2Type).size();
        double N2 = state[child2Type];

        double p1 = n1/N1;
        double p2 = n2/N2;
        double p1and2 = child1Type == child2Type
                ? p1*(n1-1)/(N2-1)
                : p1*p2;

        double p1or2 = p1 + p2 - p1and2;

        double u = Randomizer.nextDouble();

        if (u < p1and2) {
            // Coalescence

            Node child1 = activeLineages.get(child1Type).remove(Randomizer.nextInt(activeLineages.get(child1Type).size()));
            Node child2 = activeLineages.get(child2Type).remove(Randomizer.nextInt(activeLineages.get(child2Type).size()));

            Node parent = nodeFactory.newIntNode(untypedTree ? -1 : parentType, time);
            parent.addChild(child1);
            parent.addChild(child2);

            activeLineages.get(parentType).add(parent);

        } else if (u < p1) {
            // Lineage state change (child 1)

            if (child1Type == parentType)
                return;

            Node child1 = activeLineages.get(child1Type).remove(Randomizer.nextInt(activeLineages.get(child1Type).size()));
            if (untypedTree) {
                activeLineages.get(parentType).add(child1);
            } else {
                Node parent = nodeFactory.newIntNode(parentType, time);
                parent.addChild(child1);
                activeLineages.get(parentType).add(child1);
            }

        } else if (u < p1or2) {
            // Lineage state change (child 2)

            Node child2 = activeLineages.get(child2Type).remove(Randomizer.nextInt(activeLineages.get(child1Type).size()));
            if (untypedTree) {
                activeLineages.get(parentType).add(child2);
            } else {
                Node parent = nodeFactory.newIntNode(parentType, time);
                parent.addChild(child2);
                activeLineages.get(parentType).add(child2);
            }
        }
    }

    @Override
    public String toString() {
        return "CrossBirthEvent3{" +
                "parentType=" + parentType +
                ", child1Type=" + child1Type +
                ", child2Type=" + child2Type +
                ", multiplicity=" + multiplicity +
                "}";
    }

    @Override
    public String getEventCode(int nTypes) {
        return "B\t" + parentType + "\t" + (child1Type*nTypes + child2Type)
                + "\t" + multiplicity;
    }

    @Override
    public String getEventFingerprint(int nTypes) {
        return "B\t" + parentType + "\t" + (child1Type*nTypes + child2Type);
    }

    @Override
    public CrossBirthEvent3 copy() {
        return new CrossBirthEvent3(time, parentType,
                child1Type, child2Type, multiplicity);
    }
}
