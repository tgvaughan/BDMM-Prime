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

package bdmmprime.trajectories.trajevents;

import bdmmprime.trajectories.simulation.NodeFactory;
import beast.base.evolution.tree.Node;
import beast.base.util.Randomizer;

import java.util.ArrayList;
import java.util.List;

public class SamplingEvent extends TrajectoryEvent {

    int type;
    public int nRemoveSamp, nNoRemoveSamp;

    public SamplingEvent(double time, int type, int nRemoveSamp, int nNoRemoveSamp) {
        this.time = time;
        this.type = type;
        this.nRemoveSamp = nRemoveSamp;
        this.nNoRemoveSamp = nNoRemoveSamp;
        this.multiplicity = nRemoveSamp + nNoRemoveSamp;
    }

    @Override
    public void updateState(double[] state) {
        state[type] -= nRemoveSamp;
    }

    @Override
    public void reverseUpdateState(double[] state) {
        state[type] += nRemoveSamp;
    }

    @Override
    public void simulateTreeEvent(double[] state, List<List<Node>> activeLineages, NodeFactory factory,
                                  Boolean untypedTree) {

        // Add nodes corresponding to sampling WITHOUT removal

        List<Node> unsampledLineages = new ArrayList<>(activeLineages.get(type));
        double N = state[type];

        for (int i=0; i<nNoRemoveSamp; i++) {
            Node sampledNode = factory.newLeafNode(type, time);
            double pSampledAncestor = unsampledLineages.size() / N;

            if (pSampledAncestor == 1.0 || (pSampledAncestor > 0.0 && Randomizer.nextDouble() < pSampledAncestor)) {
                Node child = unsampledLineages.remove(Randomizer.nextInt(unsampledLineages.size()));
                activeLineages.get(type).remove(child);

                Node fake = factory.newIntNode(untypedTree ? -1 : type, time);
                fake.addChild(child);
                fake.addChild(sampledNode);

                activeLineages.get(type).add(fake);
            } else {
                activeLineages.get(type).add(sampledNode);
            }

            N -= 1;
        }

        // Add nodes corresponding to sampling WITH removal

        for (int i=0; i<nRemoveSamp; i++)
            activeLineages.get(type).add(factory.newLeafNode(type, time));
    }

    @Override
    public boolean isSamplingEvent() {
        return true;
    }

    @Override
    public String toString() {
        return "SamplingEvent{" +
                "type=" + type +
                ", nRemoveSamp=" + nRemoveSamp +
                ", nNoRemoveSamp=" + nNoRemoveSamp +
                '}';
    }

    @Override
    public String getEventCode(int nTypes) {
        return "S\t" + type + "\tNA\t" + (nRemoveSamp + nNoRemoveSamp);
    }

    @Override
    public String getEventFingerprint(int nTypes) {
        return "S\t" + type + "\tNA\t" + (nRemoveSamp + nNoRemoveSamp);
    }

    @Override
    public SamplingEvent copy() {
        return new SamplingEvent(time, type, nRemoveSamp, nNoRemoveSamp);
    }
}
