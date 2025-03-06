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

package bdmmprime.trajectories.simulation;

import bdmmprime.parameterization.TypeSet;
import beast.base.evolution.tree.Node;

public class NodeFactory {

    double origin;
    int nextLeafNr, nextIntNr;

    String typeLabel;
    TypeSet typeSet;

    public NodeFactory(double origin, int nSamples, String typeLabel, TypeSet typeSet) {
        this.origin = origin;
        this.nextIntNr = nSamples;
        this.nextLeafNr = 0;
        this.typeLabel = typeLabel;
        this.typeSet = typeSet;
    }

    private Node newNode(int type, double time, int nextNodeNr) {
        Node node = new Node(String.valueOf(nextLeafNr));
        node.setNr(nextNodeNr);

        node.setHeight(origin-time);

        if (type >= 0) {
            node.setMetaData(typeLabel, type);
            node.metaDataString = String.format("%s=\"%s\"", typeLabel, typeSet.getTypeName(type));
        }

        return node;
    }

    public Node newLeafNode(int type, double time) {
        Node node = newNode(type, time, nextLeafNr);
        node.setID(String.valueOf(nextLeafNr));
        nextLeafNr += 1;
        return node;
    }

    public Node newIntNode(int type, double time) {
        Node node = newNode(type, time, nextIntNr);
        nextIntNr += 1;
        return node;
    }
}
