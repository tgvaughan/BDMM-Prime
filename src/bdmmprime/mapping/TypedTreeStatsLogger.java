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

package bdmmprime.mapping;

import beast.base.evolution.tree.Node;

import java.io.PrintStream;

/**
 * Logger for generating statistics from type mapped trees.
 */
public class TypedTreeStatsLogger extends AbstractTypeTreeStatsLogger {

    int[][] countMatrix;
    double[] lengthVector;

    @Override
    public void initAndValidate() {
        super.initAndValidate();

        countMatrix = new int[nTypes][nTypes];
        lengthVector = new double[nTypes];
    }


    private void update() {
        // Zero length and count arrays
        for (int i=0; i<nTypes; i++) {
            lengthVector[i] = 0.0;
            for (int j=0; j<nTypes; j++) {
                countMatrix[i][j] = 0;
            }
        }

        Node rootNode = tree.getRoot();
        if (!includeRootEdge) {
            while (rootNode.getChildCount() == 1)
                rootNode = rootNode.getChild(0);
        }

        computeMetricsOnSubtree(rootNode);
    }

    private void computeMetricsOnSubtree(Node subtreeRoot) {

        Node node = subtreeRoot;
        int type = getType(subtreeRoot);

        while (node.getChildCount()==1) {
            Node childNode = node.getChild(0);
            int childType = getType(childNode);
            countMatrix[type][childType] += 1;
            lengthVector[childType] += node.getHeight() - childNode.getHeight();

            type = childType;
            node = childNode;
        }

        for (Node childNode : node.getChildren()) {
            int childType = getType(childNode);

            if (childType != type)
                countMatrix[type][childType] += 1;

            lengthVector[childType] += node.getHeight() - childNode.getHeight();

            computeMetricsOnSubtree(childNode);
        }
    }


    @Override
    public void init(PrintStream out) {

        String prefix = tree.getID() != null
                ? tree.getID() + "."
                : "";

        for (int type=0; type<nTypes; type++) {
            out.print(prefix + "length_" + typeSet.getTypeName(type) + "\t");

            for (int typeP=0; typeP<nTypes; typeP++) {
                if (type == typeP)
                    continue;

                out.print(prefix + "count_" + typeSet.getTypeName(type)
                        + "_to_" + typeSet.getTypeName(typeP) + "\t");

            }
        }
    }

    @Override
    public void log(long sample, PrintStream out) {
        if (tree instanceof TypeMappedTree)
            ((TypeMappedTree)tree).remapForLog(sample);

        update();

        for (int type=0; type<nTypes; type++) {
            out.print(lengthVector[type] + "\t");

            for (int typeP = 0; typeP < nTypes; typeP++) {
                if (type == typeP)
                    continue;

                out.print(countMatrix[type][typeP] + "\t");
            }
        }
    }

    @Override
    public void close(PrintStream out) { }
}
