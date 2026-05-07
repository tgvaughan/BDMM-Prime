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

package bdmmprime.util;

import beast.base.core.BEASTObject;
import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.spec.type.RealScalar;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

public class NodeAgeLogger extends BEASTObject implements Loggable {

    public Input<List<Tree>> treesInput = new Input<>("tree",
            "Tree whose node ages to log.", new ArrayList<>());

    public Input<List<RealScalar<?>>> fsosInput = new Input<>("finalSampleOffset",
            "Final sample offsets for each tree.", new ArrayList<>());

    List<Tree> trees;
    List<RealScalar<?>> fsos;

    int nNodes;

    @Override
    public void initAndValidate() {
        trees = treesInput.get();
        fsos = fsosInput.get();

        if (!fsos.isEmpty() && fsos.size() != trees.size())
            throw new IllegalArgumentException("If finalSampleOffset is provided," +
                    "there must be as many as there are trees.");

        nNodes = 0;
        for (Tree tree : trees)
            nNodes += tree.getNodeCount();
    }

    @Override
    public void init(PrintStream out) {
        int i=0;
        for (Tree tree : trees) {
            for (int j=0; j<tree.getNodeCount(); j++) {
                out.print((i++) + "\t");
            }
        }
    }

    @Override
    public void log(long sample, PrintStream out) {
        for (int treeIdx = 0; treeIdx<trees.size(); treeIdx++) {
            Tree tree = trees.get(treeIdx);
            double fso = fsos.isEmpty() ? 0.0 : fsos.get(treeIdx).get();
            for (Node node : tree.getNodesAsArray())
                out.print((node.getHeight() + fso) + "\t");
        }
    }

    @Override
    public void close(PrintStream out) {

    }
}
