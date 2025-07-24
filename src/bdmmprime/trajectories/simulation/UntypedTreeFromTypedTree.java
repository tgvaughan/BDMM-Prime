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

import beast.base.core.Input;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;

public class UntypedTreeFromTypedTree extends Tree {

    public Input<Tree> typedTreeInput = new Input<>("typedTree",
            "Typed tree used to initialize untyped tree.",
            Input.Validate.REQUIRED);

    public Input<String> typeLabelInput = new Input<>("typeLabel",
            "Type label used in typed tree.",
            "type");

    Tree typedTree;
    int nextIntNr;

    String typeLabel;

    @Override
    public void initAndValidate() {

        typedTree = typedTreeInput.get();
        typeLabel = typeLabelInput.get();

        nextIntNr = typedTree.getLeafNodeCount();

        Node newRoot = getUntypedTree(typedTree.getRoot());

        assignFromWithoutID(new Tree(newRoot));

        super.initAndValidate();
    }

    Node getUntypedTree(Node cladeRoot) {

        while (cladeRoot.getChildCount() == 1)
            cladeRoot = cladeRoot.getChild(0);

        Node newCladeRoot = new Node();
        newCladeRoot.setHeight(cladeRoot.getHeight());

        if (cladeRoot.isLeaf()) {
            newCladeRoot.setNr(cladeRoot.getNr());
            newCladeRoot.setID(cladeRoot.getID());
            newCladeRoot.setMetaData(typeLabel,
                    cladeRoot.getMetaData(typeLabel));
            newCladeRoot.metaDataString = cladeRoot.metaDataString;
        } else {
            for (Node child : cladeRoot.getChildren())
                newCladeRoot.addChild(getUntypedTree(child));

            newCladeRoot.setNr(nextIntNr++);
        }

        return newCladeRoot;
    }
}


