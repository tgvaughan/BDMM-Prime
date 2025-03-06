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

import bdmmprime.parameterization.TypeSet;
import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.CalculationNode;

public abstract class AbstractTypeTreeStatsLogger extends CalculationNode implements Loggable {

    public Input<Tree> typedTreeInput = new Input<>("typedTree",
            "Tree with type changes mapped.",
            Input.Validate.REQUIRED);

    public Input<TypeSet> typeSetInput = new Input<>("typeSet",
            "Type set specifying model types.",
            Input.Validate.REQUIRED);

    public Input<String> typeLabelInput = new Input<>("typeLabel",
            "Type label used to store type information on tree.",
            Input.Validate.REQUIRED);

    public Input<Boolean> includeRootEdgeInput = new Input<>(
            "includeRootEdge",
            "If true, include root edge in summary stats calculations.",
            false);

    Tree tree;
    TypeSet typeSet;
    int nTypes;
    String typeLabel;
    boolean includeRootEdge;

    @Override
    public void initAndValidate() {
        tree = typedTreeInput.get();
        typeSet = typeSetInput.get();
        nTypes = typeSet.getNTypes();
        typeLabel = typeLabelInput.get();
        includeRootEdge = includeRootEdgeInput.get();
    }

    protected int getType(Node node) {
        Object typeObj = node.getMetaData(typeLabel);
        if (typeObj == null)
            throw new RuntimeException("Tree does not have type metadata with the label '" + typeLabel + "'");

        if (typeObj instanceof Integer)
            return (int)typeObj;

        if (typeObj instanceof String)
            return typeSet.getTypeIndex((String)typeObj);

        throw new RuntimeException("Tree does not contain valid type metadata.");
    }
}
