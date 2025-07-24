/*
 * Copyright (C) 2022 ETH Zurich
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

package bdmmprime.util.operators;

import beast.base.core.BEASTObject;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.inference.StateNode;
import beast.base.inference.StateNodeInitialiser;
import beast.base.inference.parameter.RealParameter;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.TraitSet;
import beast.base.evolution.tree.Tree;

import java.util.List;

/**
 * Initialise tip dates for a tree with a non-zero final sample offset.
 */
public class TipDateInitialiser extends BEASTObject implements StateNodeInitialiser {


    public Input<Tree> treeInput = new Input<>("tree",
            "Tree whose tips we wish to place a prior on.",
            Input.Validate.REQUIRED);

    public Input<RealParameter> finalSampleOffsetInput = new Input<>(
            "finalSampleOffset",
            "Final sample offset",
            Input.Validate.REQUIRED);

    public Input<TraitSet> tipDatesTraitInput = new Input<>("tipDatesTrait",
            "Initial tip dates", Input.Validate.REQUIRED);

    public Input<Function> endOfSamplingTimeInput = new Input<>("endOfSamplingTime",
            "Time of the point when sampling ends.  (Necessary only " +
                    "when tip times are given forward in time.)");

    public Input<Boolean> adjustInternalNodesInput = new Input<>("adjustInternalNodes",
            "If true, original divergence times of input tree are preserved.",
            false);

    public Input<Function> internalNodeSpacingInput = new Input<>("internalNodeSpacing",
            "Default spacing used when adjusting internal node times.",
            new RealParameter("1.0"));

    Tree tree;
    RealParameter fso;
    TraitSet dateTrait;
    Function endOfSamplingTime;

    boolean traitValuesAreAges, adjustInternalNodes;

    @Override
    public void initAndValidate() {
        tree = treeInput.get();
        dateTrait = tipDatesTraitInput.get();
        fso = finalSampleOffsetInput.get();
        endOfSamplingTime = endOfSamplingTimeInput.get();
        adjustInternalNodes = adjustInternalNodesInput.get();

        traitValuesAreAges = !dateTrait.isDateTrait()
                || dateTrait.getDateType().equals(TraitSet.AGE_TRAIT)
                || dateTrait.getDateType().equals(TraitSet.DATE_BACKWARD_TRAIT);

        if (!traitValuesAreAges && endOfSamplingTime == null)
            throw new IllegalArgumentException("If tip dates are given forward " +
                    "in time, you must also provide a value for the " +
                    "endOfSamplingTime input.");

        initStateNodes();
    }

    @Override
    public void initStateNodes() {

        for (int nr=0; nr<tree.getLeafNodeCount(); nr++) {
            Node node = tree.getNode(nr);
            double nodeAge;
            if (traitValuesAreAges)
                nodeAge = dateTrait.getValue(tree.getTaxonId(node)) + dateTrait.getDate(0);
            else
                nodeAge = dateTrait.getValue(tree.getTaxonId(node)) +
                        (endOfSamplingTime.getArrayValue() - dateTrait.getDate(0));

            node.setHeight(nodeAge - fso.getValue());
            if (!adjustInternalNodes && node.getParent().getHeight() < node.getHeight())
                throw new IllegalStateException("TipDateInitialiser set child (" +
                        tree.getTaxonId(node) + ") older than parent.");
            else {
                Node nodePrime = node;
                while (nodePrime.getParent() != null && nodePrime.getParent().getHeight()<nodePrime.getHeight()) {
                    nodePrime.getParent().setHeight(nodePrime.getHeight()
                            + internalNodeSpacingInput.get().getArrayValue());
                    nodePrime = nodePrime.getParent();
                }
            }
        }

        updateHeights();
    }

    /**
     * Update node heights and final sample offset such that the smallest
     * node height is 0.
     */
    private void updateHeights() {

        // Find height of youngest leaf:
        double lowestHeight = tree.getExternalNodes().stream()
                .mapToDouble(Node::getHeight).summaryStatistics().getMin();

        if (lowestHeight == 0.0)
            return; // FSO has not changed

        double newFSO = fso.getValue() + lowestHeight;

        if (newFSO < Math.max(0.0, fso.getLower())
                || newFSO > fso.getUpper())
            throw new IllegalStateException("TipDateInitialiser tried to " +
                    "set finalSampleOffset to a value outside of its bounds.");

        for (Node node : tree.getNodesAsArray())
            node.setHeight(node.getHeight() - lowestHeight);

        fso.setValue(newFSO);
    }

    @Override
    public void getInitialisedStateNodes(List<StateNode> stateNodes) {
        stateNodes.add(treeInput.get());
        stateNodes.add(finalSampleOffsetInput.get());
    }


}
