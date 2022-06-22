package bdmmprime.util.priors;

import beast.core.Distribution;
import beast.core.Function;
import beast.core.Input;
import beast.core.State;
import beast.evolution.tree.Node;
import beast.evolution.tree.TraitSet;
import beast.evolution.tree.TreeInterface;

import java.util.List;
import java.util.Random;

public class TipDatePrior extends Distribution {

    public Input<TreeInterface> treeInput = new Input<>("tree",
            "Tree whose tips we wish to place a prior on.",
            Input.Validate.REQUIRED);

    public Input<Function> finalSampleOffsetInput = new Input<>(
            "finalSampleOffset",
            "Final sample offset",
            Input.Validate.REQUIRED);

    public Input<TraitSet> initialTipDatesInput = new Input<>("initialTipDates",
            "Initial tip dates", Input.Validate.REQUIRED);
    public Input<TraitSet> earlierBoundInput = new Input<>("earlierBound",
            "Initial tip dates", Input.Validate.REQUIRED);
    public Input<TraitSet> laterBoundInput = new Input<>("laterBound",
            "Initial tip dates", Input.Validate.REQUIRED);

    TreeInterface tree;
    TraitSet earlierBound, laterBound;

    Function fso;

    double earlyOffset, lateOffset;

    @Override
    public void initAndValidate() {
        tree = treeInput.get();

        earlierBound = earlierBoundInput.get();
        earlyOffset = getOffset(earlierBound);

        laterBound = laterBoundInput.get();
        lateOffset = getOffset(laterBound);

        fso = finalSampleOffsetInput.get();
    }

    private double getOffset(TraitSet traitSet) {
        return initialTipDatesInput.get().getDate(0) - traitSet.getDate(0);
    }

    @Override
    public double calculateLogP() {
        logP = 0.0;

        for (int nr=0; nr<tree.getLeafNodeCount(); nr++) {
            Node node = tree.getNode(nr);
            double nodeAge = node.getHeight() + fso.getArrayValue();
            double earlyAge = earlierBound.getValue(node.getID()) + earlyOffset;
            double lateAge = laterBound.getValue(node.getID()) + lateOffset;

            if (nodeAge > earlyAge || nodeAge < lateAge) {
                logP = Double.NEGATIVE_INFINITY;
                break;
            }
        }

        return logP;
    }

    @Override
    public List<String> getArguments() {
        return null;
    }

    @Override
    public List<String> getConditions() {
        return null;
    }

    @Override
    public void sample(State state, Random random) {

    }
}