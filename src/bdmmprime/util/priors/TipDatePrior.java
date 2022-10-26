package bdmmprime.util.priors;

import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.inference.Distribution;
import beast.base.inference.State;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.TraitSet;
import beast.base.evolution.tree.TreeInterface;

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
   public Input<TraitSet> earlierBoundInput = new Input<>("earlierBound",
            "Initial tip dates", Input.Validate.REQUIRED);
    public Input<TraitSet> laterBoundInput = new Input<>("laterBound",
            "Initial tip dates", Input.Validate.REQUIRED);

    TreeInterface tree;
    TraitSet earlierBound, laterBound;

    Function fso;

    /**
     * Final sample offset for the ages read from trait set representing
     * early bounds of tip times.
     */
    double earlyOffset;

    @Override
    public void initAndValidate() {
        tree = treeInput.get();

        earlierBound = earlierBoundInput.get();
        laterBound = laterBoundInput.get();

        earlyOffset = laterBound.getDate(0) - earlierBound.getDate(0);

        fso = finalSampleOffsetInput.get();
    }

    @Override
    public double calculateLogP() {
        logP = 0.0;

        for (int nr=0; nr<tree.getLeafNodeCount(); nr++) {
            Node node = tree.getNode(nr);
            double nodeAge = node.getHeight() + fso.getArrayValue();
            double earlyAge = earlierBound.getValue(node.getID()) + earlyOffset;
            double lateAge = laterBound.getValue(node.getID());

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
    public void sample(State state, Random random) { }
}
