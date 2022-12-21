package bdmmprime.util.priors;

import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.core.Log;
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

    public Input<Function> endOfSamplingTimeInput = new Input<>("endOfSamplingTime",
            "Time of the point when sampling ends.  (Necessary only " +
                    "when upper and lower bounds are given forward in time.)");

    public Input<Boolean> reportBoundsViolationsInput = new Input<>(
            "reportBoundsViolations",
            "Causes the distribution to report which taxon exceeded " +
                    "its bounds in each case.  Useful for diagnosing " +
                    "initialization problems.", false);

    TreeInterface tree;
    TraitSet earlierBound, laterBound;
    boolean boundsAreAges;

    Function fso, endOfSamplingTime;

    boolean reportBoundsViolations;

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

        boundsAreAges = !earlierBound.isDateTrait()
                || earlierBound.getDateType().equals(TraitSet.AGE_TRAIT)
                || earlierBound.getDateType().equals(TraitSet.DATE_BACKWARD_TRAIT);

        if (boundsAreAges && laterBound.isDateTrait()
                && (laterBound.getDateType().equals(TraitSet.DATE_TRAIT)
                || laterBound.getDateType().equals(TraitSet.DATE_FORWARD_TRAIT)))
            throw new IllegalArgumentException("earlierBound and " +
                    "laterBound trait sets must both be forward in time " +
                    "or both backward in time (ages), not a mixture.");

        fso = finalSampleOffsetInput.get();
        endOfSamplingTime = endOfSamplingTimeInput.get();

        if (!boundsAreAges && endOfSamplingTime == null)
            throw new IllegalArgumentException("If bounds are given forward " +
                    "in time, you must also provide a value to the " +
                    "endOfSamplingTime input.");

        reportBoundsViolations = reportBoundsViolationsInput.get();
    }

    double getBoundAge(TraitSet boundTrait, Node node) {
        if (boundsAreAges)
            return boundTrait.getValue(node.getID()) + boundTrait.getDate(0);
        else
            return boundTrait.getValue(node.getID()) +
                    (endOfSamplingTime.getArrayValue() - boundTrait.getDate(0));
    }

    @Override
    public double calculateLogP() {
        logP = 0.0;

        for (int nr=0; nr<tree.getLeafNodeCount(); nr++) {
            Node node = tree.getNode(nr);
            double nodeAge = node.getHeight() + fso.getArrayValue();
            double earlyAge = getBoundAge(earlierBound, node);
            double lateAge = getBoundAge(laterBound, node);

            if (nodeAge > earlyAge || nodeAge < lateAge) {
                if (reportBoundsViolations) {
                    Log.err.println("Taxon " + node.getID() +
                            " (" + nr + ") has an age of " + nodeAge +
                            "which is outside the allowed range of [" +
                            lateAge + "," + earlyAge + "].");
                }
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
