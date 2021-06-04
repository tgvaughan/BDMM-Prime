package bdmmprime.util.operators;

import beast.core.Description;
import beast.core.Input;
import beast.core.Operator;
import beast.core.parameter.RealParameter;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

@Description("Operator that allows ")
public class TipDateOperator extends Operator {

    public Input<Tree> treeInput = new Input<>("tree",
            "Tree on which to operate.",
            Input.Validate.REQUIRED);

    public Input<TaxonSet> taxonSetInput = new Input<>("taxonSet",
            "Taxon set containing taxa specifying leaves on which to operate");

    public Input<RealParameter> finalSampleOffsetInput = new Input<>("finalSampleOffset",
            "Final sample offset parameter from BD model.",
            Input.Validate.REQUIRED);

    public Input<Double> windowSizeInput = new Input<>("windowSize",
            "Width of window (in time) used to draw new node times from.",
            1.0);

    Tree tree;
    List<String> taxaNames;

    RealParameter finalSampleOffset;
    double windowSize;

    @Override
    public void initAndValidate() {
        tree = treeInput.get();

        if (taxonSetInput.get() != null)
            taxaNames = new ArrayList<>(taxonSetInput.get().getTaxaNames());
        else
            taxaNames = Arrays.asList(tree.getTaxaNames());

        finalSampleOffset = finalSampleOffsetInput.get();
        windowSize = windowSizeInput.get();
    }

    @Override
    public double proposal() {

        // Choose leaf on which to operate
        String taxonID = taxaNames.get(Randomizer.nextInt(taxaNames.size()));

        Node node = tree.getExternalNodes().stream().
                filter(n->n.getID().equals(taxonID)).findFirst().get();


        double oldHeight = node.getHeight();

        // Select new height:
        double newHeight = oldHeight + windowSize*(Randomizer.nextDouble() - 0.5);

        if (newHeight > node.getParent().getHeight())
            // Reflect back down off of parent height
            newHeight = node.getParent().getHeight() - (newHeight - node.getParent().getHeight());

        node.setHeight(newHeight);

        if (node.getHeight() < 0.0 || oldHeight == 0.0)
            if (updateHeights())
                return Double.NEGATIVE_INFINITY;

        return 0;
    }

    /**
     * Update node heights and final sample offset such that the smallest
     * node height is 0.
     *
     * @return false if the operation succeeds, or true if it exceed the
     * bounds of finalSampleOffset.
     */
    private boolean updateHeights() {

        // Find height of youngest leaf:
        double lowestHeight = tree.getExternalNodes().stream()
                .mapToDouble(Node::getHeight).summaryStatistics().getMin();

        if (lowestHeight == 0.0)
            return false;

        double newFSO = finalSampleOffset.getValue() + lowestHeight;

        if (newFSO < Math.max(0.0, finalSampleOffset.getLower())
            || newFSO > finalSampleOffset.getUpper())
            return true;

        tree.getExternalNodes().forEach(n -> n.setHeight(n.getHeight() - lowestHeight));
        finalSampleOffset.setValue(newFSO);

        return false;
    }
}
