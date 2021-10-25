package bdmmprime.util.operators;

import beast.core.Description;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.operators.TreeOperator;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

@Description("Operator that allows ")
public class TipEdgeScaleOperator extends TreeOperator {

    public Input<TaxonSet> taxonSetInput = new Input<>("taxonSet",
            "Taxon set containing taxa specifying leaves on which to operate");

    public Input<RealParameter> finalSampleOffsetInput = new Input<>("finalSampleOffset",
            "Final sample offset parameter from BD model.",
            Input.Validate.REQUIRED);

    public Input<Double> scaleFactorInput = new Input<>("scaleFactor",
            "Edge length is scaled by a value chosen from [sf,1/sf].",
            0.8);

    Tree tree;
    List<String> taxaNames;

    RealParameter finalSampleOffset;
    double scaleFactor;

    @Override
    public void initAndValidate() {
        tree = treeInput.get();

        if (taxonSetInput.get() != null)
            taxaNames = new ArrayList<>(taxonSetInput.get().getTaxaNames());
        else
            taxaNames = Arrays.asList(tree.getTaxaNames());

        finalSampleOffset = finalSampleOffsetInput.get();
        scaleFactor = scaleFactorInput.get();
    }

    @Override
    public double proposal() {

        // Choose leaf on which to operate
        String taxonID = taxaNames.get(Randomizer.nextInt(taxaNames.size()));

        Node node = tree.getExternalNodes().stream().
                filter(n->n.getID().equals(taxonID)).findFirst().get();

        double oldHeight = node.getHeight();
        double oldEdgeLength = node.getParent().getHeight() - oldHeight;

        double minf = Math.min(scaleFactor, 1.0/scaleFactor);
        double maxf = 1.0/minf;
        double f = minf + (maxf - minf)*Randomizer.nextDouble();

        // Select new height:
        double newEdgeLength = oldEdgeLength*f;
        double newHeight = node.getParent().getHeight() - newEdgeLength;

        node.setHeight(newHeight);

        if (node.getHeight() < 0.0 || oldHeight == 0.0)
            if (updateHeights())
                return Double.NEGATIVE_INFINITY;

        return -Math.log(f);
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

        for (Node node : tree.getNodesAsArray())
            node.setHeight(node.getHeight() - lowestHeight);

        finalSampleOffset.setValue(newFSO);

        return false;
    }
}
