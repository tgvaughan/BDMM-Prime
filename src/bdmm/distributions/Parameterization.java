package bdmm.distributions;

import beast.core.CalculationNode;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Tree;
import beast.evolution.tree.TreeInterface;

public abstract class Parameterization extends CalculationNode {

    public Input<RealParameter> originInput = new Input<>("origin",
            "Length of the BDMM process.");

    public Input<Tree> treeInput = new Input<>("treeForMRCA",
            "If provided, condition BDMM on time of MRCA of this tree.",
            Input.Validate.XOR, originInput);

    @Override
    public void initAndValidate() { }

    public double getOrigin() {
        return originInput.get().getValue();
    }

    public boolean hasOrigin() {
        return originInput.get() != null;
    }

    public double getMaxTime() {
        if (hasOrigin())
            return getOrigin();
        else
            return treeInput.get().getRoot().getHeight();
    }

    public double getRootEdgeLength(TreeInterface tree) {
        return getMaxTime() - tree.getRoot().getHeight();
    }

    @Override
    protected boolean requiresRecalculation() {
        return true;
    }
}
