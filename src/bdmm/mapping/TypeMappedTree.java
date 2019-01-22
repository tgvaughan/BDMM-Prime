package bdmm.mapping;

import bdmm.parameterization.Parameterization;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;
import beast.evolution.tree.TraitSet;
import beast.evolution.tree.Tree;

public class TypeMappedTree extends Tree {

    public Input<Parameterization> parameterizationInput = new Input<>("parameterization",
            "BDMM parameterization",
            Input.Validate.REQUIRED);

    public Input<RealParameter> frequenciesInput = new Input<>("frequencies",
            "The frequencies for each type",
            Input.Validate.REQUIRED);

    public Input<TraitSet> typeTraitSetInput = new Input<>("typeTraitSet",
            "Trait information for initializing traits " +
                    "(like node types/locations) in the tree");

    public Input<String> typeLabelInput = new Input<>("typeLabel",
            "type label in tree for initializing traits " +
                    "(like node types/locations) in the tree",
            Input.Validate.XOR, typeTraitSetInput);

    public Input<Double> relativeToleranceInput = new Input<>("relTolerance",
            "relative tolerance for numerical integration",
            1e-7);

    public Input<Double> absoluteToleranceInput = new Input<>("absTolerance",
            "absolute tolerance for numerical integration",
            1e-100 /*Double.MIN_VALUE*/);

    public Input<Tree> treeInput = new Input<>("untypedTree",
            "Tree on which to apply mapping.",
            Input.Validate.REQUIRED);

    Parameterization parameterization;

    @Override
    public void initAndValidate() {

        parameterization = parameterizationInput.get();
        Node untypedRoot = treeInput.get().getRoot();

        backwardsIntegrateSubtree(untypedRoot, parameterization.getTotalProcessLength());

//        Node typedRoot = fowardSimulation(treeInput.get().getRoot(), rootType);
//
//        assignFromWithoutID(new Tree(typedRoot));

        super.initAndValidate();
    }

    public int getLeafType(Node leafNode) {
        if (typeTraitSetInput.get() != null)
            return (int) typeTraitSetInput.get().getValue(leafNode.getID());
        else {
            Object metaData = leafNode.getMetaData(typeLabelInput.get());
            if (metaData instanceof Double)
                return (int) Math.round((double) metaData);
            else if (metaData instanceof Integer)
                return (int) metaData;
            else if (metaData instanceof String)
                return Integer.valueOf((String) metaData);
            else
                throw new IllegalArgumentException(
                        "Cannot determine type of taxon '" +
                                leafNode.getID() + "'.");
        }
    }

    enum NodeKind {LEAF, SA, INTERNAL};
    NodeKind getNodeKind(Node node) {
        if (node.isLeaf())
            return NodeKind.LEAF;

        if (node.isFake())
            return NodeKind.SA;

        return NodeKind.INTERNAL;
    }

    public boolean isRhoSampled(Node node) {
        if node.getHeight()
    }

    public void backwardsIntegrateSubtree (Node untypedSubtreeRoot,
                                           double ageOfSubtreeRootEdgeTop) {

        switch(getNodeKind(untypedSubtreeRoot)) {
            case LEAF:

            case SA:

            case INTERNAL:

            default:
                throw new RuntimeException("Node kind switch fell through!");
        }


    }

    public Node fowardSimulation() {
       Node root = new Node();

       return root;
    }

}
