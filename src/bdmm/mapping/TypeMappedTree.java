package bdmm.mapping;

import bdmm.parameterization.Parameterization;
import bdmm.util.Utils;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;
import beast.evolution.tree.TraitSet;
import beast.evolution.tree.Tree;
import org.apache.commons.math3.ode.FirstOrderIntegrator;
import org.apache.commons.math3.ode.nonstiff.DormandPrince54Integrator;

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
    Tree untypedTree;

    ODESystem odeSystem;

    @Override
    public void initAndValidate() {

        parameterization = parameterizationInput.get();
        untypedTree = treeInput.get();

        odeSystem = new ODESystem(parameterization);

        backwardsIntegrateSubtree(untypedTree.getRoot(), 0.0);

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

    boolean[] rhoSampled = null;

    /**
     * Return true if node is the result of a rho sampling event.
     *
     * @param node node about which to make query
     * @return true if node time coincides with rho sampling time.
     */
    public boolean nodeIsRhoSampled(Node node) {
        if (rhoSampled == null) {

            rhoSampled = new boolean[untypedTree.getLeafNodeCount()];

            for (int nodeNr=0; nodeNr < treeInput.get().getLeafNodeCount(); nodeNr++) {
                double nodeTime = parameterization.getNodeTime(untypedTree.getNode(nodeNr));
                for (double rhoSamplingTime : parameterization.getRhoSamplingTimes())
                    if (Utils.equalWithPrecision(nodeTime, rhoSamplingTime))
                        rhoSampled[nodeNr] = true;
            }
        }

        return rhoSampled[node.getNr()];
    }

    /**
     * Type specifying different kinds of nodes.
     */
    enum NodeKind {LEAF, SA, INTERNAL};

    /**
     * Return "kind" of node.
     *
     * @param node node to classify
     * @return node kind.
     */
    NodeKind getNodeKind(Node node) {
        if (node.isLeaf())
            return NodeKind.LEAF;

        if (node.isFake())
            return NodeKind.SA;

        return NodeKind.INTERNAL;
    }


    public void backwardsIntegrateSubtree (Node untypedSubtreeRoot,
                                           double timeOfSubtreeRootEdgeTop) {

        double[] y;

        switch(getNodeKind(untypedSubtreeRoot)) {
            case LEAF:
                y = getLeafState(untypedSubtreeRoot);
                break;

            case SA:
                y = getSAState(untypedSubtreeRoot);
                break;

            case INTERNAL:
                y = getInternalState(untypedSubtreeRoot);
                break;

            default:
                throw new RuntimeException("Node kind switch fell through!");
        }

        double rootTime = parameterization.getNodeTime(untypedSubtreeRoot);
        double edgeLength = timeOfSubtreeRootEdgeTop - rootTime;

        FirstOrderIntegrator integrator = new DormandPrince54Integrator(
                parameterization.getTotalProcessLength()/1e100,
                parameterization.getTotalProcessLength()/10.0,
                1e-100, 1e-7);

        double delta = 2*Utils.globalPrecisionThreshold;

        integrator.integrate(odeSystem,
                parameterization.getNodeTime(untypedSubtreeRoot) - delta, y,
                timeOfSubtreeRootEdgeTop+delta, y);

    }

    double[] getLeafState(Node leaf) {



        return null;
    }

    double[] getSAState(Node saNode) {
        return null;
    }

    double[] getInternalState(Node internalNode) {
        return null;
    }


    public Node fowardSimulation() {
       Node root = new Node();

       return root;
    }

}
