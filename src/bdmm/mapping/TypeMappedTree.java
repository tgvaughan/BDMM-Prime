package bdmm.mapping;

import bdmm.parameterization.Parameterization;
import bdmm.util.Utils;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;
import beast.evolution.tree.TraitSet;
import beast.evolution.tree.Tree;
import org.apache.commons.math3.ode.ContinuousOutputModel;
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
    ContinuousOutputModel[] integrationResults;

    @Override
    public void initAndValidate() {

        parameterization = parameterizationInput.get();
        untypedTree = treeInput.get();

        odeSystem = new ODESystem(parameterization);
        integrationResults = new ContinuousOutputModel[untypedTree.getNodeCount()];
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
    int[] rhoSamplingIndex = null;

    public void computeRhoSampledLeafStatus() {

        rhoSampled = new boolean[untypedTree.getLeafNodeCount()];
        rhoSamplingIndex = new int[untypedTree.getLeafNodeCount()];

        for (int nodeNr=0; nodeNr < treeInput.get().getLeafNodeCount(); nodeNr++) {
            double nodeTime = parameterization.getNodeTime(untypedTree.getNode(nodeNr));
            rhoSampled[nodeNr] = false;
            for (double rhoSamplingTime : parameterization.getRhoSamplingTimes()) {
                if (Utils.equalWithPrecision(nodeTime, rhoSamplingTime)) {
                    rhoSampled[nodeNr] = true;
                    rhoSamplingIndex[nodeNr] = parameterization.getIntervalIndex(rhoSamplingTime);
                    break;
                }
            }
        }
    }

    /**
     * Return true if node is the result of a rho sampling event.
     *
     * @param node node about which to make query
     * @return true if node time coincides with rho sampling time.
     */
    public boolean nodeIsRhoSampled(Node node) {
        if (rhoSampled == null)
            computeRhoSampledLeafStatus();

        return rhoSampled[node.getNr()];
    }

    /**
     * When node is rho sampled, return index of interval whose end time
     * corresponds to this rho sampling event.
     *
     * @param node node about which to make query
     * @return interval index
     */
    public int getRhoSamplingInterval(Node node) {
        if (nodeIsRhoSampled(node))
            return rhoSamplingIndex[node.getNr()];

        throw new IllegalArgumentException("Node is not rho sampled.");
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

    FirstOrderIntegrator getNewIntegrator() {
        return new DormandPrince54Integrator(
                parameterization.getTotalProcessLength()/1e100,
                parameterization.getTotalProcessLength()/10.0,
                1e-100, 1e-7);
    }


    public double[] backwardsIntegrateSubtree (Node untypedSubtreeRoot,
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


        ContinuousOutputModel results = new ContinuousOutputModel();

        FirstOrderIntegrator integrator = getNewIntegrator();
        integrator.addStepHandler(results);

        double delta = 2*Utils.globalPrecisionThreshold;

        integrator.integrate(odeSystem,
                parameterization.getNodeTime(untypedSubtreeRoot) - delta, y,
                timeOfSubtreeRootEdgeTop+delta, y);

        // Save integration results
        integrationResults[untypedSubtreeRoot.getNr()] = results;

        return y;
    }

    double[] getLeafState(Node leafNode) {

        double[] y = new double[parameterization.getNTypes()*2];
        for (int type=0; type<parameterization.getNTypes(); type++) {
            y[type] = 1.0;
            y[parameterization.getNTypes()+type] = 0.0;
        }

        double leafTime = parameterization.getNodeTime(leafNode);
        double T = parameterization.getTotalProcessLength();

        if (Utils.lessThanWithPrecision(leafTime, T)) {

            int finalInterval = parameterization.getIntervalIndex(T);
            for (int type=0; type<parameterization.getNTypes(); type++) {
                y[type] *= 1.0 - parameterization.getRhoValues()[finalInterval][type];
            }

            FirstOrderIntegrator integrator = getNewIntegrator();

            double delta = 2*Utils.globalPrecisionThreshold;

            odeSystem.setInterval(finalInterval);
            integrator.integrate(odeSystem, T-delta, y, leafTime+delta, y);
        }

        int leafType = getLeafType(leafNode);

        if (nodeIsRhoSampled(leafNode)) {

            int rhoSamplingInterval = getRhoSamplingInterval(leafNode);

            for (int type=0; type<parameterization.getNTypes(); type++) {
                double rho = parameterization.getRhoValues()[rhoSamplingInterval][type];
                y[type] *= 1.0 - rho;
                y[type + parameterization.getNTypes()] =
                        type==leafType
                                ? Math.log(rho)
                                : Double.NEGATIVE_INFINITY;
            }

        } else {

            int nodeInterval = parameterization.getNodeIntervalIndex(leafNode);

            for (int type=0; type<parameterization.getNTypes(); type++) {
                double psi = parameterization.getSamplingRates()[nodeInterval][type];
                double r = parameterization.getRemovalProbs()[nodeInterval][type];

                y[type + parameterization.getNTypes()] =
                        type==leafType
                                ? Math.log(psi*(r + (1.0-r)*y[type]))
                                : Double.NEGATIVE_INFINITY;
            }
        }

        return y;
    }

    double[] getSAState(Node saNode) {

        double saNodeTime = parameterization.getNodeTime(saNode);

        double[] y = backwardsIntegrateSubtree(saNode.getNonDirectAncestorChild(), saNodeTime);

        int saType = getLeafType(saNode.getDirectAncestorChild());

        if (nodeIsRhoSampled(saNode)) {

            int rhoSamplingInterval = getRhoSamplingInterval(saNode);

            for (int type=0; type<parameterization.getNTypes(); type++) {
                double rho = parameterization.getRhoValues()[rhoSamplingInterval][type];
                double r = parameterization.getRemovalProbs()[rhoSamplingInterval][type];

                y[type] *= 1.0 - rho;
                y[type+parameterization.getNTypes()] +=
                        type==saType
                                ? Math.log(rho*(1-r))
                                : Double.NEGATIVE_INFINITY;
            }

        } else {

            int nodeInterval = parameterization.getNodeIntervalIndex(saNode);

            for (int type=0; type<parameterization.getNTypes(); type++) {
                double psi = parameterization.getSamplingRates()[nodeInterval][type];
                double r = parameterization.getRemovalProbs()[nodeInterval][type];

                y[type + parameterization.getNTypes()] +=
                        type==saType
                                ? Math.log(psi*(1-r))
                                : Double.NEGATIVE_INFINITY;
            }
        }

        return y;
    }

    double[] getInternalState(Node internalNode) {

        double internalNodeTime = parameterization.getNodeTime(internalNode);

        double[] yLeft = backwardsIntegrateSubtree(internalNode.getChild(0), internalNodeTime);
        double[] yRight = backwardsIntegrateSubtree(internalNode.getChild(1), internalNodeTime);

        double[] y = new double[parameterization.getNTypes()*2];

        int nodeInterval = parameterization.getNodeIntervalIndex(internalNode);

        int N = parameterization.getNTypes();

        for (int type=0; type<parameterization.getNTypes(); type++) {
            y[type] = yLeft[type];


            double scaleFactor = Double.NEGATIVE_INFINITY;

            for (int typeOther=0; typeOther<N; typeOther++) {
                if (typeOther == type) {
                    scaleFactor = Math.max(scaleFactor, yLeft[type+N] + yRight[type+N]);
                } else {
                    scaleFactor = Math.max(scaleFactor,
                            Math.max(yLeft[type+N] + yRight[typeOther+N],
                                    yLeft[typeOther+N] + yRight[type+N]));

                }
            }

            double scaledSum = 0;

            for (int typeOther=0; typeOther<N; typeOther++) {
                if (typeOther == type) {
                    scaledSum += parameterization.getBirthRates()[nodeInterval][type]
                            * Math.exp(yLeft[type+N] + yRight[type+N] - scaleFactor);
                } else {
                    scaledSum += parameterization.getCrossBirthRates()[nodeInterval][type][typeOther]
                            * (Math.exp(yLeft[type+N] + yRight[typeOther+N] - scaleFactor)
                            + Math.exp(yLeft[typeOther+N] + yRight[type+N] - scaleFactor));
                }
            }

            y[type+N] = Math.log(scaledSum) + scaleFactor;
        }

        return y;
    }


    public Node fowardSimulation() {
       Node root = new Node();

       return root;
    }

}
