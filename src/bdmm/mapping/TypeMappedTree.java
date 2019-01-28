package bdmm.mapping;

import bdmm.parameterization.Parameterization;
import bdmm.util.Utils;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;
import beast.evolution.tree.TraitSet;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;
import org.apache.commons.math3.ode.ContinuousOutputModel;
import org.apache.commons.math3.ode.FirstOrderIntegrator;
import org.apache.commons.math3.ode.nonstiff.DormandPrince54Integrator;

import java.io.PrintStream;

/**
 * An instance of this class is a tree equivalent to untypedTree but with
 * ancestral type changes mapped according the the given multi-type birth-death
 * model.
 *
 * Note that there is a degree of duplication between the code in this class
 * and the code in BirthDeathMigrationDistribution.  Most of this is intentional:
 * the likelihood class cares a lot more about making sure the likelihood calculations
 * are accurate for large data sets, while here we avoid using SmallNumbers and
 * instead rely on dynamic scaling of integration results to prevent underflow.
 * This seems to be good enough for our purpose and allows us to simplify the
 * backwards integration stage of the SM algorithm.  This is important here
 * because unlike the likelihood computation, the SM algorithm requires recording
 * all intermediate results of the backward integration stage for use in the
 * subsequent forward-time simulation stage.
 */
public class TypeMappedTree extends Tree {

    public Input<Parameterization> parameterizationInput = new Input<>("parameterization",
            "BDMM parameterization",
            Input.Validate.REQUIRED);

    public Input<RealParameter> frequenciesInput = new Input<>("frequencies",
            "The frequencies for each type",
            Input.Validate.REQUIRED);

    public Input<TraitSet> typeTraitSetInput = new Input<>("typeTraitSet",
            "Trait information for initializing traits " +
                    "(like node types/locations) in the tree. If this is " +
                    "not provided we will try to extract this information from " +
                    "metadata on the untyped tree leaves.");

    public Input<String> typeLabelInput = new Input<>("typeLabel",
            "Type label used for traits in generated metadata.",
            Input.Validate.REQUIRED);

    public Input<Tree> treeInput = new Input<>("untypedTree",
            "Tree on which to apply mapping.",
            Input.Validate.REQUIRED);

    public Input<Boolean> remapOnLogInput = new Input<>("remapOnLog",
            "If true, mapping will be regenerated when this object " +
                    "is logged.", false);

    private Parameterization param;
    private Tree untypedTree;

    private ODESystem odeSystem;
    private ContinuousOutputModel[] integrationResults;
    double[] geScaleFactors;
    private FirstOrderIntegrator odeIntegrator;

    /**
     * Parameters for backward-time numerical integration.
     */
    private final double BACKWARD_INTEGRATION_MIN_STEP = 1e-100;
    private final double BACKWARD_INTEGRATION_MAX_STEP = 0.1;
    private final double BACKWARD_INTEGRATION_ABS_TOLERANCE = 1e-100;
    private final double BACKWARD_INTEGRATION_REL_TOLERANCE = 1e-7;
    private final int RATE_CHANGE_CHECKS_PER_EDGE = 100;
    private final double RATE_CHANGE_CHECK_CONVERGENCE = 1e-5;
    private final int RATE_CHANGE_MAX_ITERATIONS = 1000;

    /**
     * Maximum number of steps in each waiting time calculation in
     * forward simulation.
     */
    private final int FORWARD_INTEGRATION_STEPS = 100;

    @Override
    public void initAndValidate() {

        param = parameterizationInput.get();
        untypedTree = treeInput.get();

        doStochasticMapping();
    }

    /**
     * Generate new tree by stochastically mapping type changes on untyped tree.
     * Called both during initialization and at when logging.
     */
    public void doStochasticMapping() {
         // Prepare the backward-time integrator!

        odeIntegrator = new DormandPrince54Integrator(
                param.getTotalProcessLength()*BACKWARD_INTEGRATION_MIN_STEP,
                param.getTotalProcessLength()*BACKWARD_INTEGRATION_MAX_STEP,
                BACKWARD_INTEGRATION_ABS_TOLERANCE, BACKWARD_INTEGRATION_REL_TOLERANCE);

        // Prepare the ODE system and the arrays used to store the backward-time
        // integration results.

        odeSystem = new ODESystem(param);
        integrationResults = new ContinuousOutputModel[untypedTree.getNodeCount()];
        geScaleFactors = new double[untypedTree.getNodeCount()];

        // Perform the backward-time integration.
        double[] y = backwardsIntegrateSubtree(untypedTree.getRoot(), 0.0);

        // Sample starting type

        double[] startTypeProbs = new double[param.getNTypes()];

        for (int type=0; type<param.getNTypes(); type++)
            startTypeProbs[type] = y[type+param.getNTypes()]*frequenciesInput.get().getValue(type);

        // (startTypeProbs are unnormalized: this is okay for randomChoicePDF.)
        int startType = Randomizer.randomChoicePDF(startTypeProbs);

        // Simulate type changes down tree

        Node typedRoot = forwardSimulateSubtree(untypedTree.getRoot(), 0.0 , startType);

        // Ensure internal nodes are numbered correctly.  (Leaf node numbers and
        // labels are matched to those in the untyped tree during the simulation.)
        numberInternalNodesOnSubtree(typedRoot, untypedTree.getLeafNodeCount());

        assignFromWithoutID(new Tree(typedRoot));
    }

    /**
     * Obtain value of trait at leaf node.
     *
     * @param leafNode leaf node at which to obtain trait.
     * @return trait value.
     */
    private int getLeafType(Node leafNode) {
            String nodeTypeName;

            if (typeTraitSetInput.get() != null)
                nodeTypeName = typeTraitSetInput.get().getStringValue(leafNode.getID());
            else {
                Object metaData = leafNode.getMetaData(typeLabelInput.get());

                if (metaData instanceof Double)
                    nodeTypeName = String.valueOf(Math.round((double)metaData));
                else
                    nodeTypeName = metaData.toString();
            }

            return param.getTypeSet().getTypeIndex(nodeTypeName);
    }

    /**
     * Set value of trait at node.
     *
     * @param node node at which to set trait.
     * @param type numeric type index.
     */
    private void setNodeType(Node node, int type) {
        node.setMetaData(typeLabelInput.get(), type);

        node.metaDataString = String.format("%s=\"%s\"",
                typeLabelInput.get(), param.getTypeSet().getTypeName(type));
    }

    private boolean[] rhoSampled = null;
    private int[] rhoSamplingIndex = null;

    private void computeRhoSampledLeafStatus() {

        rhoSampled = new boolean[untypedTree.getLeafNodeCount()];
        rhoSamplingIndex = new int[untypedTree.getLeafNodeCount()];

        for (int nodeNr=0; nodeNr < treeInput.get().getLeafNodeCount(); nodeNr++) {
            double nodeTime = param.getNodeTime(untypedTree.getNode(nodeNr));
            rhoSampled[nodeNr] = false;
            for (double rhoSamplingTime : param.getRhoSamplingTimes()) {
                if (Utils.equalWithPrecision(nodeTime, rhoSamplingTime)) {
                    rhoSampled[nodeNr] = true;
                    rhoSamplingIndex[nodeNr] = param.getIntervalIndex(rhoSamplingTime);
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
    private boolean nodeIsRhoSampled(Node node) {
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
    private int getRhoSamplingInterval(Node node) {
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
    private NodeKind getNodeKind(Node node) {
        if (node.isLeaf())
            return NodeKind.LEAF;

        if (node.isFake())
            return NodeKind.SA;

        return NodeKind.INTERNAL;
    }

    /**
     * Integrate p0 and ge from leaves to root of subtree.  Integration results are
     * stored in the field integrationResults.
     *
     * @param untypedSubtreeRoot root node of untyped subtree
     * @param timeOfSubtreeRootEdgeTop time of top of edge above subtree
     * @return integration state at
     */
    double[] backwardsIntegrateSubtree(Node untypedSubtreeRoot,
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

        // Prepare the integrator:

        odeIntegrator.clearEventHandlers();
        odeIntegrator.clearStepHandlers();

        odeIntegrator.addStepHandler(results);

        double delta = 2*Utils.globalPrecisionThreshold;

        double timeOfSubtreeRootEdgeBottom = param.getNodeTime(untypedSubtreeRoot);

        odeIntegrator.addEventHandler(odeSystem,
                (timeOfSubtreeRootEdgeTop-timeOfSubtreeRootEdgeBottom)/RATE_CHANGE_CHECKS_PER_EDGE,
                RATE_CHANGE_CHECK_CONVERGENCE, RATE_CHANGE_MAX_ITERATIONS);

        odeSystem.setInterval(param.getIntervalIndex(timeOfSubtreeRootEdgeBottom-delta));

        // Perform the integration:

        odeIntegrator.integrate(odeSystem,
                timeOfSubtreeRootEdgeBottom - delta, y,
                timeOfSubtreeRootEdgeTop+delta, y);

        // Save integration results
        integrationResults[untypedSubtreeRoot.getNr()] = results;

        return y;
    }

    private double[] getLeafState(Node leafNode) {

        double[] y = new double[param.getNTypes()*2];
        for (int type = 0; type< param.getNTypes(); type++) {
            y[type] = 1.0;
            y[param.getNTypes()+type] = 0.0;
        }

        double leafTime = param.getNodeTime(leafNode);
        double T = param.getTotalProcessLength();

        if (Utils.lessThanWithPrecision(leafTime, T)) {

            int finalInterval = param.getIntervalIndex(T);
            for (int type = 0; type< param.getNTypes(); type++) {
                y[type] *= 1.0 - param.getRhoValues()[finalInterval][type];
            }

            double delta = 2*Utils.globalPrecisionThreshold;

            odeSystem.setInterval(finalInterval);

            odeIntegrator.clearStepHandlers();
            odeIntegrator.clearEventHandlers();
            odeIntegrator.addEventHandler(odeSystem, (T-leafTime)/100, 1e-5, 1000);
            odeIntegrator.integrate(odeSystem, T-delta, y, leafTime+delta, y);
        }

        int leafType = getLeafType(leafNode);

        if (nodeIsRhoSampled(leafNode)) {

            int rhoSamplingInterval = getRhoSamplingInterval(leafNode);

            for (int type = 0; type< param.getNTypes(); type++) {
                double rho = param.getRhoValues()[rhoSamplingInterval][type];
                y[type] *= 1.0 - rho;
                y[type + param.getNTypes()] =
                        type==leafType
                                ? rho
                                : 0.0;
            }

        } else {

            int nodeInterval = param.getNodeIntervalIndex(leafNode);

            for (int type = 0; type< param.getNTypes(); type++) {
                double psi = param.getSamplingRates()[nodeInterval][type];
                double r = param.getRemovalProbs()[nodeInterval][type];

                y[type + param.getNTypes()] =
                        type==leafType
                                ? psi*(r + (1.0-r)*y[type])
                                : 0.0;
            }
        }

        // Scale ge and record scale factor
        geScaleFactors[leafNode.getNr()] = rescaleState(y, 0.0);

        return y;
    }

    private double[] getSAState(Node saNode) {

        double saNodeTime = param.getNodeTime(saNode);

        double[] y = backwardsIntegrateSubtree(saNode.getNonDirectAncestorChild(), saNodeTime);

        int saType = getLeafType(saNode.getDirectAncestorChild());

        if (nodeIsRhoSampled(saNode)) {

            int rhoSamplingInterval = getRhoSamplingInterval(saNode);

            for (int type = 0; type< param.getNTypes(); type++) {
                double rho = param.getRhoValues()[rhoSamplingInterval][type];
                double r = param.getRemovalProbs()[rhoSamplingInterval][type];

                y[type] *= 1.0 - rho;
                y[type+ param.getNTypes()] *=
                        type==saType
                                ? rho*(1-r)
                                : 0.0;
            }

        } else {

            int nodeInterval = param.getNodeIntervalIndex(saNode);

            for (int type = 0; type< param.getNTypes(); type++) {
                double psi = param.getSamplingRates()[nodeInterval][type];
                double r = param.getRemovalProbs()[nodeInterval][type];

                y[type + param.getNTypes()] *=
                        type==saType
                                ? psi*(1-r)
                                : 0.0;
            }
        }

        // Scale ge and record scale factor
        geScaleFactors[saNode.getNr()] = rescaleState(y, geScaleFactors[saNode.getNonDirectAncestorChild().getNr()]);

        return y;
    }

    private double[] getInternalState(Node internalNode) {

        double internalNodeTime = param.getNodeTime(internalNode);

        double[] yLeft = backwardsIntegrateSubtree(internalNode.getChild(0), internalNodeTime);
        double[] yRight = backwardsIntegrateSubtree(internalNode.getChild(1), internalNodeTime);

        double logFLeft = geScaleFactors[internalNode.getChild(0).getNr()];
        double logFRight = geScaleFactors[internalNode.getChild(1).getNr()];

        double[] y = new double[param.getNTypes()*2];

        int nodeInterval = param.getNodeIntervalIndex(internalNode);

        int N = param.getNTypes();

        for (int type = 0; type< param.getNTypes(); type++) {
            y[type] = yLeft[type];
            y[N+type] = 0.0;

            for (int typeOther=0; typeOther<N; typeOther++) {
                if (typeOther == type) {
                    y[N+type] += param.getBirthRates()[nodeInterval][type]
                            *yLeft[N+type]*yRight[N+type];
                } else {
                    y[N+type] += param.getCrossBirthRates()[nodeInterval][type][typeOther]
                            *(yLeft[N+type]*yRight[N+typeOther] + yLeft[N+typeOther]*yRight[N+type]);
                }
            }
        }

        // Scale ge and record scale factor
        geScaleFactors[internalNode.getNr()] = rescaleState(y, logFLeft+logFRight);

        return y;
    }

    /**
     * Scale state vector y so that its largest element is 1.0.  Takes
     * the existing scale factor as input and returns the updated scale factor.
     *
     * @param y state vector to scale
     * @param prevLogF current (log) scaling factor
     * @return updated (log) scaling factor.
     */
    private double rescaleState(double[] y, double prevLogF) {

        double C = 0.0;
        for (int type = 0; type< param.getNTypes(); type++)
            C = Math.max(C, y[type+ param.getNTypes()]);

        for (int type = 0; type< param.getNTypes(); type++)
            y[type+ param.getNTypes()] /= C;

        return prevLogF + Math.log(C);
    }


    /**
     * Use a forward-time stochastic simulation algorithm to apply type changes
     * to a tree.  This requires that the backwards integration calculation has
     * already been performed.
     *
     * @param subtreeRoot root of (untyped) subtree to generate mapping for.
     * @param startTime time above root to start the simulation.
     * @param startType type at the start of the simulation.
     * @return root of new tree with type changes marked.
     */
    Node forwardSimulateSubtree (Node subtreeRoot, double startTime, int startType) {

        Node root = new Node();
        setNodeType(root, startType);

        Node currentNode = root;
        int currentType = startType;
        double currentTime = startTime;

        double endTime = param.getNodeTime(subtreeRoot);

        double[] rates = new double[param.getNTypes()];
        double[] ratesPrime = new double[param.getNTypes()];
        double totalRate, totalRatePrime;

        while (true) {

            // Determine time of next event

            double K = -Math.log(Randomizer.nextDouble());
            double I = 0.0;

            double t = currentTime;
            double dt = (endTime-currentTime)/ FORWARD_INTEGRATION_STEPS;
            totalRate = getTotalFowardsRate(currentType, currentTime, subtreeRoot, rates);

            int integrationStep;
            for (integrationStep=0; integrationStep< FORWARD_INTEGRATION_STEPS; integrationStep++) {
                double tprime = currentTime + (endTime-currentTime)*(integrationStep+1)/ FORWARD_INTEGRATION_STEPS;

                totalRatePrime = getTotalFowardsRate(currentType, tprime, subtreeRoot, ratesPrime);

                I += dt*(totalRate + totalRatePrime)/2.0;

                if (I >= K) {
                    currentTime = t + 0.5*dt;
                    break;
                }

                totalRate = totalRatePrime;
                double[] tmp = rates;
                rates = ratesPrime;
                ratesPrime = tmp;

                t = tprime;
            }

            if (integrationStep == FORWARD_INTEGRATION_STEPS)
                break;

            currentNode.setHeight(param.getTotalProcessLength() - currentTime);

            // Sample event type

            currentType = Randomizer.randomChoicePDF(ratesPrime);

            // Implement event in tree

            Node newNode = new Node();
            setNodeType(newNode, currentType);

            currentNode.addChild(newNode);
            currentNode = newNode;
        }

        currentNode.setHeight(param.getTotalProcessLength() - endTime);

        switch(getNodeKind(subtreeRoot)) {
            case LEAF:
                currentNode.setNr(subtreeRoot.getNr());
                currentNode.setID(subtreeRoot.getID());
                break;

            case SA:
                currentNode.addChild(forwardSimulateSubtree(subtreeRoot.getNonDirectAncestorChild(), endTime, currentType));
                break;

            case INTERNAL:
                // TODO Add support for birth among demes.
                currentNode.addChild(forwardSimulateSubtree(subtreeRoot.getChild(0), endTime, currentType));
                currentNode.addChild(forwardSimulateSubtree(subtreeRoot.getChild(1), endTime, currentType));
                break;

            default:
                throw new IllegalArgumentException("Switch fell through in forward simulation.");
        }

        return root;
    }

    private double[] getForwardsRates(int fromType, double time, Node baseNode, double[] result) {
        ContinuousOutputModel com = integrationResults[baseNode.getNr()];
        com.setInterpolatedTime(time);
        double[] y = com.getInterpolatedState();
        int interval = param.getIntervalIndex(time);

        for (int type=0; type<param.getNTypes(); type++) {
            if (type == fromType) {
                result[type] = 0.0;
                continue;
            }

            result[type] = (param.getCrossBirthRates()[interval][fromType][type] * y[type]
                        + param.getMigRates()[interval][fromType][type])
                        * y[param.getNTypes() + type];
        }

        if (y[param.getNTypes()+fromType]<=0.0) {
            // The source type prob approaches zero as the integration closes
            // in on a node with a defined type.  This causes the transition
            // rate to this type to become infinite.  What follows is a hack
            // to ensure that this important situation is handled properly.

            int maxRateIdx = -1;
            double maxRate = Double.NEGATIVE_INFINITY;
            for (int type=0; type<param.getNTypes(); type++) {
                if (result[type]>maxRate) {
                    maxRateIdx = type;
                    maxRate = result[type];
                }
            }

            for (int type=0; type<param.getNTypes(); type++)
                result[type] = type == maxRateIdx ? 1.0 : 0.0 ;

        } else {
            // Apply source type prob as rate denominator:

            for (int type=0; type<param.getNTypes(); type++)
                result[type] /= y[param.getNTypes()+fromType];

        }

        return result;
    }

    private double getTotalFowardsRate(int fromType, double time, Node baseNode, double[] working) {
        double totalRate = 0.0;
        double[] rates = getForwardsRates(fromType, time, baseNode, working);
        for (int type=0; type<param.getNTypes(); type++)
            totalRate += rates[type];

        return totalRate;
    }

    /**
     * Apply node numbers to internal nodes below and including subtreeRoot.
     * Numbers are applied postorder, so parents always have larger numbers
     * than their children and the root has the hightest number.
     *
     * @param subtreeRoot root of subtree
     * @param nextNumber next number to be used
     * @return next number to be used on another part of the tree.
     */
    private int numberInternalNodesOnSubtree(Node subtreeRoot, int nextNumber) {

        if (subtreeRoot.isLeaf())
            return nextNumber;

        for (Node child : subtreeRoot.getChildren())
            nextNumber = numberInternalNodesOnSubtree(child, nextNumber);

        subtreeRoot.setNr(nextNumber);

        return nextNumber+1;
    }

    /*
     * Loggable implementation
     */

    @Override
    public void log(long sample, PrintStream out) {

        if (remapOnLogInput.get())
            doStochasticMapping();

        Tree tree = (Tree) getCurrent();
        out.print("tree STATE_" + sample + " = ");
        // Don't sort, this can confuse CalculationNodes relying on the tree
        //tree.getRoot().sort();
        final int[] dummy = new int[1];
        final String newick = tree.getRoot().toSortedNewick(dummy, true);
        out.print(newick);
        out.print(";");
    }
}
