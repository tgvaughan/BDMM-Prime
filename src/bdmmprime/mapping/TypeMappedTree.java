/*
 * Copyright (C) 2019-2024 Tim Vaughan
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package bdmmprime.mapping;

import bdmmprime.distribution.BirthDeathMigrationDistribution;
import bdmmprime.parameterization.Parameterization;
import bdmmprime.util.Utils;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.TraitSet;
import beast.base.evolution.tree.Tree;
import beast.base.inference.parameter.RealParameter;
import beast.base.util.Randomizer;
import org.apache.commons.math3.ode.ContinuousOutputModel;
import org.apache.commons.math3.ode.FirstOrderIntegrator;
import org.apache.commons.math3.ode.nonstiff.DormandPrince54Integrator;

import java.io.PrintStream;

/**
 * <p>An instance of this class is a tree equivalent to untypedTree but with
 * ancestral type changes mapped according the the given multi-type birth-death
 * model.</p>
 *
 * <p>Note that there is a degree of duplication between the code in this class
 * and the code in BirthDeathMigrationDistribution.  Most of this is intentional:
 * the likelihood class cares a lot more about making sure the likelihood calculations
 * are accurate for large data sets, while here we avoid using SmallNumbers and
 * instead rely on dynamic scaling of integration results to prevent underflow.
 * This seems to be good enough for our purpose and allows us to simplify the
 * backwards integration stage of the SM algorithm.  This is important here
 * because unlike the likelihood computation, the SM algorithm requires recording
 * all intermediate results of the backward integration stage for use in the
 * subsequent forward-time simulation stage.</p>
 *
 * <p>Furthermore, note that the code here separately integrates the p_i(t)
 * equations up to each leaf node.  This is inefficient, as the same ODE
 * is integrated over the same time period multiple times.  In practice this
 * probably doesn't matter too much as the mapper is only called rarely and
 * this is a relatively small part of the mapping algorithm.</p>
 *
 * <p>As for the refactored BirthDeathMigrationDistribution class, the backward
 * integration strategy here is to have the integrator only handle rate shift and
 * rho sampling events which are _not_ coincident with nodes in the tree.
 * Events which _are_ coincident are handled as part of the ODE boundary
 * condition calculations done at each node.</p>
 */
public class TypeMappedTree extends Tree {

    public Input<Parameterization> parameterizationInput = new Input<>("parameterization",
            "BDMM parameterization",
            Input.Validate.REQUIRED);

    public Input<RealParameter> frequenciesInput = new Input<>("frequencies",
            "The frequencies for each type",
            Input.Validate.REQUIRED);

    public Input<Function> finalSampleOffsetInput = new Input<>("finalSampleOffset",
            "If provided, the difference in time between the final sample and the end of the BD process.",
            new RealParameter("0.0"));

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

    public Input<Boolean> mapOnInitInput = new Input<>("mapOnInit",
            "If true, mapping will be performed when object is " +
                    "first initialize.", true);

    public Input<BirthDeathMigrationDistribution> bdmmDistribInput = new Input<>("bdmmDistrib",
            "If provided, extract the parameterization from here.",
            Input.Validate.XOR, parameterizationInput);

    private Parameterization param;
    private Function finalSampleOffset;
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

        if (parameterizationInput.get() != null)
            param = parameterizationInput.get();
        else
            param = bdmmDistribInput.get().parameterizationInput.get();

        untypedTree = treeInput.get();

        finalSampleOffset = finalSampleOffsetInput.get();

        if (mapOnInitInput.get())
            doStochasticMapping();
    }

    /**
     * Generate new tree by stochastically mapping type changes on untyped tree.
     * Called both during initialization and at when logging.
     */
    private void doStochasticMapping() {
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

        // Update leaf rho sampling status:
        computeRhoSampledLeafStatus();

        // Perform the backward-time integration.
        double[] y = backwardsIntegrateSubtree(untypedTree.getRoot(), 0.0);

        // Sample starting type

        double[] startTypeProbs = new double[param.getNTypes()];

        for (int type=0; type<param.getNTypes(); type++)
            startTypeProbs[type] = y[type+param.getNTypes()]*frequenciesInput.get().getValue(type);

        // (startTypeProbs are unnormalized: this is okay for randomChoicePDF.)
        int startType = Randomizer.randomChoicePDF(startTypeProbs);

        // Simulate type changes down tree
        // As tiny numerical errors can very occasionally lead to this failing,
        // the mapping is attempted up to 3 times, with additional failures leading
        // to a runtime exception.

        int failures = 0;
        boolean success = false;
        Node typedRoot = null;
        while (!success && failures<3) {
            try {
                typedRoot = forwardSimulateSubtree(untypedTree.getRoot(), 0.0, startType);
                success = true;
            } catch (java.lang.Error ex) {
                System.err.println("Exception encountered in stochastic mapping calculation. Retrying...");
                failures += 1;
            }
        }
        if (!success)
            throw new IllegalStateException("Too many failures (3) encountered attempting stochastic trait mapping.");

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
            double nodeTime = param.getNodeTime(untypedTree.getNode(nodeNr), finalSampleOffset.getArrayValue());
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
                throw new IllegalStateException("Node kind switch fell through. (Should be impossible.)");
        }

        ContinuousOutputModel results = new ContinuousOutputModel();

        // Prepare the integrator:

        odeIntegrator.clearEventHandlers();
        odeIntegrator.clearStepHandlers();

        odeIntegrator.addStepHandler(results);

        double delta = 2*Utils.globalPrecisionThreshold;

        double timeOfSubtreeRootEdgeBottom = param.getNodeTime(untypedSubtreeRoot, finalSampleOffset.getArrayValue());

        odeIntegrator.addEventHandler(odeSystem,
                (timeOfSubtreeRootEdgeBottom-timeOfSubtreeRootEdgeTop)/RATE_CHANGE_CHECKS_PER_EDGE,
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

        double leafTime = param.getNodeTime(leafNode, finalSampleOffset.getArrayValue());
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

            int nodeInterval = param.getNodeIntervalIndex(leafNode, finalSampleOffset.getArrayValue());

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

        double saNodeTime = param.getNodeTime(saNode, finalSampleOffset.getArrayValue());

        double[] y = backwardsIntegrateSubtree(saNode.getNonDirectAncestorChild(), saNodeTime);

        int saType = getLeafType(saNode.getDirectAncestorChild());

        if (nodeIsRhoSampled(saNode.getDirectAncestorChild())) {

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

            int nodeInterval = param.getNodeIntervalIndex(saNode, finalSampleOffset.getArrayValue());

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

        double internalNodeTime = param.getNodeTime(internalNode, finalSampleOffset.getArrayValue());

        double[] yLeft = backwardsIntegrateSubtree(internalNode.getChild(0), internalNodeTime);
        double[] yRight = backwardsIntegrateSubtree(internalNode.getChild(1), internalNodeTime);

        double logFLeft = geScaleFactors[internalNode.getChild(0).getNr()];
        double logFRight = geScaleFactors[internalNode.getChild(1).getNr()];

        double[] y = new double[param.getNTypes()*2];

        int nodeInterval = param.getNodeIntervalIndex(internalNode, finalSampleOffset.getArrayValue());

        int N = param.getNTypes();

        for (int type = 0; type< param.getNTypes(); type++) {
            y[type] = yLeft[type];
            y[N+type] = 0.0;

            for (int typeOther=0; typeOther<N; typeOther++) {
                if (typeOther == type) {
                    y[N+type] += param.getBirthRates()[nodeInterval][type]
                            *yLeft[N+type]*yRight[N+type];
                } else {
                    y[N+type] += 0.5*param.getCrossBirthRates()[nodeInterval][type][typeOther]
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
    private Node forwardSimulateSubtree(Node subtreeRoot, double startTime, int startType) {

        Node root = new Node();
        setNodeType(root, startType);

        Node currentNode = root;
        int currentType = startType;
        double currentTime = startTime;

        double endTime = param.getNodeTime(subtreeRoot, finalSampleOffset.getArrayValue());

        double[] rates = new double[param.getNTypes()];
        double[] ratesPrime = new double[param.getNTypes()];
        double totalRate, totalRatePrime;

        while (true) {

            // Determine time of next rate shift

            double delta = Utils.globalPrecisionThreshold;
            int interval = param.getIntervalIndex(currentTime + 2*delta);
            double nextRateShiftTime = param.getIntervalEndTimes()[interval];
            double integrationEndTime = Math.min(endTime, nextRateShiftTime);

            // Determine time of next event

            double K = -Math.log(Randomizer.nextDouble());
            double I = 0.0;

            double t = currentTime;
            double dt = (integrationEndTime-currentTime)/ FORWARD_INTEGRATION_STEPS;
            totalRate = getTotalFowardsRate(currentType, currentTime, subtreeRoot, rates);

            int integrationStep;
            for (integrationStep=0; integrationStep< FORWARD_INTEGRATION_STEPS; integrationStep++) {
                double tprime = currentTime + (integrationEndTime-currentTime)*(integrationStep+1)/ FORWARD_INTEGRATION_STEPS;

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

            if (integrationStep == FORWARD_INTEGRATION_STEPS) {
                if (nextRateShiftTime < endTime) {
                    currentTime = nextRateShiftTime;
                    continue;
                } else {
                    break;
                }
            }

            currentNode.setHeight(param.getNodeAge(currentTime, finalSampleOffset.getArrayValue()));

            // Sample event type

            currentType = Randomizer.randomChoicePDF(ratesPrime);

            // Implement event in tree

            Node newNode = new Node();
            setNodeType(newNode, currentType);

            currentNode.addChild(newNode);
            currentNode = newNode;
        }

        currentNode.setHeight(param.getNodeAge(endTime, finalSampleOffset.getArrayValue()));

        switch(getNodeKind(subtreeRoot)) {
            case LEAF:
                currentNode.setNr(subtreeRoot.getNr());
                currentNode.setID(subtreeRoot.getID());
                break;

            case SA:
                Node oldDAChild = subtreeRoot.getDirectAncestorChild();
                Node newDAChild = new Node();

                newDAChild.setHeight(oldDAChild.getHeight());
                setNodeType(newDAChild, currentType);
                newDAChild.setNr(oldDAChild.getNr());
                newDAChild.setID((oldDAChild.getID()));

                currentNode.addChild(newDAChild);
                currentNode.addChild(forwardSimulateSubtree(subtreeRoot.getNonDirectAncestorChild(), endTime, currentType));
                break;

            case INTERNAL:
                int[] childTypes = sampleChildTypes(subtreeRoot, currentType);
                currentNode.addChild(forwardSimulateSubtree(subtreeRoot.getChild(0), endTime, childTypes[0]));
                currentNode.addChild(forwardSimulateSubtree(subtreeRoot.getChild(1), endTime, childTypes[1]));
                break;
        }

        return root;
    }

    /**
     * Retrieve (scaled) result of backward-time integration, ensuring that the query remains
     * inside the bounds of the integration.
     *
     * @param node node below edge on which integration result is to be retrieved.
     * @param time time at which state is to be retrieved
     * @return backward-time integration result at this point on the tree
     */
    private double[] getBackwardsIntegrationResult(Node node, double time) {
        double parentTime = node.isRoot() ? 0.0 : param.getNodeTime(node.getParent(), finalSampleOffset.getArrayValue());
        double adjustedTime = Math.max(time, parentTime + 2*Utils.globalPrecisionThreshold);

        ContinuousOutputModel com = integrationResults[node.getNr()];
        com.setInterpolatedTime(adjustedTime);

        double[] y = com.getInterpolatedState();

        // Trim away small negative values due to numerical integration errors
        for (int i=0; i<y.length; i++) {
            if (y[i] < 0)
                y[i] = 0.0;
        }

        return y;
    }

    private int[] sampleChildTypes(Node node, int parentType) {

        double t = param.getNodeTime(node, finalSampleOffset.getArrayValue());
        int interval = param.getIntervalIndex(t);

        double[] y1 = getBackwardsIntegrationResult(node.getChild(0), t);
        double[] y2 = getBackwardsIntegrationResult(node.getChild(1), t);

        double[][] probs = new double[param.getNTypes()][param.getNTypes()];

        double totalMass = 0.0;
        for (int type1=0; type1<param.getNTypes(); type1++) {
            for (int type2=0; type2<param.getNTypes(); type2++) {

                if (type1 != parentType && type2 != parentType) {
                    probs[type1][type2] = 0.0;
                    continue;
                }

                if (type1 == type2) {
                    probs[type1][type1] = param.getBirthRates()[interval][type1]
                            *y1[param.getNTypes()+type1]*y2[param.getNTypes()+type1];
                } else {
                    int newType = type1 != parentType ? type1 : type2;
                    probs[type1][type2] = param.getCrossBirthRates()[interval][parentType][newType]
                            * 0.5 * y1[param.getNTypes()+type1]*y2[param.getNTypes()+type2];
                }

                totalMass += probs[type1][type2];
            }
        }

        double u = Randomizer.nextDouble()*totalMass;

        for (int type1=0; type1<param.getNTypes(); type1++) {
            for (int type2 = 0; type2 < param.getNTypes(); type2++) {

                if (u < probs[type1][type2]) {
                    return new int[]{type1, type2};
                }

                u -= probs[type1][type2];
            }
        }

        throw new IllegalStateException("Internal node child type sampling loop fell through.");
    }

    /**
     * Compute forward migration rates at a particular time and from a particular
     * type. The results are stored in the provided array, a reference to which
     * is also returned.
     *
     * @param fromType current type
     * @param time time at which to compute rates
     * @param baseNode node at base of edge along which to compute rates.
     * @param result array in which results will be stored.
     * @return reference to array.
     */
    private double[] getForwardsRates(int fromType, double time, Node baseNode, double[] result) {
        double[] y = getBackwardsIntegrationResult(baseNode, time);

        int interval = param.getIntervalIndex(time);

        for (int type=0; type<param.getNTypes(); type++) {
            if (type == fromType) {
                result[type] = 0.0;
                continue;
            }

            result[type] = (param.getCrossBirthRates()[interval][fromType][type] * y[fromType]
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
                result[type] = type == maxRateIdx ? 1e10 : 0.0 ;

        } else {
            // Apply source type prob as rate denominator:

            for (int type=0; type<param.getNTypes(); type++)
                result[type] /= y[param.getNTypes()+fromType];

        }

        return result;
    }

    /**
     * Compute total forward-time transition rate.
     *
     * @param fromType starting type for transition
     * @param time time at which to compute rates
     * @param baseNode base node of edge on which to compute rates
     * @param rates array in which individual rates are stored
     * @return total forward-time transision rate.
     */
    private double getTotalFowardsRate(int fromType, double time, Node baseNode, double[] rates) {
        double totalRate = 0.0;
        getForwardsRates(fromType, time, baseNode, rates);
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

    private long lastRemapSample = -1;

    /**
     * Remap the tree.  Intended to be called by loggers requiring
     * a mapped tree.  Supplying the sample number allows the result
     * of the remapping to be cached and used for other loggers.
     *
     * @param sample sample number at log
     */
    public void remapForLog(long sample) {
        if (!remapOnLogInput.get() || sample == lastRemapSample)
            return;

        doStochasticMapping();
        lastRemapSample = sample;
    }

    @Override
    public void init(PrintStream out) {
        untypedTree.init(out);
    }

    @Override
    public void log(long sample, PrintStream out) {
        remapForLog(sample);

        Tree tree = (Tree) getCurrent();
        out.print("tree STATE_" + sample + " = ");
        // Don't sort, this can confuse CalculationNodes relying on the tree
        //tree.getRoot().sort();
        final int[] dummy = new int[1];
        final String newick = tree.getRoot().toSortedNewick(dummy, true);
        out.print(newick);
        out.print(";");
    }

    @Override
    public void close(PrintStream out) {
        untypedTree.close(out);
    }
}
