package bdmmprime.distribution;

import bdmmprime.parameterization.Parameterization;
import bdmmprime.util.Utils;
import beast.base.core.*;
import beast.base.evolution.speciation.SpeciesTreeDistribution;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.TraitSet;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeInterface;
import beast.base.inference.State;
import beast.base.inference.parameter.RealParameter;
import beast.base.util.HeapSort;
import org.apache.commons.math.special.Gamma;

import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.List;
import java.util.Random;
import java.util.concurrent.*;

/**
 * @author Denise Kuehnert
 * Date: Jul 2, 2013
 * Time: 10:28:16 AM
 */

@Citation(value = "Kuehnert D, Stadler T, Vaughan TG, Drummond AJ. (2016). " +
        "Phylodynamics with migration: \n" +
        "A computational framework to quantify population structure from genomic data. \n" +
        "Mol Biol Evol. 33(8):2102–2116."
        , DOI = "10.1093/molbev/msw064", year = 2016, firstAuthorSurname = "Kuehnert")

@Description("This model implements a multi-deme version of the BirthDeathSkylineModel " +
        "with discrete locations and migration events among demes. " +
        "This implementation also works with sampled ancestor trees.")
public class BirthDeathMigrationDistribution extends SpeciesTreeDistribution {

    public Input<Parameterization> parameterizationInput = new Input<>("parameterization",
            "BDMM parameterization",
            Input.Validate.REQUIRED);

    public Input<Function> finalSampleOffsetInput = new Input<>("finalSampleOffset",
            "If provided, the difference in time between the final sample and the end of the BD process.",
            new RealParameter("0.0"));

    public Input<RealParameter> frequenciesInput = new Input<>("frequencies",
            "The equilibrium frequencies for each type",
            new RealParameter("1.0"));

    public Input<TraitSet> typeTraitSetInput = new Input<>("typeTraitSet",
            "Trait set specifying sample trait values.");

    public Input<String> typeLabelInput = new Input<>("typeLabel",
            "Attribute key used to specify sample trait values in tree.");

    public Input<Boolean> conditionOnSurvivalInput = new Input<>("conditionOnSurvival",
            "Condition on at least one surviving lineage. (Default true.)",
            true);

    public Input<Boolean> conditionOnRootInput = new Input<>("conditionOnRoot",
            "Condition on root age, not time of origin.", false);

    public Input<Boolean> useAnalyticalSingleTypeSolutionInput = new Input<>("useAnalyticalSingleTypeSolution",
            "Use the analytical SABDSKY tree prior when the model has only one type.",
            true);

    public Input<Double> relativeToleranceInput = new Input<>("relTolerance",
            "Relative tolerance for numerical integration.",
            1e-7);

    public Input<Double> absoluteToleranceInput = new Input<>("absTolerance",
            "Absolute tolerance for numerical integration.",
            1e-100 /*Double.MIN_VALUE*/);

    public Input<Boolean> parallelizeInput = new Input<>(
            "parallelize",
            "Whether or not to parallelized the calculation of subtree likelihoods. " +
                    "(Default true.)",
            true);

    /* If a large number a cores is available (more than 8 or 10) the
    calculation speed can be increased by diminishing the parallelization
    factor. On the contrary, if only 2-4 cores are available, a slightly
    higher value (1/5 to 1/8) can be beneficial to the calculation speed. */
    public Input<Double> minimalProportionForParallelizationInput = new Input<>(
            "parallelizationFactor",
            "the minimal relative size the two children " +
                    "subtrees of a node must have to start parallel " +
                    "calculations on the children. (default: 1/10). ",
            1.0 / 10);

    public Input<Boolean> storeNodeTypesInput = new Input<>("storeNodeTypes",
            "store tip node types? this assumes that tip types cannot " +
                    "change (default false)", false);

    public Input<String> savePartialLikelihoodsToFileInput = new Input<>("savePartialLikelihoodsToFile",
            "If provided, the name of a file to which a tree annotated with partial likelihoods " +
                    "will be written.  This is useful for debugging the cause of chain initialization failure.");

    private int[] nodeStates;

    private final boolean debug = false;
//    private final boolean debug = true;

    private double[] startTypeProbs, storedStartTypeProbs;
    private boolean[] isRhoTip;

    private Parameterization parameterization;
    private Function finalSampleOffset;

    private double[][] pInitialConditions;

    private static boolean isParallelizedCalculation;
    private static double minimalProportionForParallelization;
    private double parallelizationThreshold;
    private static ThreadPoolExecutor pool;

    private double[] weightOfNodeSubTree;

    private TreeInterface tree;

    private int originalLeafCount;

    @Override
    public void initAndValidate() {
        parameterization = parameterizationInput.get();
        tree = treeInput.get();

        finalSampleOffset = finalSampleOffsetInput.get();

        if (parameterization.getNTypes() != 1 && (typeTraitSetInput.get() == null && typeLabelInput.get() == null))
            throw new RuntimeException("Error: For models with >1 type, either typeTraitSet or typeLabel must be specified.");

        if (frequenciesInput.get().getDimension() != parameterization.getNTypes())
            throw new RuntimeException("Error: dimension of equilibrium frequencies " +
                    "parameter must match number of types.");

        double freqSum = 0;
        for (double f : frequenciesInput.get().getValues()) freqSum += f;
        if (Math.abs(1.0 - freqSum) > 1e-10)
            throw new RuntimeException("Error: equilibrium frequencies must add " +
                    "up to 1 but currently add to " + freqSum + ".");

        int nLeaves = tree.getLeafNodeCount();

        weightOfNodeSubTree = new double[nLeaves * 2];

        isParallelizedCalculation = parallelizeInput.get();
        minimalProportionForParallelization = minimalProportionForParallelizationInput.get();

        if (isParallelizedCalculation) executorBootUp();

        if (storeNodeTypesInput.get()) {

            nodeStates = new int[nLeaves];

            for (Node node : tree.getExternalNodes()) {
                nodeStates[node.getNr()] = getNodeType(node, true);
            }
        }

        startTypeProbs = new double[parameterization.getNTypes()];
        storedStartTypeProbs = new double[parameterization.getNTypes()];

        // Determine which, if any, of the leaf ages correspond exactly to
        // rho sampling times.
        isRhoTip = new boolean[tree.getLeafNodeCount()];
        for (int nodeNr = 0; nodeNr < tree.getLeafNodeCount(); nodeNr++) {
            isRhoTip[nodeNr] = false;
            double nodeTime = parameterization.getNodeTime(tree.getNode(nodeNr), finalSampleOffset.getArrayValue());
            for (double rhoSampTime : parameterization.getRhoSamplingTimes()) {
                if (Utils.equalWithPrecision(rhoSampTime, nodeTime)) {
                    isRhoTip[nodeNr] = true;
                    break;
                }
            }
        }

        originalLeafCount = tree.getLeafNodeCount();
    }

    @Override
    public double calculateTreeLogLikelihood(TreeInterface dummyTree) {

        if (tree.getLeafNodeCount() != originalLeafCount)
            initAndValidate();

        Node root = tree.getRoot();

        if (Utils.lessThanWithPrecision(parameterization.getNodeTime(tree.getRoot(), finalSampleOffset.getArrayValue()), 0)) {
            if (savePartialLikelihoodsToFileInput.get() != null)
                Log.err("Tree MRCA older than start of process.");
            logP = Double.NEGATIVE_INFINITY;
            return logP;
        }

        if (conditionOnRootInput.get() && tree.getRoot().isFake()) {
            if (savePartialLikelihoodsToFileInput.get() != null)
                Log.err("Tree root is a sampled ancestor but we're conditioning on a root birth.");
            logP = Double.NEGATIVE_INFINITY;
            return logP;
        }

        // Use the exact solution in the case of a single type
        if (useAnalyticalSingleTypeSolutionInput.get() && parameterization.getNTypes() == 1
                && savePartialLikelihoodsToFileInput.get() == null) {
            logP = getSingleTypeTreeLogLikelihood();
            return logP;
        }

        P0GeSystem system = new P0GeSystem(parameterization,
                absoluteToleranceInput.get(), relativeToleranceInput.get());

        updateParallelizationThreshold();

        updateInitialConditionsForP();

        double conditionDensity = 0.0;
        double[] extinctionProb = pInitialConditions[pInitialConditions.length - 1];
        if (conditionOnRootInput.get()) {

            int intervalIndex = parameterization.getIntervalIndex(0);
            for (int type1=0; type1<parameterization.getNTypes(); type1++) {
                for (int type2=0; type2<parameterization.getNTypes(); type2++) {
                    double rate = type1 == type2
                            ? parameterization.getBirthRates()[intervalIndex][type1]
                            : parameterization.getCrossBirthRates()[intervalIndex][type1][type2];

                    if (rate == 0.0)
                        continue;

                   conditionDensity += rate*frequenciesInput.get().getArrayValue(type1)
                           * (1-extinctionProb[type1])
                           * (1-extinctionProb[type2]);
                }
            }
        } else if (conditionOnSurvivalInput.get()) {

            for (int type = 0; type < parameterization.getNTypes(); type++)
                conditionDensity += frequenciesInput.get().getArrayValue(type)
                        * (1-extinctionProb[type]);

        } else {
            conditionDensity = 1.0;
        }

        if (conditionDensity < 0)
            return Double.NEGATIVE_INFINITY;

        P0GeState finalP0Ge;
        if (conditionOnRootInput.get()) {

            // Condition on a known root time:

            finalP0Ge = new P0GeState(parameterization.getNTypes());

            Node child1 = root.getChild(0);
            Node child2 = root.getChild(1);

            P0GeState child1state = calculateSubtreeLikelihood(child1, 0,
                    parameterization.getNodeTime(child1, finalSampleOffset.getArrayValue()),
                    system, 0);
            P0GeState child2state = calculateSubtreeLikelihood(child2, 0,
                    parameterization.getNodeTime(child2, finalSampleOffset.getArrayValue()),
                    system, 0);

            int intervalIndex = parameterization.getIntervalIndex(0);
            for (int type1=0; type1<parameterization.getNTypes(); type1++) {
                finalP0Ge.p0[type1] = child1state.p0[type1];
                for (int type2=0; type2<parameterization.getNTypes(); type2++) {
                    double rate = type2 == type1
                            ? parameterization.getBirthRates()[intervalIndex][type1]
                            : parameterization.getCrossBirthRates()[intervalIndex][type1][type2];

                    if (rate == 0.0)
                        continue;

                    finalP0Ge.ge[type1] = finalP0Ge.ge[type1].addTo(
                            child1state.ge[type1]
                            .multiplyBy(child2state.ge[type2])
                                    .addTo(child1state.ge[type2].multiplyBy(child2state.ge[type1]))
                                    .scalarMultiplyBy(0.5*rate));
                }
            }

        } else {

            // Condition on origin time:

            finalP0Ge = calculateSubtreeLikelihood(root, 0,
                    parameterization.getNodeTime(tree.getRoot(), finalSampleOffset.getArrayValue()),
                    system, 0);
        }

        if (debug) System.out.print("Final state: " + finalP0Ge);

        SmallNumber PrSN = new SmallNumber(0);
        for (int startType = 0; startType < parameterization.getNTypes(); startType++) {

            SmallNumber jointProb = finalP0Ge
                    .ge[startType]
                    .scalarMultiplyBy(frequenciesInput.get().getArrayValue(startType));

            if (jointProb.getMantissa() > 0) {
                startTypeProbs[startType] = jointProb.log();
                PrSN = PrSN.addTo(jointProb);
            } else {
                startTypeProbs[startType] = Double.NEGATIVE_INFINITY;
            }
        }

        // Normalize start type probs:
        for (int startType = 0; startType < parameterization.getNTypes(); startType++) {
            startTypeProbs[startType] -= PrSN.log();
            startTypeProbs[startType] = Math.exp(startTypeProbs[startType]);
        }

        PrSN = PrSN.scalarMultiplyBy(1 / conditionDensity);

        logP = PrSN.log();

        // Convert from oriented to labeled tree probability density:
        int internalNodeCount = tree.getLeafNodeCount() - ((Tree) tree).getDirectAncestorNodeCount() - 1;
        logP += Math.log(2) * internalNodeCount - Gamma.logGamma(tree.getLeafNodeCount()+1);

        if (debug) System.out.println("\nlogP = " + logP);

        if (savePartialLikelihoodsToFileInput.get() != null) {
            Log.info("Writing partial BDMM-Prime likelihoods to file '"
                    + savePartialLikelihoodsToFileInput.get() + "'");
            try (PrintStream ps = new PrintStream(savePartialLikelihoodsToFileInput.get())) {
                ps.println(getNewickWithP0GeMetadata(tree.getRoot()));
            } catch (FileNotFoundException e) {
                Log.err("Failed to open output file '" + savePartialLikelihoodsToFileInput.get() + "'");
            }
        }

        return logP;
    }

    private int getNodeType(Node node, Boolean init) {


        if (storeNodeTypesInput.get() && !init)
            return nodeStates[node.getNr()];

        int nodeType;

        if (parameterization.getNTypes()>1) {
            String nodeTypeName;

            if (typeTraitSetInput.get() != null)
                nodeTypeName = typeTraitSetInput.get().getStringValue(node.getID());
            else {
                Object metaData = node.getMetaData(typeLabelInput.get());
                if (metaData instanceof Double)
                    nodeTypeName = String.valueOf(Math.round((double)metaData));
                else
                    nodeTypeName = metaData.toString();
            }

            nodeType = parameterization.getTypeSet().getTypeIndex(nodeTypeName);

        } else {
            nodeType = 0;
        }

        if (storeNodeTypesInput.get())
            nodeStates[node.getNr()] = nodeType;

        return nodeType;
    }


    /**
     * @param node    Node below edge.
     * @param tTop    Time of start (top) of edge.
     * @param tBottom Time of end (bottom) of edge.
     * @param system  Object describing ODEs to integrate.
     * @return State at top of edge.
     */
    private P0GeState calculateSubtreeLikelihood(Node node, double tTop, double tBottom,
                                                 P0GeSystem system, int depth) {

        if (debug) {
            debugMessage("*** Evaluating subtree for node " + node +
                            " and edge between times " + tTop + " and " + tBottom + " ...",
                    depth);
        }

        P0GeState state = new P0GeState(parameterization.getNTypes());

        int intervalIdx = parameterization.getIntervalIndex(tBottom);

        if (node.isLeaf()) { // sampling event

            // Incorporate pre-evaluated p0 values into state
            System.arraycopy(pInitialConditions[node.getNr()], 0, state.p0, 0, system.nTypes);

            int nodeType = getNodeType(node, false);

            if (nodeType == -1) { //unknown state

                //TODO test if SA model case is properly implemented (not tested!)
                for (int type = 0; type < parameterization.getNTypes(); type++) {

                    if (isRhoTip[node.getNr()]) {
                        state.ge[type] = new SmallNumber(
                                (system.r[intervalIdx][type] + state.p0[type]
                                        * (1 - system.r[intervalIdx][type]))
                                        * system.rho[intervalIdx][type]);
                    } else {
                        state.ge[type] = new SmallNumber(
                                (system.r[intervalIdx][type] + state.p0[type] * (1 - system.r[intervalIdx][type]))
                                        * system.s[intervalIdx][type]);
                        // with SA: ψ_i(r + (1 − r)p_i(τ))
                    }
                }
            } else {

                if (isRhoTip[node.getNr()]) {

                    state.ge[nodeType] = new SmallNumber(
                            (system.r[intervalIdx][nodeType] + state.p0[nodeType]
                                    * (1 - system.r[intervalIdx][nodeType]))
                                    * system.rho[intervalIdx][nodeType]);
                } else {
                    state.ge[nodeType] = new SmallNumber(
                            (system.r[intervalIdx][nodeType] + state.p0[nodeType]
                                    * (1 - system.r[intervalIdx][nodeType]))
                                    * system.s[intervalIdx][nodeType]);
                    // with SA: ψ_i(r + (1 − r)p_i(τ))
                }

            }

            // Incorporate rho sampling if we're on a boundary:
            if (isRhoTip[node.getNr()]) {
                for (int type = 0; type < parameterization.getNTypes(); type++) {
                    state.p0[type] *= (1 - system.rho[intervalIdx][type]);
                }
            }

            if (debug) debugMessage("Sampling at time " + tBottom, depth);

        } else if (node.getChildCount() == 2) {  // birth / infection event or sampled ancestor

            if (node.getChild(0).isDirectAncestor() || node.getChild(1).isDirectAncestor()) {   // found a sampled ancestor

                int childIndex = 0;

                if (node.getChild(childIndex).isDirectAncestor())
                    childIndex = 1;

                P0GeState g = calculateSubtreeLikelihood(
                        node.getChild(childIndex), tBottom,
                        parameterization.getNodeTime(node.getChild(childIndex), finalSampleOffset.getArrayValue()),
                        system, depth + 1);

                int saNodeType = getNodeType(node.getChild(childIndex ^ 1), false); // get state of direct ancestor, XOR operation gives 1 if childIndex is 0 and vice versa

                //TODO test if properly implemented (not tested!)
                if (saNodeType == -1) { // unknown state
                    for (int type = 0; type < parameterization.getNTypes(); type++) {
                        if (!isRhoTip[node.getChild(childIndex ^ 1).getNr()]) {

                            state.p0[type] = g.p0[type];
                            state.ge[type] = g.ge[type].scalarMultiplyBy(system.s[intervalIdx][type]
                                    * (1 - system.r[intervalIdx][type]));

                        } else {
                            // TODO COME BACK AND CHANGE (can be dealt with with getAllPInitialConds)
                            state.p0[type] = g.p0[type] * (1 - system.rho[intervalIdx][type]);
                            state.ge[type] = g.ge[type].scalarMultiplyBy(system.rho[intervalIdx][type]
                                    * (1 - system.r[intervalIdx][type]));

                        }
                    }
                } else {
                    if (!isRhoTip[node.getChild(childIndex ^ 1).getNr()]) {

                        state.p0[saNodeType] = g.p0[saNodeType];
                        state.ge[saNodeType] = g.ge[saNodeType]
                                .scalarMultiplyBy(system.s[intervalIdx][saNodeType]
                                        * (1 - system.r[intervalIdx][saNodeType]));

                    } else {
                        // TODO COME BACK AND CHANGE (can be dealt with with getAllPInitialConds)
                        state.p0[saNodeType] = g.p0[saNodeType]
                                * (1 - system.rho[intervalIdx][saNodeType]);
                        state.ge[saNodeType] = g.ge[saNodeType]
                                .scalarMultiplyBy(system.rho[intervalIdx][saNodeType]
                                        * (1 - system.r[intervalIdx][saNodeType]));

                    }
                }
            } else {   // birth / infection event

                int indexFirstChild = 0;
                if (node.getChild(1).getNr() > node.getChild(0).getNr())
                    indexFirstChild = 1; // always start with the same child to avoid numerical differences

                int indexSecondChild = Math.abs(indexFirstChild - 1);

                P0GeState childState1 = null, childState2 = null;

                // evaluate if the next step in the traversal should be split between one new thread and the currrent thread and run in parallel.

                if (isParallelizedCalculation
                        && weightOfNodeSubTree[node.getChild(indexFirstChild).getNr()] > parallelizationThreshold
                        && weightOfNodeSubTree[node.getChild(indexSecondChild).getNr()] > parallelizationThreshold) {

                    try {
                        // start a new thread to take care of the second subtree
                        Future<P0GeState> secondChildTraversal = pool.submit(
                                new TraversalService(node.getChild(indexSecondChild), tBottom,
                                        parameterization.getNodeTime(node.getChild(indexSecondChild), finalSampleOffset.getArrayValue()),
                                        depth + 1));

                        childState1 = calculateSubtreeLikelihood(
                                node.getChild(indexFirstChild), tBottom,
                                parameterization.getNodeTime(node.getChild(indexFirstChild), finalSampleOffset.getArrayValue()),
                                system, depth + 1);
                        childState2 = secondChildTraversal.get();
                    } catch (InterruptedException | ExecutionException e) {
                        e.printStackTrace();

                        System.exit(1);
                    }

                } else {
                    childState1 = calculateSubtreeLikelihood(node.getChild(
                            indexFirstChild), tBottom,
                            parameterization.getNodeTime(node.getChild(indexFirstChild), finalSampleOffset.getArrayValue()),
                            system, depth + 1);
                    childState2 = calculateSubtreeLikelihood(node.getChild(indexSecondChild), tBottom,
                            parameterization.getNodeTime(node.getChild(indexSecondChild), finalSampleOffset.getArrayValue()),
                            system, depth + 1);
                }


                if (debug) debugMessage("Infection at time " + tBottom, depth);

                for (int childType = 0; childType < parameterization.getNTypes(); childType++) {

                    state.p0[childType] = childState1.p0[childType];
                    state.ge[childType] = childState1.ge[childType]
                            .multiplyBy(childState2.ge[childType])
                            .scalarMultiplyBy(system.b[intervalIdx][childType]);

                    for (int otherChildType = 0; otherChildType < parameterization.getNTypes(); otherChildType++) {
                        if (otherChildType == childType)
                            continue;

                        state.ge[childType] = state.ge[childType]
                                .addTo((childState1.ge[childType].multiplyBy(childState2.ge[otherChildType]))
                                        .addTo(childState1.ge[otherChildType].multiplyBy(childState2.ge[childType]))
                                        .scalarMultiplyBy(0.5 * system.b_ij[intervalIdx][childType][otherChildType]));
                    }


                    if (Double.isInfinite(state.p0[childType])) {
                        throw new RuntimeException("infinite likelihood");
                    }
                }
            }

        }

        if (debug) debugMessage("State at base of edge: " + state, depth);
        if (savePartialLikelihoodsToFileInput.get() != null) {
            node.setMetaData("p0", getP0MetadataString(state));
            node.setMetaData("ge", getGeMetadataString(state));
            node.setMetaData("geZero", partialLikelihodZero(state) ? "true" : "false");
            node.setMetaData("interval", String.valueOf(intervalIdx));
        }

        integrateP0Ge(node, tTop, state, system);

        if (debug)
            debugMessage("State at top of edge: " + state + "\n", depth);

        return state;
    }

    /**
     * Print message to stdout with given indentation depth.
     *
     * @param message debug message
     * @param depth   indentation level
     */
    private void debugMessage(String message, int depth) {
        for (int i = 0; i < depth; i++)
            System.out.print("  ");

        System.out.println(message);
    }

    private String getGeMetadataString(P0GeState state) {
        StringBuilder sb = new StringBuilder();

        sb.append("{");
        for (int i=0; i<state.dimension; i++) {
            if (i > 0)
                sb.append(",");
            sb.append(state.ge[i].toString());
        }
        sb.append("}");

        return sb.toString();
    }

    private String getP0MetadataString(P0GeState state) {
        StringBuilder sb = new StringBuilder();

        sb.append("{");
        for (int i=0; i<state.dimension; i++) {
            if (i > 0)
                sb.append(",");
            sb.append(state.p0[i]);
        }
        sb.append("}");

        return sb.toString();
    }

    private boolean partialLikelihodZero(P0GeState state) {

        for (int i=0; i<state.dimension; i++) {
            if (state.ge[i].getMantissa() != 0.0)
                return false;
        }

        return true;
    }

    private String getNewickWithP0GeMetadata(Node node) {
        StringBuilder sb = new StringBuilder();

        if (!node.isLeaf()) {
            sb.append("(");
            boolean isFirst = true;
            for (Node child : node.getChildren()) {
                if (isFirst)
                    isFirst = false;
                else
                    sb.append(",");
                sb.append(getNewickWithP0GeMetadata(child));
            }
            sb.append(")");
        }

        if (node.getID() != null)
            sb.append(node.getID());

        String p0Metadata = (String) node.getMetaData("p0");
        String geMetadata = (String) node.getMetaData("ge");
        String geZero = (String) node.getMetaData("geZero");

        if (p0Metadata != null && geMetadata != null && geZero != null)
            sb.append("[&p0=").append(p0Metadata)
                    .append(",ge=").append(geMetadata)
                    .append(",geZero=").append(node.getMetaData("geZero"))
                    .append(",interval=").append(node.getMetaData("interval"))
                    .append("]");

        sb.append(":");
        if (node.isRoot())
            sb.append("0.0;");
        else
            sb.append(node.getParent().getHeight()-node.getHeight());

        return sb.toString();
    }

    /**
     * @return retrieve current set of root type probabilities.
     */
    double[] getStartTypeProbs() {
        return startTypeProbs;
    }

    /**
     * Compute all initial conditions for all future integrations on p0 equations.
     */
    private void updateInitialConditionsForP() {

        P0System p0System = new P0System(parameterization,
                absoluteToleranceInput.get(), relativeToleranceInput.get());

        int leafCount = tree.getLeafNodeCount();
        double[] leafTimes = new double[leafCount];
        int[] indicesSortedByLeafTime = new int[leafCount];

        for (int i = 0; i < leafCount; i++) { // get all leaf times
            leafTimes[i] = parameterization.getNodeTime(tree.getNode(i), finalSampleOffset.getArrayValue());
//            leafTimes[i] = parameterization.getTotalProcessLength() - tree.getNode(i).getHeight();
            indicesSortedByLeafTime[i] = i;
        }

        HeapSort.sort(leafTimes, indicesSortedByLeafTime);
        //"sort" sorts in ascending order, so we have to be careful since the
        // integration starts from the leaves at time T and goes up to the
        // root at time 0 (or >0)


        // The initial value is zero, so that all modifications can be expressed
        // as products.
        P0State p0State = new P0State(p0System.nTypes);
        for (int type=0; type<p0System.nTypes; type++)
            p0State.p0[type] = 1.0;

        pInitialConditions = new double[leafCount + 1][p0System.nTypes];

        double tprev = p0System.totalProcessLength;

        for (int i = leafCount - 1; i >= 0; i--) {
            double t = leafTimes[indicesSortedByLeafTime[i]];

            //If the next higher leaf is actually at the same height, store previous results and skip iteration
            if (Utils.equalWithPrecision(t, tprev)) {
                tprev = t;
                if (i < leafCount-1) {
                    pInitialConditions[indicesSortedByLeafTime[i]] =
                            pInitialConditions[indicesSortedByLeafTime[i + 1]];
                } else {
                    System.arraycopy(p0State.p0, 0,
                            pInitialConditions[indicesSortedByLeafTime[i]], 0,
                            p0System.nTypes);
                }
                continue;
            }

            // Only include rho contribution when starting integral to earlier times.
            // This means that the value of pInitialConditions will always require the
            // inclusion of a (1-rho) factor if it lies on an interval boundary, just
            // as for the Ge calculation.
            int prevIndex = parameterization.getIntervalIndex(tprev);
            if (Utils.equalWithPrecision(parameterization.getIntervalEndTimes()[prevIndex], tprev)) {
                for (int type = 0; type < parameterization.getNTypes(); type++) {
                    p0State.p0[type] *= (1 - parameterization.getRhoValues()[prevIndex][type]);
                }
            }

            integrateP0(tprev, t, p0State, p0System);

            System.arraycopy(p0State.p0, 0,
                    pInitialConditions[indicesSortedByLeafTime[i]], 0, p0System.nTypes);
            tprev = t;
        }

        if (Utils.greaterThanWithPrecision(tprev , 0.0)) {
            int prevIndex = parameterization.getIntervalIndex(tprev);
            if (Utils.equalWithPrecision(parameterization.getIntervalEndTimes()[prevIndex], tprev)) {
                for (int type = 0; type < parameterization.getNTypes(); type++) {
                    p0State.p0[type] *= (1 - parameterization.getRhoValues()[prevIndex][type]);
                }
            }
        }

        integrateP0(tprev, 0, p0State, p0System);
        System.arraycopy(p0State.p0, 0,
                pInitialConditions[leafCount], 0, p0System.nTypes);
    }

    /**
     * Perform integration on differential equations p
     */
    private void integrateP0(double tStart, double tEnd, P0State state, P0System system) {

        double thisTime = tStart;
        int thisInterval = parameterization.getIntervalIndex(thisTime);
        int endInterval = parameterization.getIntervalIndex(tEnd);

        while (thisInterval > endInterval) {

            double nextTime = system.intervalEndTimes[thisInterval-1];

            if (Utils.lessThanWithPrecision(nextTime , thisTime)) {
                system.setInterval(thisInterval);
                system.integrate(state, thisTime, nextTime);
            }

            if (Utils.greaterThanWithPrecision(nextTime, tEnd)) {
                for (int i = 0; i < system.nTypes; i++)
                    state.p0[i] *= (1 - system.rho[thisInterval - 1][i]);
            }

            thisTime = nextTime;
            thisInterval -= 1;

        }

        if (Utils.greaterThanWithPrecision(thisTime, tEnd)) {
            system.setInterval(thisInterval);
            system.integrate(state, thisTime, tEnd);
        }
    }

    /**
     * Integrate state along edge above baseNode until time tTop according to system.
     *
     * @param baseNode node at base of edge
     * @param tTop     time at top of edge
     * @param state    ODE variables at bottom of edge
     * @param system   ODE system to integrate
     */
    private void integrateP0Ge(Node baseNode, double tTop, P0GeState state, P0GeSystem system) {

        // pgScaled contains the set of initial conditions scaled made to fit
        // the requirements on the values 'double' can represent. It also
        // contains the factor by which the numbers were multiplied.
        ScaledNumbers pgScaled = state.getScaledState();

        double thisTime = parameterization.getNodeTime(baseNode, finalSampleOffset.getArrayValue());
        int thisInterval = parameterization.getIntervalIndex(thisTime);
        int endInterval = parameterization.getIntervalIndex(tTop);
        double oneMinusRho;

        system.setInterval(thisInterval);

        while (thisInterval > endInterval) {
            double nextTime = system.intervalEndTimes[thisInterval-1];

            if (Utils.lessThanWithPrecision(nextTime, thisTime)) {
                pgScaled = system.safeIntegrate(pgScaled, thisTime, nextTime);

                state.setFromScaledState(pgScaled.getEquation(), pgScaled.getScalingFactor());

                if (Utils.greaterThanWithPrecision(nextTime, tTop)) {
                    for (int i = 0; i < parameterization.getNTypes(); i++) {
                        oneMinusRho = 1 - system.rho[thisInterval - 1][i];
                        state.p0[i] *= oneMinusRho;
                        state.ge[i] = state.ge[i].scalarMultiplyBy(oneMinusRho);
                    }
                }

                // 'rescale' the results of the last integration to prepare for the next integration step
                pgScaled = state.getScaledState();
            }

            thisTime = nextTime;
            thisInterval -= 1;

            system.setInterval(thisInterval);
        }

         // solve PG , store solution temporarily integrationResults
        if (Utils.greaterThanWithPrecision(thisTime, tTop))
            pgScaled = system.safeIntegrate(pgScaled, thisTime, tTop);

        // 'unscale' values in integrationResults so as to retrieve accurate values after the integration.
        state.setFromScaledState(pgScaled.getEquation(), pgScaled.getScalingFactor());
    }


    /**
     * Perform an initial traversal of the tree to get the 'weights' (sum of all its edges lengths) of all sub-trees
     * Useful for performing parallelized calculations on the tree.
     * The weights of the subtrees tell us the depth at which parallelization should stop, so as to not parallelize on subtrees that are too small.
     * Results are stored in 'weightOfNodeSubTree' array
     *
     * @param tree tree whose subtree to compute weights of
     */
    private void getAllSubTreesWeights(TreeInterface tree) {
        Node root = tree.getRoot();
        double weight = 0;
        for (final Node child : root.getChildren()) {
            weight += getSubTreeWeight(child);
        }
        weightOfNodeSubTree[root.getNr()] = weight;
    }

    /**
     * Perform an initial traversal of the subtree to get its 'weight': sum of all its edges.
     *
     * @param node root of subtree
     * @return weight
     */
    private double getSubTreeWeight(Node node) {

        // if leaf, stop recursion, get length of branch above and return
        if (node.isLeaf()) {
            weightOfNodeSubTree[node.getNr()] = node.getLength();
            return node.getLength();
        }

        // else, iterate over the children of the node
        double weight = 0;
        for (final Node child : node.getChildren()) {
            weight += getSubTreeWeight(child);
        }
        // add length of parental branch
        weight += node.getLength();
        // store the value
        weightOfNodeSubTree[node.getNr()] = weight;

        return weight;
    }

    private void updateParallelizationThreshold() {
        if (isParallelizedCalculation) {
            getAllSubTreesWeights(tree);
            // set 'parallelizationThreshold' to a fraction of the whole tree weight.
            // The size of this fraction is determined by a tuning parameter. This parameter should be adjusted (increased) if more computation cores are available
            parallelizationThreshold = weightOfNodeSubTree[tree.getRoot().getNr()] * minimalProportionForParallelization;
        }
    }

    private static void executorBootUp() {
        ExecutorService executor = Executors.newCachedThreadPool();
        pool = (ThreadPoolExecutor) executor;
    }

    class TraversalService implements Callable<P0GeState> {

        int depth;
        Node rootSubtree;
        protected double from;
        protected double to;
        P0GeSystem PG;

        TraversalService(Node root, double from, double to, int depth) {
            this.rootSubtree = root;
            this.from = from;
            this.to = to;
            this.depth = depth;

            PG = new P0GeSystem(parameterization,
                    absoluteToleranceInput.get(),
                    relativeToleranceInput.get());
        }

        @Override
        public P0GeState call() {
            // traverse the tree in a potentially-parallelized way
            return calculateSubtreeLikelihood(rootSubtree, from, to, PG, depth);
        }
    }

    /* --- Exact calculation for single type case --- */

    private double get_p_i(double lambda, double mu, double psi, double A, double B, double t_i, double t) {

        if (lambda > 0.0) {
            double v = Math.exp(A * (t_i - t)) * (1 + B);
            return (lambda + mu + psi - A * (v - (1 - B)) / (v + (1 - B)))
                    / (2 * lambda);
        } else {
            // The limit of p_i as lambda -> 0
            return 0.5;
        }
    }

    private double get_q_i(double A, double B, double t_i, double t) {
        double v = Math.exp(A * (t_i - t));
        return 4 * v / Math.pow(v*(1+B) + (1-B), 2.0);
    }

    private void computeConstants(double[] A, double[] B) {

        for (int i=parameterization.getTotalIntervalCount()-1; i>=0; i--) {

            double p_i_prev;
            if (i + 1 < parameterization.getTotalIntervalCount()) {
                p_i_prev = get_p_i(parameterization.getBirthRates()[i+1][0],
                        parameterization.getDeathRates()[i+1][0],
                        parameterization.getSamplingRates()[i+1][0],
                        A[i+1], B[i+1],
                        parameterization.getIntervalEndTimes()[i+1],
                        parameterization.getIntervalEndTimes()[i]);
            } else {
                p_i_prev = 1.0;
            }

            double rho_i = parameterization.getRhoValues()[i][0];
            double lambda_i = parameterization.getBirthRates()[i][0];
            double mu_i = parameterization.getDeathRates()[i][0];
            double psi_i = parameterization.getSamplingRates()[i][0];

            A[i] = Math.sqrt((lambda_i-mu_i-psi_i)*(lambda_i-mu_i-psi_i) + 4*lambda_i*psi_i);
            B[i] = ((1 - 2*(1-rho_i)*p_i_prev)*lambda_i + mu_i + psi_i)/A[i];
        }
    }

    private double getSingleTypeTreeLogLikelihood() {

        double[] A = new double[parameterization.getTotalIntervalCount()];
        double[] B = new double[parameterization.getTotalIntervalCount()];

        computeConstants(A, B);

        double logP;

        if (conditionOnRootInput.get()) {

            double t_root = parameterization.getNodeTime(tree.getRoot(), finalSampleOffset.getArrayValue());
            logP = getSingleTypeSubtreeLogLikelihood(tree.getRoot().getChild(0), t_root, A, B)
                    + getSingleTypeSubtreeLogLikelihood(tree.getRoot().getChild(1), t_root, A, B)
                    + Math.log(2);

            int i_root = parameterization.getIntervalIndex(t_root);
            logP += 2*Math.log(get_q_i(A[i_root], B[i_root], parameterization.getIntervalEndTimes()[i_root], 0.0));

        } else {

            logP = getSingleTypeSubtreeLogLikelihood(tree.getRoot(), 0.0, A, B);

            int i = parameterization.getIntervalIndex(0.0);
            logP += Math.log(get_q_i(A[i], B[i], parameterization.getIntervalEndTimes()[i], 0.0));

        }

        if (conditionOnSurvivalInput.get() || conditionOnRootInput.get()) {
            int i = parameterization.getIntervalIndex(0.0);
            double p_i = get_p_i(parameterization.getBirthRates()[i][0],
                    parameterization.getDeathRates()[i][0],
                    parameterization.getSamplingRates()[i][0],
                    A[i], B[i],
                    parameterization.getIntervalEndTimes()[i], 0.0);

            if (p_i == 1)
                return Double.NEGATIVE_INFINITY; // Following BDSKY's behaviour

//            logP -= Math.log(1.0 - p_i);
            if (conditionOnRootInput.get())
                logP -= 2.0 * Math.log(1.0 - p_i);
            else
                logP -= Math.log(1.0 - p_i);
        }

        // Account for possible label permutations
        logP += -Gamma.logGamma(tree.getLeafNodeCount() + 1);

        return logP;
    }

    private double getSingleTypeSubtreeLogLikelihood(Node subtreeRoot, double timeOfSubtreeRootEdgeTop,
                                                     double[] A, double[] B) {

        double t_node = parameterization.getNodeTime(subtreeRoot, finalSampleOffset.getArrayValue());
        int i = parameterization.getIntervalIndex(t_node);
        double t_i = parameterization.getIntervalEndTimes()[i];

        double rho_i = parameterization.getRhoValues()[i][0];
        double lambda_i = parameterization.getBirthRates()[i][0];
        double mu_i = parameterization.getDeathRates()[i][0];
        double psi_i = parameterization.getSamplingRates()[i][0];
        double r_i = parameterization.getRemovalProbs()[i][0];
        double r_iplus1 = i+1 < parameterization.getTotalIntervalCount()
                ? parameterization.getRemovalProbs()[i+1][0]
                : 1.0;

        double logP;

        if (subtreeRoot.isLeaf()) {
            // Leaf Node

            logP = 0.0;

            if (isRhoTip[subtreeRoot.getNr()]) {

                double p_iplus1 = i + 1 < parameterization.getTotalIntervalCount()
                        ? get_p_i(parameterization.getBirthRates()[i + 1][0],
                        parameterization.getDeathRates()[i + 1][0],
                        parameterization.getSamplingRates()[i + 1][0],
                        A[i + 1], B[i + 1], parameterization.getIntervalEndTimes()[i+1], t_node)
                        : 1.0;

                logP += Math.log(rho_i*(r_iplus1 + (1 - r_iplus1) * p_iplus1));

            } else {

                logP += Math.log(psi_i)
                        + Math.log(r_i + (1-r_i)*get_p_i(lambda_i, mu_i, psi_i, A[i], B[i], t_i, t_node))
                        - Math.log(get_q_i(A[i], B[i], t_i, t_node));

            }


        } else if (subtreeRoot.isFake()) {
            // SA node

            logP = getSingleTypeSubtreeLogLikelihood(subtreeRoot.getNonDirectAncestorChild(), t_node, A, B);

            if (isRhoTip[subtreeRoot.getDirectAncestorChild().getNr()]) {

                double q_iplus1 = i + 1 < parameterization.getTotalIntervalCount()
                        ? get_q_i(A[i+1], B[i+1], parameterization.getIntervalEndTimes()[i+1], t_node)
                        : 1.0;

                logP += Math.log(rho_i*(1-r_iplus1)*q_iplus1);

            } else {

                logP += Math.log(psi_i*(1-r_i));

            }

        } else {
            // Internal node

            logP = getSingleTypeSubtreeLogLikelihood(subtreeRoot.getChild(0), t_node, A, B)
                    + getSingleTypeSubtreeLogLikelihood(subtreeRoot.getChild(1), t_node, A, B);

            logP += Math.log(2*lambda_i);

            double q_i = get_q_i(A[i], B[i], t_i, t_node);
            logP -= Math.log(q_i);

            if (Utils.equalWithPrecision(t_i, t_node)) {
                double q_iplus1 = i + 1 < parameterization.getTotalIntervalCount()
                        ? get_q_i(A[i+1], B[i+1], parameterization.getIntervalEndTimes()[i+1], t_node)
                        : 1.0;

                logP += 2*Math.log((1-rho_i)*q_iplus1);
            } else {

                logP += 2*Math.log(q_i);
            }
        }

        // Compute contributions from intervals along edge

        while (i >= 0 && Utils.greaterThanWithPrecision(parameterization.getIntervalEndTimes()[i], timeOfSubtreeRootEdgeTop)) {

            if (Utils.lessThanWithPrecision(parameterization.getIntervalEndTimes()[i], t_node)) {
                double q_iplus1 = i + 1 < parameterization.getTotalIntervalCount()
                        ? get_q_i(A[i + 1], B[i + 1],
                        parameterization.getIntervalEndTimes()[i + 1],
                        parameterization.getIntervalEndTimes()[i])
                        : 1.0;

                logP += Math.log((1-parameterization.getRhoValues()[i][0])*q_iplus1);
            }

            i -= 1;
        }

        return logP;
    }

    /* StateNode implementation */

    @Override
    public List<String> getArguments() {
        return null;
    }


    @Override
    public List<String> getConditions() {
        return null;
    }

    @Override
    public void sample(State state, Random random) {
    }

    @Override
    public boolean requiresRecalculation() {
        return true;
    }

    @Override
    public void store() {
        super.store();

        System.arraycopy(startTypeProbs, 0, storedStartTypeProbs, 0, parameterization.getNTypes());
    }

    @Override
    public void restore() {
        super.restore();

        double[] tmp = storedStartTypeProbs;
        startTypeProbs = storedStartTypeProbs;
        storedStartTypeProbs = tmp;
    }
}
