package bdmm.distributions;

import bdmm.parameterization.Parameterization;
import bdmm.util.Utils;
import beast.core.Citation;
import beast.core.Description;
import beast.core.Input;
import beast.core.State;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.evolution.speciation.SpeciesTreeDistribution;
import beast.evolution.tree.Node;
import beast.evolution.tree.TraitSet;
import beast.evolution.tree.Tree;
import beast.evolution.tree.TreeInterface;
import beast.util.HeapSort;
import org.apache.commons.math3.ode.FirstOrderIntegrator;
import org.apache.commons.math3.ode.nonstiff.DormandPrince54Integrator;

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

@Description("This model implements a multi-deme version of the BirthDeathSkylineModel with discrete locations and migration events among demes. " +
        "This implementation also works with sampled ancestor trees.")
public class BirthDeathMigrationDistribution extends SpeciesTreeDistribution {

    public Input<Parameterization> parameterizationInput = new Input<>("parameterization",
            "BDMM parameterization",
            Input.Validate.REQUIRED);

    public Input<RealParameter> frequenciesInput = new Input<>("frequencies",
            "The frequencies for each type",
            Input.Validate.REQUIRED);

    public Input<TraitSet> tiptypes = new Input<>("tiptypes",
            "trait information for initializing traits " +
                    "(like node types/locations) in the tree");

    public Input<String> typeLabel = new Input<>("typeLabel",
            "type label in tree for initializing traits " +
                    "(like node types/locations) in the tree");

    public Input<IntegerParameter> tipTypeArray = new Input<>("tipTypeArray",
            "integer array of traits (like node types/locations) " +
                    "in the tree, index corresponds to node number in tree");

    public Input<Integer> maxEvaluations = new Input<>("maxEvaluations",
            "The maximum number of evaluations for ODE solver",
            1000000);

    public Input<Boolean> conditionOnSurvival = new Input<>("conditionOnSurvival",
            "condition on at least one survival? Default true.",
            true);

    public Input<Double> relativeTolerance = new Input<>("relTolerance",
            "relative tolerance for numerical integration",
            1e-7);

    public Input<Double> absoluteTolerance = new Input<>("absTolerance",
            "absolute tolerance for numerical integration",
            1e-100 /*Double.MIN_VALUE*/);

    public Input<Boolean> checkRho = new Input<>("checkRho",
            "check if rho is set if multiple tips are " +
                    "given at present (default true)",
            true);

    public Input<Boolean> parallelizeInput = new Input<>(
            "parallelize",
            "is the calculation parallelized on sibling subtrees " +
                    "or not (default true)",
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

    public Input<Boolean> storeNodeTypes = new Input<>("storeNodeTypes",
            "store tip node types? this assumes that tip types cannot " +
                    "change (default false)", false);

    private int[] nodeStates;

    //	private final boolean debug = false;
    private final boolean debug = true;

    private double[] rootTypeProbs, storedRootTypeProbs;
    private boolean[] isRhoTip;

    private Parameterization parameterization;

    private P0System P;
    private P0GeSystem PG;
    private static Double minstep, maxstep;


    private static double[][] pInitialConditions;

    private static boolean isParallelizedCalculation;
    private static double minimalProportionForParallelization;
    private final static double globalPrecisionThreshold = 1e-10;
    private double parallelizationThreshold;
    private static ThreadPoolExecutor pool;

    private double[] weightOfNodeSubTree;

    private TreeInterface tree;

    @Override
    public void initAndValidate() {
        if ((tiptypes.get() == null ? 0 : 1) + (typeLabel.get() == null ? 0 : 1) + (tipTypeArray.get() == null ? 0 : 1) != 1)
            throw new RuntimeException("Tip types need to be specified exactly once using either tiptypes OR typeLabel OR tipTypeArray.");

        parameterization = parameterizationInput.get();

        tree = treeInput.get();

        double freqSum = 0;
        for (double f : frequenciesInput.get().getValues()) freqSum += f;
        if (Math.abs(1.0 - freqSum) > 1e-10)
            throw new RuntimeException("Error: frequencies must add up to 1 but currently add to " + freqSum + ".");


        int nLeaves = tree.getLeafNodeCount();

        weightOfNodeSubTree = new double[nLeaves * 2];

        isParallelizedCalculation = parallelizeInput.get();
        minimalProportionForParallelization = minimalProportionForParallelizationInput.get();

        if (isParallelizedCalculation) executorBootUp();

        setupIntegrators();

        if (storeNodeTypes.get()) {

            nodeStates = new int[nLeaves];

            for (Node node : tree.getExternalNodes()) {
                nodeStates[node.getNr()] = getNodeType(node, true);
            }
        }

        rootTypeProbs = new double[parameterization.getNTypes()];
        storedRootTypeProbs = new double[parameterization.getNTypes()];

        // Determine which, if any, of the leaf ages correspond exactly to
        // rho sampling times.
        isRhoTip = new boolean[tree.getLeafNodeCount()];
        for (int nodeNr = 0; nodeNr < tree.getLeafNodeCount(); nodeNr++) {
            isRhoTip[nodeNr] = false;
            double nodeTime = parameterization.getOrigin() - tree.getNode(nodeNr).getHeight();
            for (double rhoSampTime : parameterization.getRhoSamplingTimes()) {
                if (rhoSampTime == nodeTime) {
                    isRhoTip[nodeNr] = true;
                    break;
                }
            }
        }
    }

    @Override
    public double calculateTreeLogLikelihood(TreeInterface tree) {

        Node root = tree.getRoot();

        if (parameterization.getOrigin() < tree.getRoot().getHeight()) {
            logP = Double.NEGATIVE_INFINITY;
            return logP;
        }

        // update the threshold for parallelization
        //TODO only do it if tree shape changed
        updateParallelizationThreshold();

        updateInitialConditionsForP(tree);

        double probNoSample = 0;
        if (conditionOnSurvival.get()) {

            double[] noSampleExistsProp = pInitialConditions[pInitialConditions.length - 1];
            if (debug) System.out.println("\nnoSampleExistsProp = "
                    + noSampleExistsProp[0] + ", " + noSampleExistsProp[1]);

            for (int rootType = 0; rootType < parameterization.getNTypes(); rootType++) {
                probNoSample += frequenciesInput.get().getArrayValue(rootType) * noSampleExistsProp[rootType];
            }

            if (probNoSample < 0 || probNoSample > 1)
                return Double.NEGATIVE_INFINITY;

        }

        P0GeState finalP0Ge = calculateSubtreeLikelihood(root, 0,
                parameterization.getOrigin() - tree.getRoot().getHeight(), PG, 0);

        if (debug) System.out.print("Final state: " + finalP0Ge);

        SmallNumber PrSN = new SmallNumber(0);
        for (int rootType = 0; rootType < parameterization.getNTypes(); rootType++) {

            SmallNumber jointProb = finalP0Ge
                    .ge[rootType]
                    .scalarMultiplyBy(frequenciesInput.get().getArrayValue(rootType));

            if (jointProb.getMantissa() > 0) {
                rootTypeProbs[rootType] = jointProb.log();
                PrSN = PrSN.addTo(jointProb);
            } else {
                rootTypeProbs[rootType] = Double.NEGATIVE_INFINITY;
            }
        }

        // Normalize root type probs:
        for (int rootType = 0; rootType < parameterization.getNTypes(); rootType++) {
            rootTypeProbs[rootType] -= PrSN.log();
            rootTypeProbs[rootType] = Math.exp(rootTypeProbs[rootType]);
        }

        if (conditionOnSurvival.get()) {
            PrSN = PrSN.scalarMultiplyBy(1 / (1 - probNoSample));
        }

        logP = PrSN.log();

        if (debug) System.out.println("\nlogP = " + logP);

        // TGV: Why is this necessary?
        if (Double.isInfinite(logP)) logP = Double.NEGATIVE_INFINITY;

        // Weird correction:
        // TGV: Why is this only applied when the removal prob is less than 1??
        if (parameterization.getRemovalProbs()[0][0] != 1.0) {
            int internalNodeCount = tree.getLeafNodeCount() - ((Tree) tree).getDirectAncestorNodeCount() - 1;
            logP += Math.log(2) * internalNodeCount;
        }

        return logP;
    }

    private int getNodeType(Node node, Boolean init) {

        try {

            if (!storeNodeTypes.get() || init) {

                int nodestate = -1;

                if (tiptypes.get() != null) {

                    nodestate = (int) tiptypes.get().getValue((node.getID()));

                } else if (typeLabel.get() != null) {

                    Object d = node.getMetaData(typeLabel.get());

                    if (d instanceof Integer) nodestate = (Integer) d;
                    else if (d instanceof Double)
                        nodestate = ((Double) d).intValue();
                    else if (d instanceof int[]) nodestate = ((int[]) d)[0];
                    else if (d instanceof String)
                        nodestate = Integer.valueOf((String) d);
                    else
                        throw new RuntimeException("Error interpreting as type index: " + d);

                } else {

                    nodestate = (int) tipTypeArray.get().getArrayValue((node.getNr()));
                }

                if (nodestate < -1)
                    throw new ConstraintViolatedException("State assignment failed.");

                return nodestate;

            } else return nodeStates[node.getNr()];

        } catch (Exception e) {
            throw new ConstraintViolatedException("Something went wrong with the assignment of types to the nodes (node ID=" + node.getID() + "). Please check your XML file!");
        }
    }

    public int getLeafStateForLogging(Node node, long sampleNr) {
        if (!node.isLeaf()) {
            throw new IllegalArgumentException("Node should be a leaf");
        }
        return getNodeType(node, sampleNr == 0);
    }


    /**
     * @param node    Node below edge.
     * @param tTop    Time of start (top) of edge.
     * @param tBottom Time of end (bottom) of edge.
     * @param system  Object describing ODEs to integrate.
     * @return State at top of edge.
     */
    private P0GeState calculateSubtreeLikelihood(Node node, double tTop, double tBottom, P0GeSystem system, int depth) {

        if (debug) {
            debugMessage("*** Evaluating subtree for node " + node +
                            " and edge between times " + tTop + " and " + tBottom + " ...",
                    depth);
        }

        P0GeState state = new P0GeState(parameterization.getNTypes());

        int intervalIdx = Utils.getIntervalIndex(tBottom, system.intervalStartTimes);

        if (node.isLeaf()) { // sampling event

            int nodeType = getNodeType(node, false);

            if (nodeType == -1) { //unknown state

                //TODO test if SA model case is properly implemented (not tested!)
                for (int type = 0; type < parameterization.getNTypes(); type++) {

                    if (!isRhoTip[node.getNr()]) {
                        state.ge[type] = new SmallNumber(
                                (system.r[intervalIdx][type] + pInitialConditions[node.getNr()][type] * (1 - system.r[intervalIdx][type]))
                                        * system.s[intervalIdx][type]);
                        // with SA: ψ_i(r + (1 − r)p_i(τ))
                    } else {
                        state.ge[type] = new SmallNumber(
                                (system.r[intervalIdx][type] + pInitialConditions[node.getNr()][type]
                                        / (1 - system.rho[intervalIdx][type]) * (1 - system.r[intervalIdx][type]))
                                        * system.rho[intervalIdx][type]);
                    }
                }
            } else {

                if (!isRhoTip[node.getNr()]) {
                    state.ge[nodeType] = new SmallNumber(
                            (system.r[intervalIdx][nodeType] + pInitialConditions[node.getNr()][nodeType]
                                    * (1 - system.r[intervalIdx][nodeType]))
                                    * system.s[intervalIdx][nodeType]);
                    // with SA: ψ_i(r + (1 − r)p_i(τ))
                } else {
                    state.ge[nodeType] = new SmallNumber(
                            (system.r[intervalIdx][nodeType] + pInitialConditions[node.getNr()][nodeType]
                                    / (1 - system.rho[intervalIdx][nodeType]) * (1 - system.r[intervalIdx][nodeType]))
                                    * system.rho[intervalIdx][nodeType]);
                }

            }

            // Incorporate pre-evaluated p0 values into state
            System.arraycopy(pInitialConditions[node.getNr()], 0, state.p0, 0, system.nTypes);

            if (debug) debugMessage("Sampling at time " + tBottom, depth);

        } else if (node.getChildCount() == 2) {  // birth / infection event or sampled ancestor

            if (node.getChild(0).isDirectAncestor() || node.getChild(1).isDirectAncestor()) {   // found a sampled ancestor

                int childIndex = 0;

                if (node.getChild(childIndex).isDirectAncestor())
                    childIndex = 1;

                P0GeState g = calculateSubtreeLikelihood(
                        node.getChild(childIndex), tBottom,
                        parameterization.getOrigin() - node.getChild(childIndex).getHeight(),
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

//					System.out.println("SA but not rho sampled");

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

                P0GeState childState1 = new P0GeState();
                P0GeState childState2 = new P0GeState();

                // evaluate if the next step in the traversal should be split between one new thread and the currrent thread and run in parallel.

                if (isParallelizedCalculation
                        && weightOfNodeSubTree[node.getChild(indexFirstChild).getNr()] > parallelizationThreshold
                        && weightOfNodeSubTree[node.getChild(indexSecondChild).getNr()] > parallelizationThreshold) {

                    try {
                        // start a new thread to take care of the second subtree
                        Future<P0GeState> secondChildTraversal = pool.submit(
                                new TraversalServiceUncoloured(node.getChild(indexSecondChild), tBottom,
                                        parameterization.getOrigin() - node.getChild(indexSecondChild).getHeight(),
                                        depth + 1));

                        childState1 = calculateSubtreeLikelihood(
                                node.getChild(indexFirstChild), tBottom,
                                parameterization.getOrigin() - node.getChild(indexFirstChild).getHeight(),
                                system, depth + 1);
                        childState2 = secondChildTraversal.get();
                    } catch (InterruptedException | ExecutionException e) {
                        e.printStackTrace();

                        System.exit(1);
                    }

                } else {
                    childState1 = calculateSubtreeLikelihood(node.getChild(
                            indexFirstChild), tBottom,
                            parameterization.getOrigin() - node.getChild(indexFirstChild).getHeight(),
                            system, depth + 1);
                    childState2 = calculateSubtreeLikelihood(node.getChild(indexSecondChild), tBottom,
                            parameterization.getOrigin() - node.getChild(indexSecondChild).getHeight(),
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

        P0GeState statePrime = getG(node, tTop, state, system);

        if (debug)
            debugMessage("State at top of edge: " + statePrime + "\n", depth);

        return statePrime;
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

    /**
     * Used to indicate that the state assignment went wrong.
     */
    protected class ConstraintViolatedException extends RuntimeException {
        public ConstraintViolatedException(String s) {
            super(s);
        }

    }

    class TraversalServiceUncoloured extends TraversalService {

        int depth;

        public TraversalServiceUncoloured(Node root, double from, double to, int depth) {
            super(root, from, to);
            this.depth = depth;
        }

        @Override
        protected P0GeState calculateSubtreeLikelihoodInThread() {
            return calculateSubtreeLikelihood(rootSubtree, from, to, PG, depth);
        }

    }

    /**
     * @return retrieve current set of root type probabilities.
     */
    public double[] getRootTypeProbs() {
        return rootTypeProbs;
    }


    /**
     * Perform integration on differential equations p
     *
     * @param tEnd
     * @return
     */
    private P0State getP(double tStart, double tEnd, P0State state, P0System system) {

        if (Math.abs(system.origin - tEnd) < globalPrecisionThreshold
                || Math.abs(tStart - tEnd) < globalPrecisionThreshold)
            return state;

        double thisTime = tStart;
        int thisInterval = Utils.getIntervalIndex(thisTime, system.intervalStartTimes);
        int endInterval = Utils.getIntervalIndex(tEnd, system.intervalStartTimes);

        while (thisInterval > endInterval) {

            double nextTime = system.intervalStartTimes[thisInterval];

            if (nextTime < thisTime) {
                p_integrator.integrate(system, thisTime, state.p0, nextTime, state.p0);

                for (int i = 0; i < system.nTypes; i++)
                    state.p0[i] *= (1 - system.rho[thisInterval][i]);
            }

            thisTime = nextTime;
            thisInterval -= 1;
        }

        p_integrator.integrate(P, thisTime, state.p0, tEnd, state.p0); // solve diffEquationOnP, store solution in y

        // TO DO
        // check that both rateChangeTimes are really overlapping
        // but really not sure that this is enough, i have to build appropriate tests
        if (Math.abs(tEnd - system.intervalStartTimes[endInterval]) < globalPrecisionThreshold) {
            for (int i = 0; i < system.nTypes; i++)
                state.p0[i] = 1 - system.rho[endInterval][i];
        }

        return state;
    }

    /**
     * Integrate state along edge above baseNode until time tTop according to system.
     *
     * @param baseNode node at base of edge
     * @param tTop     time at top of edge
     * @param state    ODE variables at bottom of edge
     * @param system   ODE system to integrate
     * @return
     */
    public P0GeState getG(Node baseNode, double tTop, P0GeState state, P0GeSystem system) {


        // pgScaled contains the set of initial conditions scaled made to fit
        // the requirements on the values 'double' can represent. It also
        // contains the factor by which the numbers were multiplied.
        ScaledNumbers pgScaled = state.getScaledState();

        double thisTime = system.origin - baseNode.getHeight();
        int thisInterval = Utils.getIntervalIndex(thisTime, system.intervalStartTimes);
        int endInterval = Utils.getIntervalIndex(tTop, system.intervalStartTimes);
        double oneMinusRho;

        while (thisInterval > endInterval) {
            double nextTime = system.intervalStartTimes[thisInterval];

            if (nextTime < thisTime) {
                pgScaled = safeIntegrate(system, thisTime-1e-10, pgScaled, nextTime+1e-10);

                state.setFromScaledState(pgScaled.getEquation(), pgScaled.getScalingFactor());

                for (int i = 0; i < parameterization.getNTypes(); i++) {
                    oneMinusRho = 1 - system.rho[thisInterval][i];
                    state.p0[i] *= oneMinusRho;
                    state.ge[i] = state.ge[i].scalarMultiplyBy(oneMinusRho);
                }

                // 'rescale' the results of the last integration to prepare for the next integration step
                pgScaled = state.getScaledState();
            }

            thisTime = nextTime;
            thisInterval -= 1;
        }

         // solve PG , store solution temporarily integrationResults
        pgScaled = safeIntegrate(system, thisTime-1e-10, pgScaled, tTop+1e-10);

        // 'unscale' values in integrationResults so as to retrieve accurate values after the integration.
        state.setFromScaledState(pgScaled.getEquation(), pgScaled.getScalingFactor());

        return state;
    }


    /**
     * Perform an initial traversal of the tree to get the 'weights' (sum of all its edges lengths) of all sub-trees
     * Useful for performing parallelized calculations on the tree.
     * The weights of the subtrees tell us the depth at which parallelization should stop, so as to not parallelize on subtrees that are too small.
     * Results are stored in 'weightOfNodeSubTree' array
     *
     * @param tree
     */
    public void getAllSubTreesWeights(TreeInterface tree) {
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
     * @param node
     * @return
     */
    public double getSubTreeWeight(Node node) {

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

    private FirstOrderIntegrator p_integrator;

    private void setupIntegrators() {   // set up ODE's and integrators

        //TODO set these minstep and maxstep to be a class field
        if (minstep == null) minstep = parameterization.getOrigin() * 1e-100;
        if (maxstep == null) maxstep = parameterization.getOrigin() / 10;

        P = new P0System(parameterization);
        PG = new P0GeSystem(parameterization);

        p_integrator = new DormandPrince54Integrator(minstep, maxstep, absoluteTolerance.get(), relativeTolerance.get());
    }

    /**
     * Perform the integration of PG with initial conds in pgScaled between to and from
     * Use an adaptive-step-size integrator
     * "Safe" because it divides the integration interval in two
     * if the interval is (arbitrarily) judged to be too big to give reliable results
     *
     * @param PG
     * @param tStart
     * @param pgScaled
     * @param tEnd
     * @return result of integration
     */
    private ScaledNumbers safeIntegrate(P0GeSystem PG, double tStart, ScaledNumbers pgScaled, double tEnd) {

        // if the integration interval is too small, nothing is done (to prevent infinite looping)
        if (Math.abs(tEnd - tStart) < globalPrecisionThreshold /*(T * 1e-20)*/)
            return pgScaled;

        //TODO make threshold a class field
        if (PG.origin > 0 && Math.abs(tEnd - tStart) > PG.origin / 6) {
            pgScaled = safeIntegrate(PG, tStart, pgScaled, tEnd + (tStart - tEnd) / 2);
            pgScaled = safeIntegrate(PG, tEnd + (tStart - tEnd) / 2, pgScaled, tEnd);
        } else {

            //setup of the relativeTolerance and absoluteTolerance input of the adaptive integrator
            //TODO set these two as class fields
            double relativeToleranceConstant = 1e-7;
            double absoluteToleranceConstant = 1e-100;
            double[] absoluteToleranceVector = new double[2 * PG.nTypes];
            double[] relativeToleranceVector = new double[2 * PG.nTypes];

            for (int i = 0; i < PG.nTypes; i++) {
                absoluteToleranceVector[i] = absoluteToleranceConstant;
                if (pgScaled.getEquation()[i + PG.nTypes] > 0) { // adapt absoluteTolerance to the values stored in pgScaled
                    absoluteToleranceVector[i + PG.nTypes] = Math.max(1e-310, pgScaled.getEquation()[i + PG.nTypes] * absoluteToleranceConstant);
                } else {
                    absoluteToleranceVector[i + PG.nTypes] = absoluteToleranceConstant;
                }
                relativeToleranceVector[i] = relativeToleranceConstant;
                relativeToleranceVector[i + PG.nTypes] = relativeToleranceConstant;
            }

            double[] integrationResults = new double[pgScaled.getEquation().length];
            int a = pgScaled.getScalingFactor(); // store scaling factor
            int n = pgScaled.getEquation().length / 2; // dimension of the ODE system


            FirstOrderIntegrator integrator = new DormandPrince54Integrator(minstep, maxstep, absoluteToleranceVector, relativeToleranceVector);
            integrator.integrate(PG, tStart, pgScaled.getEquation(), tEnd, integrationResults); // perform the integration step

            double[] pConditions = new double[n];
            SmallNumber[] geConditions = new SmallNumber[n];
            for (int i = 0; i < n; i++) {
                pConditions[i] = integrationResults[i];
                geConditions[i] = new SmallNumber(integrationResults[i + n]);
            }
            pgScaled = (new P0GeState(pConditions, geConditions)).getScaledState();
            pgScaled.augmentFactor(a);
        }

        return pgScaled;
    }


    /**
     * Compute all initial conditions for all future integrations on p0 equations.
     *
     * @param tree
     */
    public void updateInitialConditionsForP(TreeInterface tree) {

        int leafCount = tree.getLeafNodeCount();
        double[] leafTimes = new double[leafCount];
        int[] indicesSortedByLeafTime = new int[leafCount];

        for (int i = 0; i < leafCount; i++) { // get all leaf times
            leafTimes[i] = parameterization.getOrigin() - tree.getNode(i).getHeight();
            indicesSortedByLeafTime[i] = i;
        }

        HeapSort.sort(leafTimes, indicesSortedByLeafTime);
        //"sort" sorts in ascending order, so we have to be careful since the
        // integration starts from the leaves at time T and goes up to the
        // root at time 0 (or >0)


        double t = leafTimes[indicesSortedByLeafTime[leafCount - 1]];

        P0State p0State = new P0State(PG.nTypes);

        int startIndex = Utils.getIntervalIndex(parameterization.getOrigin(),
                parameterization.getIntervalStartTimes());
        for (int i = 0; i < parameterization.getNTypes(); i++) {
            p0State.p0[i] = (1 - parameterization.getRhoValues()[startIndex][i]);    // initial condition: y_i[T]=1-rho_i
        }

        pInitialConditions = new double[leafCount + 1][P.nTypes];
        System.arraycopy(getP(P.origin, t, p0State, P).p0, 0,
                pInitialConditions[indicesSortedByLeafTime[leafCount-1]], 0, P.nTypes);
        double tprev = t;

        if (leafCount > 1) {
            for (int i = leafCount - 2; i > -1; i--) {
                t = leafTimes[indicesSortedByLeafTime[i]];

                //If the next higher leaf is actually at the same height, store previous results and skip iteration
                if (Math.abs(t - tprev) < globalPrecisionThreshold) {
                    tprev = t;
                    pInitialConditions[indicesSortedByLeafTime[i]] = pInitialConditions[indicesSortedByLeafTime[i + 1]];
                    continue;
                } else {
					/* TODO the integration performed in getP is done before all
                       the other potentially-parallelized getG, so should not
                       matter that it has its own integrator, but if it does
                       (or to simplify the code), take care of passing an integrator
                       as a local variable. */
                    System.arraycopy(getP(tprev, t, p0State, P).p0, 0,
                            pInitialConditions[indicesSortedByLeafTime[i]], 0, P.nTypes);
                    tprev = t;
                }

            }
        }

        System.arraycopy(getP(tprev, 0, p0State, P).p0, 0,
                pInitialConditions[leafCount], 0, P.nTypes);
    }

    static void executorBootUp() {
        ExecutorService executor = Executors.newCachedThreadPool();
        pool = (ThreadPoolExecutor) executor;
    }

    static void executorShutdown() {
        pool.shutdown();
    }

    abstract class TraversalService implements Callable<P0GeState> {

        protected Node rootSubtree;
        protected double from;
        protected double to;
        protected P0GeSystem PG;
        protected FirstOrderIntegrator pg_integrator;

        public TraversalService(Node root, double from, double to) {
            this.rootSubtree = root;
            this.from = from;
            this.to = to;
            this.setupODEs();
        }

        private void setupODEs() {  // set up ODE's and integrators

            //TODO set minstep and maxstep to be PiecewiseBDDistr fields
            if (minstep == null)
                minstep = parameterization.getOrigin() * 1e-100;
            if (maxstep == null) maxstep = parameterization.getOrigin() / 10;

            PG = new P0GeSystem(parameterization);

            pg_integrator = new DormandPrince54Integrator(minstep, maxstep, absoluteTolerance.get(), relativeTolerance.get());
        }

        abstract protected P0GeState calculateSubtreeLikelihoodInThread();

        @Override
        public P0GeState call() throws Exception {
            // traverse the tree in a potentially-parallelized way
            return calculateSubtreeLikelihoodInThread();
        }
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

        for (int i = 0; i < parameterization.getNTypes(); i++)
            storedRootTypeProbs[i] = rootTypeProbs[i];
    }

    @Override
    public void restore() {
        super.restore();

        double[] tmp = storedRootTypeProbs;
        rootTypeProbs = storedRootTypeProbs;
        storedRootTypeProbs = tmp;
    }


}
