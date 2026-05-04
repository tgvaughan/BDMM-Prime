package bdmmflow;

import bdmmflow.extinctionSystem.ExtinctionProbabilities;
import bdmmflow.extinctionSystem.ExtinctionProbabilitiesODESystem;
import bdmmflow.flowSystems.*;
import bdmmflow.intervals.Interval;
import bdmmflow.intervals.IntervalODESystem;
import bdmmflow.intervals.IntervalUtils;
import bdmmflow.utils.Result;
import bdmmflow.utils.Utils;
import bdmmprime.parameterization.Parameterization;
import beast.base.core.*;
import beast.base.evolution.speciation.SpeciesTreeDistribution;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.TraitSet;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeInterface;
import beast.base.inference.parameter.RealParameter;
import org.apache.commons.math.special.Gamma;
import org.apache.commons.math3.exception.NumberIsTooSmallException;
import org.apache.commons.math3.linear.SingularMatrixException;
import org.apache.commons.math3.ode.ContinuousOutputModel;

import java.util.Arrays;
import java.util.List;
import java.util.Objects;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.CompletionException;
import java.util.concurrent.ForkJoinPool;
import java.util.stream.DoubleStream;

@Citation(value = "Kuehnert D, Stadler T, Vaughan TG, Drummond AJ. (2016). " +
        "A General and Efficient Algorithm for the Likelihood of Diversification and Discrete-Trait Evolutionary Models, \n" +
        "Systematic Biology, Volume 69, Issue 3, May 2020, Pages 545–556."
        , DOI = "10.1093/sysbio/syz055", year = 2020, firstAuthorSurname = "Louca")

@Description("This model implements a multi-deme version of the BirthDeathSkylineModel \" +\n" +
        "        \"with discrete locations and migration events among demes. \" +\n" +
        "        \"This implementation uses the Flow representation of the probability  \" +\n" +
        "        \"ODE for better performance. \" +\n" +
        "        \"It can be used as a drop-in replacement of the BDMM-Prime package."
)
public class BirthDeathMigrationDistribution extends SpeciesTreeDistribution {

    public Input<Parameterization> parameterizationInput = new Input<>(
            "parameterization",
            "BDMM parameterization",
            Input.Validate.REQUIRED
    );

    public Input<Function> finalSampleOffsetInput = new Input<>(
            "finalSampleOffset",
            "If provided, the difference in time between the final sample and the end of the BD process.",
            new RealParameter("0.0")
    );

    public Input<RealParameter> startTypePriorProbsInput = new Input<>(
            "startTypePriorProbs",
            "The prior probabilities for the type of the first individual",
            new RealParameter("1.0")
    );

    public Input<String> typeLabelInput = new Input<>(
            "typeLabel",
            "Attribute key used to specify sample trait values in tree."
    );

    public Input<TraitSet> typeTraitSetInput = new Input<>("typeTraitSet",
            "Trait set specifying sample trait values.");

    public Input<Boolean> conditionOnSurvivalInput = new Input<>("conditionOnSurvival",
            "Condition on at least one surviving lineage. (Default true.)",
            true);

    public Input<Boolean> conditionOnRootInput = new Input<>("conditionOnRoot",
            "Condition on root age, not time of origin.", false);

    public Input<Double> relativeToleranceInput = new Input<>(
            "relTolerance",
            "Relative tolerance for numerical integration.",
            1e-7
    );

    public Input<Double> absoluteToleranceInput = new Input<>(
            "absTolerance",
            "Absolute tolerance for numerical integration.",
            1e-100
    );

    public Input<String> initialMatrixStrategyInput = new Input<>(
            "initialMatrixStrategy",
            "The strategy to use to get the initial flow state. Either 'random', 'heuristic', or 'identity'.",
            "identity"
    );

    public Input<Boolean> useInverseFlowInput = new Input<>(
            "useInverseFlow",
            "Whether to use the inverse flow algorithm. It is faster, but can lead to higher numerical instability.",
            false
    );

    public Input<Integer> seedInput = new Input<>(
            "seed",
            "The random seed used in the analysis.",
            3215
    );

    public Input<Boolean> parallelizeInput = new Input<>(
            "parallelize",
            "Whether or not parallelize the computation.",
            true
    );

    public Input<Integer> minimalSubtreeSizeForParallelizationInput = new Input<>(
            "minimalSubtreeSizeForParallelization",
            "the minimal absolute size the two children " +
                    "subtrees of a node must have to start parallel " +
                    "calculations on the children. ",
            64);

    public Input<Double> maxConditioningNumberInput = new Input<>(
            "maxConditioningNumber",
            "The maximal conditioning number to reach until an interval is split.",
            1e8
    );

    public Input<Boolean> useLoucaPennellIntervalsInput = new Input<>(
            "useLoucaPennellIntervals",
            "Whether to use the interval heruistic introduced by Louca and Pennell.",
            false
    );

    private Parameterization parameterization;

    private String initialMatrixStrategy;

    private double finalSampleOffset;
    private TreeInterface tree;
    private String typeLabel;
    private TraitSet typeTraitSet;

    double[] startTypePriorProbs;

    boolean conditionOnRoot;
    boolean conditionOnSurvival;

    double absoluteTolerance;
    double relativeTolerance;

    double maxConditioningNumber;
    boolean useLoucaPennellIntervals;

    boolean useInverseFlow;
    int seed;

    boolean parallelize;
    double minimalProportionForParallelization = 0.05;
    int minimalSubtreeSizeForParallelization;

    ForkJoinPool forkJoinPool;
    int[] subtreeSizes;
    double parallelizeSubtreeSizeThreshold;

    int numTypes;

    double[] logScalingFactors;
    boolean[] isRhoSampled;

    int totalNumEvaluations = 0;
    int numEvaluationsSinceReset = 0;
    int numFailedEvaluationsSinceReset = 0;
    int numDeviations = 0;
    double sumDeviation = 0;

    bdmmprime.distribution.BirthDeathMigrationDistribution bdmmPrime;

    IFlow storedFlow;
    ExtinctionProbabilities storedExtinctionProbabilities;

    IFlow currentFlow;
    ExtinctionProbabilities currentExtinctionProbabilities;

    @Override
    public void initAndValidate() {
        // unpack input values

        this.parameterization = this.parameterizationInput.get();
        this.finalSampleOffset = this.finalSampleOffsetInput.get().getArrayValue();
        this.tree = this.treeInput.get();
        this.typeLabel = this.typeLabelInput.get();
        this.typeTraitSet = this.typeTraitSetInput.get();
        this.startTypePriorProbs = this.startTypePriorProbsInput.get().getDoubleValues();
        this.conditionOnRoot = this.conditionOnRootInput.get();
        this.conditionOnSurvival = this.conditionOnSurvivalInput.get();
        this.absoluteTolerance = this.absoluteToleranceInput.get();
        this.relativeTolerance = this.relativeToleranceInput.get();
        this.numTypes = this.parameterization.getNTypes();
        this.initialMatrixStrategy = this.initialMatrixStrategyInput.get();
        this.seed = this.seedInput.get();
        this.parallelize = this.parallelizeInput.get();
        this.minimalSubtreeSizeForParallelization = minimalSubtreeSizeForParallelizationInput.get();
        this.useInverseFlow = this.useInverseFlowInput.get();
        this.maxConditioningNumber = this.maxConditioningNumberInput.get();
        this.useLoucaPennellIntervals = this.useLoucaPennellIntervalsInput.get();

        // reset cache

        this.storedExtinctionProbabilities = null;
        this.currentExtinctionProbabilities = null;
        this.storedFlow = null;
        this.currentFlow = null;

        // validate type label

        if (this.numTypes != 1 && this.typeLabel == null && this.typeTraitSet == null) {
            throw new RuntimeException(
                    "Error: For models with >1 type, typeLabel or typeTraitSet must be specified."
            );
        }

        // validate start type prior probabilities

        if (this.startTypePriorProbs.length != numTypes) {
            throw new RuntimeException(
                    "Error: dimension of start type prior probabilities must match number of types."
            );
        }

        double probSum = 0.0;
        for (double f : this.startTypePriorProbs) {
            probSum += f;
        }
        if (!bdmmprime.util.Utils.equalWithPrecision(probSum, 1.0)) {
            throw new RuntimeException(
                    "Error: start type prior probabilities must add up to 1 but currently add to %f.".formatted(probSum)
            );
        }

        // check that we don't have birth events with two different birth types

//        if (this.parameterization.hasCrossBirthRates3()) {
//            throw new RuntimeException(
//                    "Error: BDMM-Flow does not support birth events with two different child types. Use BDMM-Prime instead."
//            );
//        }

        // initialize utils

        this.logScalingFactors = new double[this.tree.getNodeCount()];
        this.initializeIsRhoSampled();
        this.forkJoinPool = new ForkJoinPool();

        // initialize bdmm prime as a fallback if we detect numerical issues

        this.bdmmPrime = new bdmmprime.distribution.BirthDeathMigrationDistribution();
        this.bdmmPrime.initByName(
                "tree", this.treeInput.get(),
                "parameterization", this.parameterizationInput.get(),
                "finalSampleOffset", this.finalSampleOffsetInput.get(),
                "startTypePriorProbs", this.startTypePriorProbsInput.get(),
                "typeTraitSet", this.typeTraitSetInput.get(),
                "typeLabel", this.typeLabelInput.get(),
                "conditionOnSurvival", this.conditionOnSurvivalInput.get(),
                "conditionOnRoot", this.conditionOnRootInput.get()
        );
    }

    /**
     * Initializes the isRhoSampled array. The array contains a boolean for every node indicating
     * whether it was rho-sampled or not.
     */
    private void initializeIsRhoSampled() {
        this.isRhoSampled = new boolean[this.tree.getLeafNodeCount()];

        for (int i = 0; i < this.tree.getLeafNodeCount(); i++) {
            this.isRhoSampled[i] = false;

            double nodeTime = this.parameterization.getNodeTime(this.tree.getNode(i), this.finalSampleOffset);

            for (double rhoSamplingTime : this.parameterization.getRhoSamplingTimes()) {
                if (bdmmprime.util.Utils.equalWithPrecision(nodeTime, rhoSamplingTime)) {
                    this.isRhoSampled[i] = true;
                    break;
                }
            }
        }
    }

    /**
     * Initialized the `subtreeSizes` array with the number of nodes in the subtree of
     * each node.
     */
    private void initializeSubtreeSizes() {
        this.subtreeSizes = new int[this.tree.getNodeCount()];
        initializeSubtreeSizes(tree.getRoot());
        this.parallelizeSubtreeSizeThreshold = Math.max(
                this.subtreeSizes[tree.getRoot().getNr()] * this.minimalProportionForParallelization,
                this.minimalSubtreeSizeForParallelization
        );
    }

    /**
     * Initialized the `subtreeSizes` array with the number of nodes in the subtree of
     * `node`.
     */
    private int initializeSubtreeSizes(Node node) {
        int subtreeSize = 1;
        for (Node child : node.getChildren()) {
            subtreeSize += this.initializeSubtreeSizes(child);
        }
        this.subtreeSizes[node.getNr()] = subtreeSize;
        return subtreeSize;
    }

    /**
     * Calculates the log tree likelihood.
     *
     * @param dummyTree a dummyTree that is not considered, a BEAST implementation detail.
     * @return the calculated log tree likelihood for the tree specified in the given parameterization.
     */
    @Override
    public double calculateTreeLogLikelihood(TreeInterface dummyTree) {
        warnAboutNumericalIssuesIfNecessary();
        this.numEvaluationsSinceReset++;
        this.totalNumEvaluations++;

        // validate input values

        if (inputValuesHaveZeroDensity()) {
            return Double.NEGATIVE_INFINITY;
        }

        // set up subtrees for parallelization

        this.initializeSubtreeSizes();

        // set up intervals

        List<Interval> intervals = IntervalUtils.getIntervals(this.parameterization);

        // integrate over the extinction probabilities ODE and the flow ODE

        ExtinctionProbabilities extinctionProbabilities = null;
        IFlow flow = null;
        try {
            extinctionProbabilities = this.calculateExtinctionProbabilities(intervals);
            flow = this.calculateFlow(intervals, extinctionProbabilities);
        } catch (NumberIsTooSmallException | SingularMatrixException | IllegalStateException e) {
            this.numFailedEvaluationsSinceReset++;
            this.resetCache();
            return this.bdmmPrime.calculateTreeLogLikelihood(dummyTree);
        }

        // recursively traverse the tree to calculate the root likelihood per state

        Node root = this.tree.getRoot();
        double[] rootLikelihoodPerState;

        try {
            rootLikelihoodPerState = this.calculateSubTreeLikelihood(
                    root,
                    0,
                    this.parameterization.getNodeTime(root, this.finalSampleOffset),
                    flow,
                    extinctionProbabilities
            );
        } catch (CompletionException | SingularMatrixException | IllegalStateException exception) {
            this.numFailedEvaluationsSinceReset++;
            return this.bdmmPrime.calculateTreeLogLikelihood(dummyTree);
        }

        // get tree likelihood by a weighted average of the root likelihood per state

        double treeLikelihood = 0.0;
        for (int i = 0; i < this.parameterization.getNTypes(); i++) {
            treeLikelihood += rootLikelihoodPerState[i] * this.startTypePriorProbs[i];
        }

        if (treeLikelihood <= 0) {
            return Double.NEGATIVE_INFINITY;
        };

        // consider different ways to condition the tree

        double conditionDensity = this.calculateConditionDensityFactor(extinctionProbabilities);
        treeLikelihood /= conditionDensity;

        // turn the likelihood into log likelihood and correct for scaling

        double logTreeLikelihood = Math.log(treeLikelihood) + this.logScalingFactors[root.getNr()];

        // convert from oriented to labeled tree likelihood

        int internalNodeCount = tree.getLeafNodeCount() - ((Tree) tree).getDirectAncestorNodeCount() - 1;
        logTreeLikelihood += Math.log(2) * internalNodeCount - Gamma.logGamma(tree.getLeafNodeCount() + 1);

        // periodically compare with BDMMPrime

        this.periodicallyCompareToBDMMPrime(dummyTree, logTreeLikelihood);

        return logTreeLikelihood;
    }

    /**
     * Returns true if the input values are invalid or have a density of 0.
     */
    private boolean inputValuesHaveZeroDensity() {
        // force update of parameterization if needed
        this.parameterization.getIntervalEndTimes();

        if (!this.parameterization.valuesAreValid()) {
            return true;
        }

        if (bdmmprime.util.Utils.lessThanWithPrecision(
                parameterization.getNodeTime(tree.getRoot(), this.finalSampleOffset),
                0)) {
            // tree MRCA older than the start of the process
            return true;
        }

        if (conditionOnRootInput.get() && tree.getRoot().isFake()) {
            // tree root is a sampled ancestor, but we're conditioning on a root birth.
            return true;
        }

        // all good :)

        return false;
    }

    /**
     * Periodically computes the BDMM-Prime likelihood and compares it to the one we get. Prints a warning in case
     * we detect a big deviation.
     */
    private void periodicallyCompareToBDMMPrime(TreeInterface dummyTree, double bdmmFlowLikelihood) {
        if (this.totalNumEvaluations < 1_000) return;
        if (this.totalNumEvaluations % 2_000 != 0) return;

        double bdmmPrimeLikelihood = this.bdmmPrime.calculateTreeLogLikelihood(dummyTree);
        double deviation = Math.abs(Math.abs(bdmmFlowLikelihood - bdmmPrimeLikelihood) / bdmmPrimeLikelihood);

        if (Double.isFinite(deviation)) {
            this.sumDeviation += deviation;
            this.numDeviations++;
        }

        if (deviation > 1e-2) {
            Log.warning("Found relative deviation of " + 100*deviation + "% to BDMM-Prime. Consider using BDMM-Prime instead of BDMM-Flow.");
        }

        // we log the deviation every 20_000 steps
        if (this.numDeviations % 10 == 0) {
            double meanDeviation = this.sumDeviation / this.numDeviations;
            Log.warning("Mean deviation was " + meanDeviation + " (sum " + this.sumDeviation + ", num " + this.numDeviations + ")");
        }
    }

    /**
     * Periodically computes the fraction of evaluations with numerical issues. If this failure rate exceeds a certain
     * value, we print a warning.
     */
    private void warnAboutNumericalIssuesIfNecessary() {
        if (this.totalNumEvaluations < 1_000) return;

        // we only update the values every 10_000 iterations except in the very beginning
        if (this.totalNumEvaluations % 1_000 != 0 || (10_000 < this.totalNumEvaluations && this.totalNumEvaluations % 10_000 != 0))
            return;

        double failureRate = 1.0 * this.numFailedEvaluationsSinceReset / this.numEvaluationsSinceReset;

        if (failureRate > 0.05) {
            Log.warning("Failure rate was " + failureRate + ". Consider using BDMM-Prime instead of BDMM-Flow.");
        }

        // reset counters

        this.numFailedEvaluationsSinceReset = 0;
        this.numEvaluationsSinceReset = 0;
    }

    /**
     * Integrates over the extinction probabilities ODE.
     *
     * @return a wrapper class that allows to query the extinction probabilities at any given time.
     */
    private ExtinctionProbabilities calculateExtinctionProbabilities(List<Interval> intervals) {
        if (!this.parameterization.isDirtyCalculation() && this.currentExtinctionProbabilities != null) {
            // the parameterization hasn't changed, which means the extinction probabilities are still the same
            return this.currentExtinctionProbabilities;
        }

        // initialize ODE system

        IntervalODESystem system = new ExtinctionProbabilitiesODESystem(
                this.parameterization,
                intervals,
                this.absoluteTolerance,
                this.relativeTolerance / 10.0
        );

        // create the initial states

        int endInterval = this.parameterization.getTotalIntervalCount() - 1;

        double[] initialState = new double[this.parameterization.getNTypes()];
        for (int i = 0; i < this.parameterization.getNTypes(); i++) {
            initialState[i] = 1 - this.parameterization.getRhoValues()[endInterval][i];
        }

        List<double[]> initialStates = List.of(initialState);

        // integrate

        ContinuousOutputModel[] integrationResults = system.integrateBackwards(
                initialStates, intervals, false, parallelize
        );

        ExtinctionProbabilities extinctionProbabilities = new ExtinctionProbabilities(integrationResults, this.parameterization.getNTypes());
        this.currentExtinctionProbabilities = extinctionProbabilities;
        return extinctionProbabilities;
    }

    /**
     * Precomputes the flow ODE.
     *
     * @param intervals
     * @param extinctionProbabilities the precomputed extinction probabilities.
     * @return a wrapper class that allows to query the flow at any given time.
     */
    private IFlow calculateFlow(List<Interval> intervals, ExtinctionProbabilities extinctionProbabilities) {
        if (!this.parameterization.isDirtyCalculation() && this.currentFlow != null) {
            // the parameterization hasn't changed, which means the flow is still the same
            return this.currentFlow;
        }

        IFlowODESystem system;

        // we use the sum of heights as seed, this makes it deterministic for identical trees
        DoubleStream heights = Arrays.stream(this.tree.getNodesAsArray()).mapToDouble(node -> node.getHeight());
        int heightSum = (int) Math.floor(10_000 * heights.sum());

        if (this.useInverseFlow) {
            system = new InverseFlowODESystem(
                    this.parameterization,
                    extinctionProbabilities,
                    intervals,
                    this.absoluteTolerance,
                    this.relativeTolerance,
                    heightSum,
                    this.maxConditioningNumber,
                    this.useLoucaPennellIntervals
            );
        } else {
            system = new FlowODESystem(
                    this.parameterization,
                    extinctionProbabilities,
                    intervals,
                    this.absoluteTolerance,
                    this.relativeTolerance,
                    heightSum,
                    this.maxConditioningNumber,
                    this.useLoucaPennellIntervals
            );
        }

        extinctionProbabilities.validateProbabilities(true);
        IFlow flow = system.calculateFlowIntegral(
                initialMatrixStrategy,
                this.parallelize
        );
        extinctionProbabilities.validateProbabilities(false);

        this.currentFlow = flow;
        return flow;
    }

    /**
     * Calculates the probability density factor due to the way the tree is conditioned.
     * <p>
     * See Tanja Stadler, How Can We Improve Accuracy of Macroevolutionary Rate Estimates?,
     * Systematic Biology, Volume 62, Issue 2, March 2013, Pages 321–329,
     * https://doi.org/10.1093/sysbio/sys073
     *
     * @param extinctionProbabilities the calculated extinction probabilities integral.
     * @return the factor due to tree conditioning.
     */
    private double calculateConditionDensityFactor(ExtinctionProbabilities extinctionProbabilities) {
        double conditionDensity = 0.0;

        if (this.conditionOnRoot) {
            double[] extinctionAtRoot = extinctionProbabilities.getProbability(0);

            int startInterval = this.parameterization.getIntervalIndex(0);

            for (int type1 = 0; type1 < parameterization.getNTypes(); type1++) {
                for (int type2 = 0; type2 < parameterization.getNTypes(); type2++) {
                    double rate = type1 == type2
                            ? parameterization.getBirthRates()[startInterval][type1]
                            : parameterization.getCrossBirthRates()[startInterval][type1][type2];

                    conditionDensity += rate * this.startTypePriorProbs[type1]
                            * (1 - extinctionAtRoot[type1])
                            * (1 - extinctionAtRoot[type2]);
                }
            }
        } else if (this.conditionOnSurvival) {
            double[] extinctionAtRoot = extinctionProbabilities.getProbability(0);

            for (int type = 0; type < parameterization.getNTypes(); type++) {
                conditionDensity += this.startTypePriorProbs[type] * (1 - extinctionAtRoot[type]);
            }
        } else {
            conditionDensity = 1.0;
        }

        return conditionDensity;
    }

    /**
     * Calculates the per-type likelihood of the subtree of the given node including the edge leading to the
     * node.
     */
    private double[] calculateSubTreeLikelihood(
            Node node,
            double timeEdgeStart,
            double timeEdgeEnd,
            IFlow flow,
            ExtinctionProbabilities extinctionProbabilities
    ) {
        double[] likelihoodEdgeEnd;

        if (node.isLeaf()) {
            likelihoodEdgeEnd = calculateLeafLikelihood(node, timeEdgeEnd, extinctionProbabilities);
        } else if (node.getChild(0).isDirectAncestor() || node.getChild(1).isDirectAncestor()) {
            likelihoodEdgeEnd = calculateDirectAncestorWithChildLikelihood(node, timeEdgeEnd, flow, extinctionProbabilities);
        } else {
            likelihoodEdgeEnd = calculateInternalEdgeLikelihood(node, timeEdgeEnd, flow, extinctionProbabilities);
        }

        IntegrationResult likelihoodEdgeStart = flow.integrateUsingFlow(
                timeEdgeStart,
                timeEdgeEnd,
                likelihoodEdgeEnd
        );

        // make sure all likelihoods are positive

        for (int i = 0; i < likelihoodEdgeStart.result().length; i++) {
            if (likelihoodEdgeStart.result()[i] < 0) {
                throw new IllegalStateException("Negative probability detected.");
            }
        }

        this.logScalingFactors[node.getNr()] += likelihoodEdgeStart.logScalingFactor();
        return likelihoodEdgeStart.result();
    }

    /**
     * Calculates the likelihood of a single leaf node including the edge leading to it.
     */
    private double[] calculateLeafLikelihood(
            Node node,
            double timeEdgeEnd,
            ExtinctionProbabilities extinctionProbabilities
    ) {

        int intervalEdgeEnd = this.parameterization.getIntervalIndex(timeEdgeEnd);
        double[] extinctionProbabilityEdgeEnd = extinctionProbabilities.getProbability(timeEdgeEnd);

        int nodeType = this.getNodeType(node);

        double[] likelihoodEdgeEnd = new double[this.parameterization.getNTypes()];

        if (parameterization.getTypeSet().isAmbiguousTypeIndex(nodeType)) {
            // this is an ambiguous state, we set the end likelihoods for all states
            // TODO: test if SA model case is properly implemented

            for (int type = 0; type < parameterization.getNTypes(); type++) {
                if (parameterization.getTypeSet().ambiguityExcludesType(nodeType, type))
                    continue;

                if (isRhoSampled[node.getNr()]) {
                    likelihoodEdgeEnd[type] = this.parameterization.getRhoValues()[intervalEdgeEnd][type];
                    // in this case, the other boundary conditions are handled by the ODE system in
                    // FlowODESystem and ExtinctionProbabilitiesODESystem
                } else {
                    likelihoodEdgeEnd[type] = this.parameterization.getSamplingRates()[intervalEdgeEnd][type] *
                            (this.parameterization.getRemovalProbs()[intervalEdgeEnd][type]
                                    + (1 - this.parameterization.getRemovalProbs()[intervalEdgeEnd][type])
                                    * extinctionProbabilityEdgeEnd[type]);
                }
            }

        } else {
            // we know the state and only set its end likelihood

            if (isRhoSampled[node.getNr()]) {
                likelihoodEdgeEnd[nodeType] = this.parameterization.getRhoValues()[intervalEdgeEnd][nodeType];
                // in this case, the other boundary conditions are handled by the ODE system in
                // FlowODESystem and ExtinctionProbabilitiesODESystem
            } else {
                likelihoodEdgeEnd[nodeType] = this.parameterization.getSamplingRates()[intervalEdgeEnd][nodeType] *
                        (this.parameterization.getRemovalProbs()[intervalEdgeEnd][nodeType]
                                + (1 - this.parameterization.getRemovalProbs()[intervalEdgeEnd][nodeType])
                                * extinctionProbabilityEdgeEnd[nodeType]);
            }
        }

        this.logScalingFactors[node.getNr()] = Utils.rescale(likelihoodEdgeEnd);

        return likelihoodEdgeEnd;
    }

    /**
     * Calculates the likelihood of a node that has a direct ancestor as a child; including the edge leading to it.
     */
    private double[] calculateDirectAncestorWithChildLikelihood(
            Node node,
            double timeEdgeEnd,
            IFlow flow,
            ExtinctionProbabilities extinctionProbabilities
    ) {
        int intervalEdgeEnd = this.parameterization.getIntervalIndex(timeEdgeEnd);

        // find the direct ancestor and the child

        Node directAncestor = node.getChild(0).isDirectAncestor() ?
                node.getChild(0) : node.getChild(1);
        Node child = node.getChild(0).isDirectAncestor() ?
                node.getChild(1) : node.getChild(0);

        // calculate the subtree likelihood of the child

        double[] likelihoodChild = this.calculateSubTreeLikelihood(
                child,
                timeEdgeEnd,
                this.parameterization.getNodeTime(child, this.finalSampleOffset),
                flow,
                extinctionProbabilities
        );

        // calculate the likelihood at the edge end

        double[] likelihoodEdgeEnd = new double[this.parameterization.getNTypes()];

        int daNodeType = this.getNodeType(directAncestor);

        if (parameterization.getTypeSet().isAmbiguousTypeIndex(daNodeType)) {
            // the direct ancestor is in an ambiguous state, we set the end likelihoods for all states
            // TODO: test if SA model case is properly implemented

            for (int type = 0; type < parameterization.getNTypes(); type++) {
                if (parameterization.getTypeSet().ambiguityExcludesType(daNodeType, type))
                    continue;

                if (isRhoSampled[directAncestor.getNr()]) {
                    likelihoodEdgeEnd[type] = this.parameterization.getRhoValues()[intervalEdgeEnd][type];
                } else {
                    likelihoodEdgeEnd[type] = this.parameterization.getSamplingRates()[intervalEdgeEnd][type];
                }

                likelihoodEdgeEnd[type] *= (1 - this.parameterization.getRemovalProbs()[intervalEdgeEnd][type])
                        * likelihoodChild[type];
            }

        } else {
            // we know the direct ancestor state and set the likelihood edge enf only for this type

            if (isRhoSampled[directAncestor.getNr()]) {
                likelihoodEdgeEnd[daNodeType] = this.parameterization.getRhoValues()[intervalEdgeEnd][daNodeType];
            } else {
                likelihoodEdgeEnd[daNodeType] = this.parameterization.getSamplingRates()[intervalEdgeEnd][daNodeType];
            }

            likelihoodEdgeEnd[daNodeType] *= (1 - this.parameterization.getRemovalProbs()[intervalEdgeEnd][daNodeType])
                    * likelihoodChild[daNodeType];
        }

        this.logScalingFactors[node.getNr()] = Utils.rescale(likelihoodEdgeEnd, this.logScalingFactors[child.getNr()]);

        return likelihoodEdgeEnd;
    }

    /**
     * Calculates the likelihood of the subtree of the given internal node including the edge leading to it.
     */
    private double[] calculateInternalEdgeLikelihood(
            Node node,
            double timeEdgeEnd,
            IFlow flow,
            ExtinctionProbabilities extinctionProbabilities
    ) {
        int intervalEdgeEnd = this.parameterization.getIntervalIndex(timeEdgeEnd);

        Node child1 = node.getChild(0);
        Node child2 = node.getChild(1);

        // calculate the likelihood of the two subtrees

        double[] likelihoodChild1;
        double[] likelihoodChild2;

        if (parallelize && subtreeSizes[child1.getNr()] > this.parallelizeSubtreeSizeThreshold
                && subtreeSizes[child2.getNr()] > this.parallelizeSubtreeSizeThreshold) {
            CompletableFuture<Result<double[]>> futureLikelihoodChild1 = CompletableFuture.supplyAsync(() ->
                    Result.of(() -> this.calculateSubTreeLikelihood(
                        child1,
                        timeEdgeEnd,
                        this.parameterization.getNodeTime(child1, this.finalSampleOffset),
                        flow,
                        extinctionProbabilities
                    )), this.forkJoinPool
            );
            likelihoodChild2 = this.calculateSubTreeLikelihood(
                    child2,
                    timeEdgeEnd,
                    this.parameterization.getNodeTime(child2, this.finalSampleOffset),
                    flow,
                    extinctionProbabilities
            );
            likelihoodChild1 = futureLikelihoodChild1.join().getOrThrow();
        } else {
            likelihoodChild1 = this.calculateSubTreeLikelihood(
                    child1,
                    timeEdgeEnd,
                    this.parameterization.getNodeTime(child1, this.finalSampleOffset),
                    flow,
                    extinctionProbabilities
            );
            likelihoodChild2 = this.calculateSubTreeLikelihood(
                    child2,
                    timeEdgeEnd,
                    this.parameterization.getNodeTime(child2, this.finalSampleOffset),
                    flow,
                    extinctionProbabilities
            );
        }

        // combine the child likelihoods to get the likelihood at the edge end

        double[] likelihoodEdgeEnd = new double[this.parameterization.getNTypes()];
        for (int i = 0; i < this.parameterization.getNTypes(); i++) {
            likelihoodEdgeEnd[i] += this.parameterization.getBirthRates()[intervalEdgeEnd][i] * (
                    likelihoodChild1[i] * likelihoodChild2[i]
            );

            for (int j = 0; j < parameterization.getNTypes(); j++) {
                if (i == j) {
                    continue;
                }

                likelihoodEdgeEnd[i] += 0.5 * this.parameterization.getCrossBirthRates()[intervalEdgeEnd][i][j] * (
                        likelihoodChild1[i] * likelihoodChild2[j] + likelihoodChild1[j] * likelihoodChild2[i]
                );
            }
        }

        this.logScalingFactors[node.getNr()] = Utils.rescale(
                likelihoodEdgeEnd,
                this.logScalingFactors[child1.getNr()] + this.logScalingFactors[child2.getNr()]
        );

        return likelihoodEdgeEnd;
    }

    /**
     * Returns the type of the given node.
     *
     * @param node the node to return the type of.
     * @return the type of the node.
     */
    private int getNodeType(Node node) {
        if (parameterization.getNTypes() == 1) {
            return 0;
        }

        String nodeTypeName;

        if (this.typeTraitSet != null)
            nodeTypeName = this.typeTraitSet.getStringValue(node.getID());
        else {
            Object metaData = node.getMetaData(this.typeLabel);
            if (metaData instanceof Double) {
                nodeTypeName = String.valueOf(Math.round((double) metaData));
            } else {
                nodeTypeName = metaData.toString();
            }
        }

        return parameterization.getTypeSet().getTypeIndex(nodeTypeName);
    }

    /** Caching **/

    @Override
    public boolean requiresRecalculation() {
        return true;
    }

    @Override
    public void accept() {
        this.storedExtinctionProbabilities = this.currentExtinctionProbabilities;
        this.storedFlow = this.currentFlow;
    }

    @Override
    public void restore() {
        this.currentExtinctionProbabilities = this.storedExtinctionProbabilities;
        this.currentFlow = this.storedFlow;
    }

    public void resetCache() {
        this.currentExtinctionProbabilities = null;
        this.storedExtinctionProbabilities = null;
        this.currentFlow = null;
        this.storedFlow = null;
    }

    @Override
    public boolean isStochastic() {
        return Objects.equals(this.initialMatrixStrategy, "random");
    }

}
