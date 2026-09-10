package bdmmprime.facade;

import bdmmprime.flow.flowSystems.InitialMatrixStrategy;
import bdmmprime.parameterization.Parameterization;
import beast.base.core.*;
import beast.base.evolution.speciation.SpeciesTreeDistribution;
import beast.base.evolution.tree.TraitSet;
import beast.base.evolution.tree.TreeInterface;
import beast.base.spec.domain.NonNegativeReal;
import beast.base.spec.inference.parameter.RealScalarParam;
import beast.base.spec.inference.parameter.SimplexParam;
import beast.base.spec.type.RealScalar;
import beast.base.spec.type.Simplex;

@Citation(value = "Kuehnert D, Stadler T, Vaughan TG, Drummond AJ. (2016). Phylodynamics with migration: " +
        "A computational framework to quantify population structure from genomic data. " +
        "Mol Biol Evol. 33(8):2102-2116.",
        DOI = "10.1093/molbev/msw064", year = 2016, firstAuthorSurname = "Kuehnert")

@Citation(value = "Louca S, Pennell MW. (2020). A general and efficient algorithm for the likelihood of " +
        "diversification and discrete-trait evolutionary models. " +
        "Syst Biol. 69(3):545-556.",
        DOI = "10.1093/sysbio/syz055", year = 2020, firstAuthorSurname = "Louca")

@Description("This model implements a multi-deme version of the BirthDeathSkylineModel " +
        "with discrete locations and migration events among demes. " +
        "This class supports both the classic implementation " +
        "(bdmmprime.distribution.BirthDeathMigrationDistribution) " +
        "and the flow implementation " +
        "(bdmmprime.flow.BirthDeathMigrationDistribution) depending on " +
        "the 'method' input.")
public class BirthDeathMigrationDistribution extends SpeciesTreeDistribution {

    // engine selection

    public Input<Method> methodInput = new Input<>(
            "method",
            "Which likelihood engine to use: 'auto', 'flow', or 'classic'. 'auto' uses the classic implementation " +
                    "when there is a single type or when there are fewer than 100 samples.",
            Method.auto,
            Method.values()
    );

    // inputs shared by both engines

    public Input<Parameterization> parameterizationInput = new Input<>(
            "parameterization",
            "BDMM parameterization",
            Input.Validate.REQUIRED
    );

    public Input<RealScalar<? extends NonNegativeReal>> finalSampleOffsetInput = new Input<>(
            "finalSampleOffset",
            "If provided, the difference in time between the final sample and the end of the BD process.",
            new RealScalarParam<>(0.0, NonNegativeReal.INSTANCE)
    );

    public Input<Simplex> startTypePriorProbsInput = new Input<>(
            "startTypePriorProbs",
            "The prior probabilities for the type of the first individual",
            new SimplexParam(new double[] {1.0})
    );

    public Input<String> typeLabelInput = new Input<>(
            "typeLabel",
            "Attribute key used to specify sample trait values in tree."
    );

    public Input<TraitSet> typeTraitSetInput = new Input<>(
            "typeTraitSet",
            "Trait set specifying sample trait values."
    );

    public Input<Boolean> conditionOnSurvivalInput = new Input<>(
            "conditionOnSurvival",
            "Condition on at least one surviving lineage. (Default true.)",
            true
    );

    public Input<Boolean> conditionOnRootInput = new Input<>(
            "conditionOnRoot",
            "Condition on root age, not time of origin.",
            false
    );

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

    public Input<Boolean> parallelizeInput = new Input<>(
            "parallelize",
            "Whether or not to parallelize the computation. (Default true.)",
            true
    );

    // inputs specific to the classic (bdmmprime.distribution) engine

    public Input<Boolean> useAnalyticalSingleTypeSolutionInput = new Input<>(
            "useAnalyticalSingleTypeSolution",
            "Classic engine only: use the analytical SABDSKY tree prior when the model has only one type.",
            true
    );

    public Input<Double> minimalProportionForParallelizationInput = new Input<>(
            "parallelizationFactor",
            "Classic engine only: the minimal relative size the two children subtrees of a node must have to " +
                    "start parallel calculations on the children. (default: 1/10).",
            1.0 / 10
    );

    public Input<String> savePartialLikelihoodsToFileInput = new Input<>(
            "savePartialLikelihoodsToFile",
            "Classic engine only: if provided, the name of a file to which a tree annotated with partial " +
                    "likelihoods will be written."
    );

    public Input<Boolean> saveIntegrationResultsInput = new Input<>(
            "storeIntegrationResults",
            "Classic engine only: if true, save results of ge(t) integration for use by other models.",
            false
    );

    // inputs specific to the flow (bdmmprime.flow) engine

    public Input<InitialMatrixStrategy> initialMatrixStrategyInput = new Input<>(
            "initialMatrixStrategy",
            "Flow engine only: strategy for the initial flow state.",
            InitialMatrixStrategy.average_inverse,
            InitialMatrixStrategy.values()
    );

    public Input<Boolean> useInverseFlowInput = new Input<>(
            "useInverseFlow",
            "Flow engine only: whether to use the inverse flow algorithm. It is faster, but can lead to higher " +
                    "numerical instability.",
            false
    );

    public Input<Integer> seedInput = new Input<>(
            "seed",
            "Flow engine only: the random seed used in the analysis.",
            3215
    );

    public Input<Integer> minimalSubtreeSizeForParallelizationInput = new Input<>(
            "minimalSubtreeSizeForParallelization",
            "Flow engine only: the minimal absolute size the two children subtrees of a node must have to start " +
                    "parallel calculations on the children.",
            64
    );

    public Input<Double> maxConditioningNumberInput = new Input<>(
            "maxConditioningNumber",
            "Flow engine only: the maximal conditioning number to reach until an interval is split.",
            1e8
    );

    public Input<Boolean> useLoucaPennellIntervalsInput = new Input<>(
            "useLoucaPennellIntervals",
            "Flow engine only: whether to use the interval heuristic introduced by Louca and Pennell.",
            false
    );

    private BirthDeathMigrationLikelihoodEngine engine;

    @Override
    public void initAndValidate() {
        super.initAndValidate();

        Method method = this.methodInput.get();
        Parameterization param = this.parameterizationInput.get();
        TreeInterface tree = this.treeInput.get();

        this.engine = switch (method) {
            case Method.auto -> {
                boolean useClassic = param.getNTypes() == 1 || tree.getLeafNodeCount() < 100;
                yield useClassic ? this.buildClassicEngine() : this.buildFlowEngine();
            }
            case Method.flow -> this.buildFlowEngine();
            case Method.classic -> this.buildClassicEngine();
        };
    }

    /**
     * Builds the flow engine, forwarding the shared inputs and the flow-specific ones.
     */
    private BirthDeathMigrationLikelihoodEngine buildFlowEngine() {
        Log.info("Using the flow implementation of BDMM-Prime.");

        bdmmprime.flow.BirthDeathMigrationDistribution impl = new bdmmprime.flow.BirthDeathMigrationDistribution();

        // shared inputs

        forward(impl.treeInput, this.treeInput, impl);
        forward(impl.parameterizationInput, this.parameterizationInput, impl);
        forward(impl.finalSampleOffsetInput, this.finalSampleOffsetInput, impl);
        forward(impl.startTypePriorProbsInput, this.startTypePriorProbsInput, impl);
        forward(impl.typeLabelInput, this.typeLabelInput, impl);
        forward(impl.typeTraitSetInput, this.typeTraitSetInput, impl);
        forward(impl.conditionOnSurvivalInput, this.conditionOnSurvivalInput, impl);
        forward(impl.conditionOnRootInput, this.conditionOnRootInput, impl);
        forward(impl.relativeToleranceInput, this.relativeToleranceInput, impl);
        forward(impl.absoluteToleranceInput, this.absoluteToleranceInput, impl);
        forward(impl.parallelizeInput, this.parallelizeInput, impl);

        // flow-specific inputs

        forward(impl.initialMatrixStrategyInput, this.initialMatrixStrategyInput, impl);
        forward(impl.useInverseFlowInput, this.useInverseFlowInput, impl);
        forward(impl.seedInput, this.seedInput, impl);
        forward(impl.minimalSubtreeSizeForParallelizationInput, this.minimalSubtreeSizeForParallelizationInput, impl);
        forward(impl.maxConditioningNumberInput, this.maxConditioningNumberInput, impl);
        forward(impl.useLoucaPennellIntervalsInput, this.useLoucaPennellIntervalsInput, impl);

        impl.initAndValidate();
        return impl;
    }

    /**
     * Builds the classic engine, forwarding the shared inputs and the classic-specific ones.
     */
    private BirthDeathMigrationLikelihoodEngine buildClassicEngine() {
        Log.info("Using the classic implementation of BDMM-Prime.");

        bdmmprime.distribution.BirthDeathMigrationDistribution impl = new bdmmprime.distribution.BirthDeathMigrationDistribution();

        // shared inputs

        forward(impl.treeInput, this.treeInput, impl);
        forward(impl.parameterizationInput, this.parameterizationInput, impl);
        forward(impl.finalSampleOffsetInput, this.finalSampleOffsetInput, impl);
        forward(impl.startTypePriorProbsInput, this.startTypePriorProbsInput, impl);
        forward(impl.typeLabelInput, this.typeLabelInput, impl);
        forward(impl.typeTraitSetInput, this.typeTraitSetInput, impl);
        forward(impl.conditionOnSurvivalInput, this.conditionOnSurvivalInput, impl);
        forward(impl.conditionOnRootInput, this.conditionOnRootInput, impl);
        forward(impl.relativeToleranceInput, this.relativeToleranceInput, impl);
        forward(impl.absoluteToleranceInput, this.absoluteToleranceInput, impl);
        forward(impl.parallelizeInput, this.parallelizeInput, impl);

        // classic-specific inputs

        forward(impl.useAnalyticalSingleTypeSolutionInput, this.useAnalyticalSingleTypeSolutionInput, impl);
        forward(impl.minimalProportionForParallelizationInput, this.minimalProportionForParallelizationInput, impl);
        forward(impl.savePartialLikelihoodsToFileInput, this.savePartialLikelihoodsToFileInput, impl);
        forward(impl.saveIntegrationResultsInput, this.saveIntegrationResultsInput, impl);

        impl.initAndValidate();
        return impl;
    }

    /**
     * Copies the value of a facade input into the matching engine input.
     * Inputs left unset on the facade (e.g. the optional {@code typeLabel}) are skipped so the engine keeps its default.
     */
    private static <T> void forward(Input<T> target, Input<T> source, BEASTInterface owner) {
        T value = source.get();
        if (value != null) {
            target.setTypedValue(value, owner);
        }
    }

    @Override
    public double calculateTreeLogLikelihood(TreeInterface dummyTree) {
        return this.engine.calculateTreeLogLikelihood(this.treeInput.get());
    }

    public double[] getStartTypePosteriorProbs() {
        return this.engine.getStartTypePosteriorProbs();
    }

    @Override
    public boolean requiresRecalculation() {
        return this.engine.requiresRecalculation();
    }

    @Override
    public void store() {
        super.store();
        this.engine.store();
    }

    @Override
    public void restore() {
        super.restore();
        this.engine.restore();
    }

    @Override
    public void accept() {
        super.accept();
        this.engine.accept();
    }

    @Override
    public boolean isStochastic() {
        return this.engine.isStochastic();
    }

}
