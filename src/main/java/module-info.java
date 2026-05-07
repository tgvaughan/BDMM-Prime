open module bdmmprime {
    requires beast.pkgmgmt;
    requires beast.base;
    requires beast.fx;
    requires org.apache.commons.statistics.distribution;
    requires org.antlr.antlr4.runtime;
    requires com.google.common;
    requires commons.math3;
    requires javafx.controls;
    requires colt;

    exports bdmmprime.distribution;
    exports bdmmprime.mapping;
    exports bdmmprime.parameterization;
    exports bdmmprime.trajectories;
    exports bdmmprime.util;

    provides beast.base.core.BEASTInterface with
        bdmmprime.mapping.TypeMappedTree,
        bdmmprime.mapping.TypedTreeStatsLogger,
        bdmmprime.mapping.TypedNodeTreeLogger,
        bdmmprime.mapping.TypedTipTreeLogger,
        bdmmprime.mapping.TransitionTimeLogger,
        bdmmprime.mapping.TypeMappingDurationLogger,
        bdmmprime.mapping.SkylineNodeTreeLogger,
        bdmmprime.util.TipDatesFromTree,
        bdmmprime.util.priors.SmartZeroExcludingRealIID,
        bdmmprime.util.priors.ZeroExcludingRealIID,
        bdmmprime.util.priors.OUSkyGridPrior,
        bdmmprime.util.priors.SkyGridPrior,
        bdmmprime.util.priors.TipDatePrior,
        bdmmprime.util.priors.MultivariateLogNormal,
        bdmmprime.util.operators.SmartScaleOperator,
        bdmmprime.util.operators.ChangeTimeOperator,
        bdmmprime.util.operators.TipEdgeScaleOperator,
        bdmmprime.util.operators.TipDateInitialiser,
        bdmmprime.util.operators.TipDateOperator,
        bdmmprime.util.operators.SmartRealRandomWalkOperator,
        bdmmprime.util.operators.JointSkylineScaleOperator,
        bdmmprime.util.InitializedTraitSet,
        bdmmprime.util.Slice,
        bdmmprime.util.ProcessLength,
        bdmmprime.util.NodeAgeLogger,
        bdmmprime.util.OptionalLogger,
        bdmmprime.util.StartTypePriorProbsLogger,
        bdmmprime.util.StartTypePosteriorProbsLogger,
        bdmmprime.parameterization.TimedParameter,
        bdmmprime.parameterization.CanonicalParameterization,
        bdmmprime.parameterization.SkylineMatrixParameter,
        bdmmprime.parameterization.SkylineVectorParameter,
        bdmmprime.parameterization.EpiParameterizationMod,
        bdmmprime.parameterization.EpiParameterization,
        bdmmprime.parameterization.FBDParameterization,
        bdmmprime.parameterization.TypeSet,
        bdmmprime.distribution.BirthDeathMigrationDistribution,
        bdmmprime.trajectories.simulation.SimulatedTrajectoryLogger,
        bdmmprime.trajectories.simulation.UnconditionedTrajectoryLogger,
        bdmmprime.trajectories.simulation.SimulatedTree,
        bdmmprime.trajectories.simulation.UntypedTreeFromTypedTree,
        bdmmprime.trajectories.SampledTrajectory,
        bdmmprime.trajectories.TreeProbEstimateLogger,
        bdmmprime.trajectories.TrajSamplingDurationLogger;

    provides beastfx.app.inputeditor.InputEditor with
        bdmmprime.beauti.TimedParameterInputEditor,
        bdmmprime.beauti.ParameterizationInputEditor,
        bdmmprime.beauti.SkylineVectorInputEditor,
        bdmmprime.beauti.SkylineMatrixInputEditor,
        bdmmprime.beauti.TypeTraitSetInputEditor,
        bdmmprime.beauti.ProcessLengthInputEditor,
        bdmmprime.beauti.SimplexParamInputEditor;
}
