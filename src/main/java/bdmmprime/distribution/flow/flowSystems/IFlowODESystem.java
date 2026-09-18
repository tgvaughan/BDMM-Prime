package bdmmprime.distribution.flow.flowSystems;

public interface IFlowODESystem {
    IFlow calculateFlowIntegral(
            InitialMatrixStrategy initialMatrixStrategy,
            boolean parallelize
    );
}
