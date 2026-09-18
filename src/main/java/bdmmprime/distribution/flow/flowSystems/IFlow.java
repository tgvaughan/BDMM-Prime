package bdmmprime.distribution.flow.flowSystems;

public interface IFlow {
    IntegrationResult integrateUsingFlow(
            double timeStart,
            double timeEnd,
            double[] endState
    );
}
