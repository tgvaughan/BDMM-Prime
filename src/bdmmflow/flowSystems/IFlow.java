package bdmmflow.flowSystems;

import org.apache.commons.math3.linear.RealMatrix;

public interface IFlow {
    IntegrationResult integrateUsingFlow(
            double timeStart,
            double timeEnd,
            double[] endState
    );
}
