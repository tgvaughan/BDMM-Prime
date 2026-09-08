package bdmmflow.flowSystems;

import bdmmflow.intervals.Interval;
import org.apache.commons.math3.linear.RealMatrix;

import java.util.List;

public interface IFlowODESystem {
    IFlow calculateFlowIntegral(
            String initialMatrixStrategy,
            boolean parallelize
    );
}
