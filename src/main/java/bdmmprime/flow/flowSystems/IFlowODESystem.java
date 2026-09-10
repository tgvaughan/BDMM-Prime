package bdmmprime.flow.flowSystems;

import bdmmprime.flow.intervals.Interval;
import org.apache.commons.math3.linear.RealMatrix;

import java.util.List;

public interface IFlowODESystem {
    IFlow calculateFlowIntegral(
            InitialMatrixStrategy initialMatrixStrategy,
            boolean parallelize
    );
}
