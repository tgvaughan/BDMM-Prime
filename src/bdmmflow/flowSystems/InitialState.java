package bdmmflow.flowSystems;

import org.apache.commons.math3.linear.RealMatrix;

/**
 * Stores an initial state (preconditioner) and its inverse matrix.
 */
public record InitialState(
        double[] initialState,
        RealMatrix inverse
) { }
