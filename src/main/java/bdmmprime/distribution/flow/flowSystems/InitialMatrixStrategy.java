package bdmmprime.distribution.flow.flowSystems;

/**
 * The strategy used to pick the initial state (preconditioner) of the flow integration on each interval.
 */
public enum InitialMatrixStrategy {

    identity,
    random,
    average_inverse,

}
