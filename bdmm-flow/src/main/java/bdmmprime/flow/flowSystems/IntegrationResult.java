package bdmmflow.flowSystems;

/**
 * Stores an integration result as well as a scaling factor.
 * The actual result is result * exp(logScalingFactor).
 */
public record IntegrationResult(double[] result, double logScalingFactor) {
}
