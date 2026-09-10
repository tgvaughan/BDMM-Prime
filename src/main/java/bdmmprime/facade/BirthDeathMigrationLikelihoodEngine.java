package bdmmprime.facade;

import beast.base.evolution.tree.TreeInterface;

/**
 * Interface for a concrete BDMM tree-likelihood engine.
 */
public interface BirthDeathMigrationLikelihoodEngine {

    /**
     * Computes the log tree likelihood for the tree the engine was configured with.
     */
    double calculateTreeLogLikelihood(TreeInterface tree);

    /**
     * Returns the posterior probabilities for the type of the first individual, as computed during
     * the most recent call to {@link #calculateTreeLogLikelihood(TreeInterface)}. The returned array
     * has one entry per type and is normalized to sum to one.
     */
    double[] getStartTypePosteriorProbs();

    /* BEAST StateNode methods */
    void store();
    void restore();
    void accept();
    boolean requiresRecalculation();
    boolean isStochastic();

}
