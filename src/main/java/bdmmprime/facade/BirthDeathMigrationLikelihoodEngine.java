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

    /* BEAST StateNode methods */
    void store();
    void restore();
    void accept();
    boolean requiresRecalculation();
    boolean isStochastic();

}
