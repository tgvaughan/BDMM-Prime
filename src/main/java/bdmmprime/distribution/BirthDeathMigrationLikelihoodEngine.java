/*
 * Copyright (c) 2017-2026 ETH Zürich
 *
 * This file is part of bdmm-prime.
 *
 * bdmm-prime is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or (at
 * your option) any later version.
 *
 * bdmm-prime is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with bdmm-prime. If not, see <https://www.gnu.org/licenses/>.
 */

package bdmmprime.distribution;

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
