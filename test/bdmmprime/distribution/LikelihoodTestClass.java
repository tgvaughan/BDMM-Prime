/*
 * Copyright (C) 2019-2025 ETH Zurich
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package bdmmprime.distribution;

import beast.base.evolution.tree.Tree;
import org.apache.commons.math.special.Gamma;

public abstract class LikelihoodTestClass {

    /**
     * The original tests were developed assuming BDSKY/BDMM-like behaviour, i.e. return an oriented
     * tree probability unless r!=1 in which case return an un-oriented and unlabeled tree probability.
     * In contrast, BDMM-Prime always returns a labeled tree probability.
     *
     * This method exists to convert BDSKY/BDMM test probabilities to be labeled tree probabilities,
     * allowing comparison with BDMM-Prime.
     *
     * @param density BDMM-prime probability density object
     * @return conversion factor
     */
    protected double labeledTreeConversionFactor(BirthDeathMigrationDistribution density) {
        Tree tree = (Tree)density.treeInput.get();
        boolean SAmodel = density.parameterizationInput.get().getRemovalProbs()[0][0] != 1.0;
        double factor = - Gamma.logGamma(tree.getLeafNodeCount() +1);

        if (!SAmodel)
            factor += Math.log(2) * (tree.getLeafNodeCount() - tree.getDirectAncestorNodeCount() - 1);

        return factor;
    }
}
