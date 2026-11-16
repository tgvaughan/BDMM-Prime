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

package bdmmprime.distribution.flow;

import bdmmprime.parameterization.*;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import beast.base.spec.domain.NonNegativeReal;
import beast.base.spec.domain.UnitInterval;
import beast.base.spec.inference.parameter.RealScalarParam;
import beast.base.spec.inference.parameter.RealVectorParam;
import beast.base.spec.inference.parameter.SimplexParam;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.ValueSource;

import static org.junit.jupiter.api.Assertions.assertEquals;

/**
 * Tests that the flow engine reuses the likelihoods of unchanged subtrees and that reusing them
 * does not change the likelihood.
 */
public class SubtreeCachingTest {

    private static final String NEWICK = "(((a[&type=0]:1.0,b[&type=1]:0.7):1.2,(c[&type=1]:0.5,d[&type=0]:1.4):0.9):1.1,"
            + "((e[&type=0]:0.8,f[&type=1]:1.3):1.0,(g[&type=1]:0.6,h[&type=0]:1.1):0.7):1.5);";

    private static final String SAMPLED_ANCESTOR_NEWICK = "(((a[&type=0]:1.0,b[&type=1]:0.0):1.2,c[&type=1]:0.9):1.1,"
            + "(d[&type=0]:0.8,e[&type=1]:1.3):1.5);";

    private static final double TOLERANCE = 1e-10;

    @ParameterizedTest
    @ValueSource(booleans = {false, true})
    public void testFirstCalculationUpdatesAllNodes(boolean parallelize) {
        Tree tree = new TreeParser(NEWICK, false);
        BirthDeathMigrationDistribution density = this.buildDensity(tree, 1.0, parallelize);

        density.calculateLogP();

        assertEquals(tree.getNodeCount(), density.getNumUpdatedNodes());
    }

    @ParameterizedTest
    @ValueSource(booleans = {false, true})
    public void testUnchangedTreeUpdatesNothing(boolean parallelize) {
        Tree tree = new TreeParser(NEWICK, false);
        BirthDeathMigrationDistribution density = this.buildDensity(tree, 1.0, parallelize);

        double logP = density.calculateLogP();

        assertEquals(logP, density.calculateLogP());
        assertEquals(0, density.getNumUpdatedNodes());
    }

    @ParameterizedTest
    @ValueSource(booleans = {false, true})
    public void testDisabledCachingUpdatesAllNodes(boolean parallelize) {
        Tree tree = new TreeParser(NEWICK, false);
        BirthDeathMigrationDistribution density = this.buildDensity(tree, 1.0, parallelize, false);

        double logP = density.calculateLogP();

        assertEquals(logP, density.calculateLogP());
        assertEquals(tree.getNodeCount(), density.getNumUpdatedNodes());
    }

    @ParameterizedTest
    @ValueSource(booleans = {false, true})
    public void testInternalNodeHeightChange(boolean parallelize) {
        Tree tree = new TreeParser(NEWICK, false);
        BirthDeathMigrationDistribution density = this.buildDensity(tree, 1.0, parallelize);
        density.calculateLogP();

        Node node = this.getLeaf(tree, "a").getParent();
        node.setHeight(node.getHeight() + 0.1);

        double logP = density.calculateLogP();

        // the node, its two children and its two ancestors
        assertEquals(5, density.getNumUpdatedNodes());
        assertEquals(this.buildDensity(tree, 1.0, parallelize).calculateLogP(), logP, TOLERANCE);
    }

    @ParameterizedTest
    @ValueSource(booleans = {false, true})
    public void testLeafHeightChange(boolean parallelize) {
        Tree tree = new TreeParser(NEWICK, false);
        BirthDeathMigrationDistribution density = this.buildDensity(tree, 1.0, parallelize);
        density.calculateLogP();

        Node leaf = this.getLeaf(tree, "c");
        leaf.setHeight(leaf.getHeight() + 0.2);

        double logP = density.calculateLogP();

        // the leaf and its three ancestors
        assertEquals(4, density.getNumUpdatedNodes());
        assertEquals(this.buildDensity(tree, 1.0, parallelize).calculateLogP(), logP, TOLERANCE);
    }

    @ParameterizedTest
    @ValueSource(booleans = {false, true})
    public void testSubtreeSwap(boolean parallelize) {
        Tree tree = new TreeParser(NEWICK, false);
        BirthDeathMigrationDistribution density = this.buildDensity(tree, 1.0, parallelize);
        density.calculateLogP();

        // swap the leaf c with its uncle, the parent of a and b

        Node child = this.getLeaf(tree, "c");
        Node parent = child.getParent();
        Node grandparent = parent.getParent();
        Node uncle = this.getLeaf(tree, "a").getParent();

        parent.removeChild(child);
        grandparent.removeChild(uncle);
        parent.addChild(uncle);
        grandparent.addChild(child);

        double logP = density.calculateLogP();

        // the two swapped nodes, the parent, the grandparent and the root
        assertEquals(5, density.getNumUpdatedNodes());
        assertEquals(this.buildDensity(tree, 1.0, parallelize).calculateLogP(), logP, TOLERANCE);
    }

    @ParameterizedTest
    @ValueSource(booleans = {false, true})
    public void testSampledAncestors(boolean parallelize) {
        Tree tree = new TreeParser(SAMPLED_ANCESTOR_NEWICK, false);
        BirthDeathMigrationDistribution density = this.buildDensity(tree, 0.5, parallelize);
        density.calculateLogP();

        // the sampled ancestor is handled as part of its parent
        assertEquals(tree.getNodeCount() - 1, density.getNumUpdatedNodes());

        Node leaf = this.getLeaf(tree, "c");
        leaf.setHeight(leaf.getHeight() + 0.2);

        double logP = density.calculateLogP();

        // the leaf and its two ancestors
        assertEquals(3, density.getNumUpdatedNodes());
        assertEquals(this.buildDensity(tree, 0.5, parallelize).calculateLogP(), logP, TOLERANCE);
    }

    @ParameterizedTest
    @ValueSource(booleans = {false, true})
    public void testRestoreAfterRejectedChange(boolean parallelize) {
        Tree tree = new TreeParser(NEWICK, false);
        BirthDeathMigrationDistribution density = this.buildDensity(tree, 1.0, parallelize);

        double logP = density.calculateLogP();
        density.accept();

        // propose a change

        Node node = this.getLeaf(tree, "a").getParent();
        double height = node.getHeight();

        density.store();
        node.setHeight(height + 0.1);
        density.calculateLogP();

        // reject the change

        node.setHeight(height);
        density.restore();

        assertEquals(logP, density.calculateLogP());
        assertEquals(0, density.getNumUpdatedNodes());

        // propose and accept another change

        density.store();
        node.setHeight(height - 0.1);
        double logPAfterChange = density.calculateLogP();
        density.accept();

        assertEquals(5, density.getNumUpdatedNodes());
        assertEquals(this.buildDensity(tree, 1.0, parallelize).calculateLogP(), logPAfterChange, TOLERANCE);
    }

    private Node getLeaf(Tree tree, String id) {
        for (Node leaf : tree.getExternalNodes()) {
            if (id.equals(leaf.getID())) return leaf;
        }

        throw new IllegalArgumentException("Unknown leaf " + id);
    }

    /**
     * Builds a flow density with subtree caching for the given tree.
     */
    private BirthDeathMigrationDistribution buildDensity(Tree tree, double removalProb, boolean parallelize) {
        return this.buildDensity(tree, removalProb, parallelize, true);
    }

    /**
     * Builds a flow density with two types and migration for the given tree. When parallelization is
     * enabled, the subtree size needed to start parallel calculations is lowered such that the small
     * test trees are traversed in parallel.
     */
    private BirthDeathMigrationDistribution buildDensity(Tree tree, double removalProb, boolean parallelize, boolean useSubtreeCaching) {
        Parameterization parameterization = new CanonicalParameterization();
        parameterization.initByName(
                "typeSet", new TypeSet(2),
                "processLength", new RealScalarParam<>(5.0, NonNegativeReal.INSTANCE),
                "birthRate", new SkylineVectorParameter(
                        null,
                        new RealVectorParam<>(new double[] {2.0, 1.5}, NonNegativeReal.INSTANCE)),
                "deathRate", new SkylineVectorParameter(
                        null,
                        new RealVectorParam<>(new double[] {1.0, 0.8}, NonNegativeReal.INSTANCE)),
                "birthRateAmongDemes", new SkylineMatrixParameter(
                        null,
                        new RealVectorParam<>(new double[] {0.0}, NonNegativeReal.INSTANCE), 2),
                "migrationRate", new SkylineMatrixParameter(
                        null,
                        new RealVectorParam<>(new double[] {0.2, 0.7}, NonNegativeReal.INSTANCE)),
                "samplingRate", new SkylineVectorParameter(
                        null,
                        new RealVectorParam<>(new double[] {0.5, 0.3}, NonNegativeReal.INSTANCE)),
                "removalProb", new SkylineVectorParameter(
                        null,
                        new RealVectorParam<>(new double[] {removalProb}, UnitInterval.INSTANCE), 2));

        BirthDeathMigrationDistribution density = new BirthDeathMigrationDistribution();
        density.initByName(
                "parameterization", parameterization,
                "startTypePriorProbs", new SimplexParam(new double[] {0.5, 0.5}),
                "conditionOnSurvival", false,
                "tree", tree,
                "typeLabel", "type",
                "parallelize", parallelize,
                "minimalSubtreeSizeForParallelization", 1,
                "useSubtreeCaching", useSubtreeCaching
        );

        return density;
    }

}
