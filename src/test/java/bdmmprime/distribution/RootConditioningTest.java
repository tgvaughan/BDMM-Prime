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

import bdmmprime.parameterization.*;
import bdmmprime.testclasses.RealVectorParamFromString;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import beast.base.spec.domain.NonNegativeReal;
import beast.base.spec.domain.UnitInterval;
import beast.base.spec.inference.parameter.SimplexParam;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;

public class RootConditioningTest extends LikelihoodTestClass {

    @Test
    public void testSingleType() {

        Tree tree = new TreeParser(
                "((3[&type=0] : 1.5, 4[&type=1] : 0.5) : 1 , (1[&type=1] : 2, 2[&type=0] : 1) : 3);",
                false);

        Parameterization parameterization = new CanonicalParameterization();
        parameterization.initByName(
                "processLength", "5.0",
                "birthRate", new SkylineVectorParameter(
                        null,
                        new RealVectorParamFromString<>("2.0", NonNegativeReal.INSTANCE)),
                "deathRate", new SkylineVectorParameter(
                        null,
                        new RealVectorParamFromString<>("1.5", NonNegativeReal.INSTANCE)),
                "samplingRate", new SkylineVectorParameter(
                        null,
                        new RealVectorParamFromString<>("0.1", NonNegativeReal.INSTANCE)),

                "removalProb", new SkylineVectorParameter(
                        null,
                        new RealVectorParamFromString<>("1.0", UnitInterval.INSTANCE)));

        BirthDeathMigrationDistribution density = new BirthDeathMigrationDistribution();
        density.initByName("parameterization", parameterization,
                "conditionOnRoot", true,
                "tree", tree,
                "typeLabel", "type",
                "useAnalyticalSingleTypeSolution", false);

        double numericalLogP = density.calculateLogP();

        BirthDeathMigrationDistribution densityExact = new BirthDeathMigrationDistribution();
        densityExact.initByName("parameterization", parameterization,
                "conditionOnRoot", true,
                "tree", tree,
                "typeLabel", "type",
                "useAnalyticalSingleTypeSolution", true);

        double analyticalLogP = densityExact.calculateLogP();

        assertEquals(numericalLogP, analyticalLogP, 1e-5);
    }

    @Test
    public void testOriginException() {

        Tree tree = new TreeParser(
                "((3[&type=0] : 1.5, 4[&type=1] : 0.5) : 1 , (1[&type=1] : 2, 2[&type=0] : 1) : 3);",
                false);

        Parameterization parameterization = new CanonicalParameterization();
        parameterization.initByName(
                "processLength", "6.0",
                "birthRate", new SkylineVectorParameter(
                        null,
                        new RealVectorParamFromString<>("2.0", NonNegativeReal.INSTANCE)),
                "deathRate", new SkylineVectorParameter(
                        null,
                        new RealVectorParamFromString<>("1.5", NonNegativeReal.INSTANCE)),
                "samplingRate", new SkylineVectorParameter(
                        null,
                        new RealVectorParamFromString<>("0.1", NonNegativeReal.INSTANCE)),

                "removalProb", new SkylineVectorParameter(
                        null,
                        new RealVectorParamFromString<>("1.0", UnitInterval.INSTANCE)));

        BirthDeathMigrationDistribution density = new BirthDeathMigrationDistribution();

        assertThrows(RuntimeException.class,
                () -> density.initByName("parameterization",parameterization,
                        "conditionOnRoot", true,
                        "tree", tree,
                        "typeLabel", "type",
                        "useAnalyticalSingleTypeSolution", true));

        assertThrows(RuntimeException.class,
                () -> density.initByName("parameterization",parameterization,
                        "conditionOnRoot", true,
                        "tree", tree,
                        "typeLabel", "type",
                        "useAnalyticalSingleTypeSolution", false));
    }

    @Test
    public void testSAAvoidance() {

        Tree tree = new TreeParser(
                "(4[&type=1] : 0.0,((3[&type=0] : 1.5, 5[&type=1] : 0.5) : 1 , (1[&type=1] : 2, 2[&type=0] : 1) : 3):0.5):0.0;",
                false);

        Parameterization parameterization = new CanonicalParameterization();
        parameterization.initByName(
                "processLength", "5.5",
                "birthRate", new SkylineVectorParameter(
                        null,
                        new RealVectorParamFromString<>("2.0", NonNegativeReal.INSTANCE)),
                "deathRate", new SkylineVectorParameter(
                        null,
                        new RealVectorParamFromString<>("1.5", NonNegativeReal.INSTANCE)),
                "samplingRate", new SkylineVectorParameter(
                        null,
                        new RealVectorParamFromString<>("0.1", NonNegativeReal.INSTANCE)),

                "removalProb", new SkylineVectorParameter(
                        null,
                        new RealVectorParamFromString<>("0.0", UnitInterval.INSTANCE)));

        BirthDeathMigrationDistribution density = new BirthDeathMigrationDistribution();

        density.initByName("parameterization",parameterization,
                "conditionOnRoot", true,
                "tree", tree,
                "typeLabel", "type",
                "useAnalyticalSingleTypeSolution", false);
        assertEquals(Double.NEGATIVE_INFINITY, density.calculateLogP());

        density.initByName("parameterization",parameterization,
                "conditionOnRoot", true,
                "tree", tree,
                "typeLabel", "type",
                "useAnalyticalSingleTypeSolution", true);
        assertEquals(Double.NEGATIVE_INFINITY, density.calculateLogP());

    }

    @Test
    public void testBasicConditioning() {

        Tree tree = new TreeParser(
                "((3[&type=0] : 1.5, 4[&type=1] : 0.5) : 1 , (1[&type=1] : 2, 2[&type=0] : 1) : 3);",
                false);

        Parameterization parameterization = new CanonicalParameterization();
        parameterization.initByName(
                "processLength", "5.0",
                "typeSet", new TypeSet(2),
                "birthRate", new SkylineVectorParameter(
                        null,
                        new RealVectorParamFromString<>("2.0 1.5", NonNegativeReal.INSTANCE)),
                "deathRate", new SkylineVectorParameter(
                        null,
                        new RealVectorParamFromString<>("1.5 1.0", NonNegativeReal.INSTANCE)),
                "samplingRate", new SkylineVectorParameter(
                        null,
                        new RealVectorParamFromString<>("0.2 0.5", NonNegativeReal.INSTANCE)),
                "migrationRate", new SkylineMatrixParameter(
                        null,
                        new RealVectorParamFromString<>("0.2 0.1", NonNegativeReal.INSTANCE)),

                "removalProb", new SkylineVectorParameter(
                        null,
                        new RealVectorParamFromString<>("0.5", UnitInterval.INSTANCE), 2));

        BirthDeathMigrationDistribution density = new BirthDeathMigrationDistribution();
        density.initByName("parameterization", parameterization,
                "startTypePriorProbs", new SimplexParam(new double[] {0.5, 0.5}),
                "conditionOnRoot", true,
                "tree", tree,
                "typeLabel", "type",
                "parallelize", false);



        // TODO: Check!!
        assertEquals(-19.1786415, density.calculateLogP(), 1e-5);
    }
}
