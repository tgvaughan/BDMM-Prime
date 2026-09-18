package bdmmprime.facade;

import bdmmprime.distribution.flow.flowSystems.InitialMatrixStrategy;
import bdmmprime.parameterization.*;
import bdmmprime.util.ProcessLength;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import beast.base.spec.domain.NonNegativeReal;
import beast.base.spec.domain.UnitInterval;
import beast.base.spec.inference.parameter.RealScalarParam;
import beast.base.spec.inference.parameter.RealVectorParam;
import beast.base.spec.inference.parameter.SimplexParam;
import org.apache.commons.math3.special.Gamma;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.MethodSource;

import java.util.Arrays;
import java.util.Collection;

import static org.junit.jupiter.api.Assertions.assertEquals;

/**
 * These tests were taken from the original BDMM-Prime distribution tests.
 */
public class BirthDeathMigrationLikelihoodTest {

    /**
     * Creates the different test parameterizations.
     * Format is { method, initial matrix, use inverse flow, parallelize }.
     */
    public static Collection<Object[]> data() {
        return Arrays.asList(new Object[][] {
            { Engine.classic, InitialMatrixStrategy.identity, false,  false },
            { Engine.classic, InitialMatrixStrategy.identity, false,  true },
            { Engine.flow, InitialMatrixStrategy.identity, false,  false },
            { Engine.flow, InitialMatrixStrategy.identity, false,  true },
            { Engine.flow, InitialMatrixStrategy.identity, true,  false },
            { Engine.flow, InitialMatrixStrategy.identity, true,  true },
            { Engine.flow, InitialMatrixStrategy.random,  false, false },
            { Engine.flow, InitialMatrixStrategy.random,  false, true },
            { Engine.flow, InitialMatrixStrategy.random, true,  false },
            { Engine.flow, InitialMatrixStrategy.random, true,  true },
            { Engine.flow, InitialMatrixStrategy.average_inverse,  false, false },
            { Engine.flow, InitialMatrixStrategy.average_inverse,  false, true },
            { Engine.flow, InitialMatrixStrategy.average_inverse, true,  false },
            { Engine.flow, InitialMatrixStrategy.average_inverse, true,  true }
        });
    }

    /**
     * The original tests were developed assuming BDSKY/BDMM-like behaviour, i.e. return an oriented
     * tree probability unless r!=1 in which case return an un-oriented and unlabeled tree probability.
     * In contrast, BDMM-Prime always returns a labeled tree probability.
     * <p>
     * This method exists to convert BDSKY/BDMM test probabilities to be labeled tree probabilities,
     * allowing comparison with BDMM-Prime.
     *
     * @param density BDMM-prime probability density object
     * @return conversion factor
     */
    private double labeledTreeConversionFactor(BirthDeathMigrationDistribution density) {
        Tree tree = (Tree) density.treeInput.get();
        boolean SAmodel = density.parameterizationInput.get().getRemovalProbs()[0][0] != 1.0;
        double factor = -Gamma.logGamma(tree.getLeafNodeCount() + 1);

        if (!SAmodel)
            factor += Math.log(2) * (tree.getLeafNodeCount() - tree.getDirectAncestorNodeCount() - 1);

        return factor;
    }

    /**
     * Basic test for migration rate change
     * Reference from BDMM itself
     * Canonical parameterization
     */
    @ParameterizedTest
    @MethodSource("data")
    public void testLikelihoodMigRateChangeBasicCanonical(Engine method, InitialMatrixStrategy initialStateStrategy, boolean useInverseFlow, boolean parallelize) {

        // Test for uncoloured tree

        String newick = "(t1[&state=0] : 1.5, t2[&state=1] : 0.5);";

        RealScalarParam<NonNegativeReal> originParam = new RealScalarParam<>(2.5, NonNegativeReal.INSTANCE);

        Parameterization parameterization = new CanonicalParameterization();
        parameterization.initByName(
                "typeSet", new TypeSet(2),
                "processLength", originParam,
                "birthRate", new SkylineVectorParameter(
                        null,
                        new RealVectorParam<>(new double[] {2.0}, NonNegativeReal.INSTANCE), 2),
                "deathRate", new SkylineVectorParameter(
                        null,
                        new RealVectorParam<>(new double[] {1.0}, NonNegativeReal.INSTANCE), 2),
                "birthRateAmongDemes", new SkylineMatrixParameter(
                        null,
                        new RealVectorParam<>(new double[] {0.0}, NonNegativeReal.INSTANCE), 2),
                "migrationRate", new SkylineMatrixParameter(
                        new RealVectorParam<>(new double[] {1.0}, NonNegativeReal.INSTANCE),
                        new RealVectorParam<>(new double[] {0.1, 0.2}, NonNegativeReal.INSTANCE), 2),
                "samplingRate", new SkylineVectorParameter(
                        null,
                        new RealVectorParam<>(new double[] {0.5}, NonNegativeReal.INSTANCE), 2),
                "removalProb", new SkylineVectorParameter(
                        null,
                        new RealVectorParam<>(new double[] {1.0}, NonNegativeReal.INSTANCE), 2));


        BirthDeathMigrationDistribution density = new BirthDeathMigrationDistribution();

        density.initByName(
                "method", method,
                "parameterization", parameterization,
                "startTypePriorProbs", new SimplexParam(new double[] {0.5, 0.5}),
                "tree", new TreeParser(newick,
                        false, false,
                        true, 0),
                "conditionOnSurvival", false,
                "typeLabel", "state",
                "initialMatrixStrategy", initialStateStrategy,
                "useInverseFlow", useInverseFlow,
                "parallelize", parallelize
        );

        double logL = density.calculateLogP();

        System.out.println("Birth-death result: " + logL + "\t- Test LikelihoodMigRateChange 1");

        // Reference BDMM (version 0.2.0) 22/06/2017
        assertEquals(-6.7022069383966025 - labeledTreeConversionFactor(density),
                logL, 1e-5);
    }

    /**
     * Basic test for migration rate change
     * Reference from BDMM itself
     * Epi parameterization
     */
    @ParameterizedTest
    @MethodSource("data")
    public void testLikelihoodMigRateChangeBasicEpi(Engine method, InitialMatrixStrategy initialStateStrategy, boolean useInverseFlow, boolean parallelize) {

        // Test for uncoloured tree

        String newick = "(t1[&state=0] : 1.5, t2[&state=1] : 0.5);";

        RealScalarParam<NonNegativeReal> originParam = new RealScalarParam<>(2.5, NonNegativeReal.INSTANCE);

        Parameterization parameterization = new EpiParameterization();
        parameterization.initByName(
                "typeSet", new TypeSet(2),
                "processLength", originParam,
                "Re", new SkylineVectorParameter(
                        null,
                        new RealVectorParam<>(new double[] {4.0 / 3.0, 4.0 / 3.0}, NonNegativeReal.INSTANCE)),
                "becomeUninfectiousRate", new SkylineVectorParameter(
                        null,
                        new RealVectorParam<>(new double[] {1.5, 1.5}, NonNegativeReal.INSTANCE)),
                "ReAmongDemes", new SkylineMatrixParameter(
                        null,
                        new RealVectorParam<>(new double[] {0.0, 0.0}, NonNegativeReal.INSTANCE)),
                "migrationRate", new SkylineMatrixParameter(
                        new RealVectorParam<>(new double[] {1.0}, NonNegativeReal.INSTANCE),
                        new RealVectorParam<>(new double[] {0.1, 0.1, 0.2, 0.2}, NonNegativeReal.INSTANCE)),
                "samplingProportion", new SkylineVectorParameter(
                        null,
                        new RealVectorParam<>(new double[] {1.0 / 3.0, 1.0 / 3.0}, UnitInterval.INSTANCE)),
                "removalProb", new SkylineVectorParameter(
                        null,
                        new RealVectorParam<>(new double[] {1.0, 1.0}, UnitInterval.INSTANCE)));


        BirthDeathMigrationDistribution density = new BirthDeathMigrationDistribution();

        density.initByName(
                "method", method,
                "parameterization", parameterization,
                "startTypePriorProbs", new SimplexParam(new double[] {0.5, 0.5}),
                "tree", new TreeParser(newick,
                        false, false,
                        true, 0),
                "conditionOnSurvival", false,
                "typeLabel", "state",
                "initialMatrixStrategy", initialStateStrategy,
                "useInverseFlow", useInverseFlow,
                "parallelize", parallelize
        );

        double logL = density.calculateLogP();

        System.out.println("Birth-death result: " + logL + "\t- Test LikelihoodMigRateChange 1");

        // Reference BDMM (version 0.2.0) 22/06/2017
        assertEquals(-6.7022069383966025, logL, 1e-5);
    }

    /**
     * Basic test for removal-probability rate change
     * Reference from BDMM itself
     *
     * @throws Exception
     */
    @ParameterizedTest
    @MethodSource("data")
    public void testLikelihoodRemovalProbChangeBasic(Engine method, InitialMatrixStrategy initialStateStrategy, boolean useInverseFlow, boolean parallelize) {

        String newick = "((1[&type=0]: 1.5, 2[&type=0]: 0)3[&type=0]: 3.5, 4[&type=0]: 4) ;";

        Parameterization parameterization = new EpiParameterization();
        parameterization.initByName(
                "typeSet", new TypeSet(1),
                "processLength", new RealScalarParam<>(6.0, NonNegativeReal.INSTANCE),
                "Re", new SkylineVectorParameter(
                        null,
                        new RealVectorParam<>(new double[] {4.0/3.0}, NonNegativeReal.INSTANCE)),
                "becomeUninfectiousRate", new SkylineVectorParameter(
                        null,
                        new RealVectorParam<>(new double[] {1.5}, NonNegativeReal.INSTANCE)),
                "ReAmongDemes", new SkylineMatrixParameter(null, null),
                "migrationRate", new SkylineMatrixParameter(null, null),
                "samplingProportion", new SkylineVectorParameter(
                        null,
                        new RealVectorParam<>(new double[] {1.0/3.0}, UnitInterval.INSTANCE)),
                "removalProb", new SkylineVectorParameter(
                        new RealVectorParam<>(new double[] {1.0}, NonNegativeReal.INSTANCE),
                        new RealVectorParam<>(new double[] {0.3, 0.7}, UnitInterval.INSTANCE)));

        BirthDeathMigrationDistribution density = new BirthDeathMigrationDistribution();
        density.initByName(
                "method", method,
                "parameterization", parameterization,
                "startTypePriorProbs", new SimplexParam(new double[] {1.0}),
                "tree", new TreeParser(newick, false, false, true, 0),
                "conditionOnSurvival", false,
                "typeLabel", "type",
                "initialMatrixStrategy", initialStateStrategy,
                "useInverseFlow", useInverseFlow,
                "parallelize", parallelize
        );

        double logL = density.calculateLogP();

        // Reference BDMM (version 0.2.0) 22/06/2017
        assertEquals(-21.25413884159791 + labeledTreeConversionFactor(density),
                logL, 1e-5);

        BirthDeathMigrationDistribution densityExact = new BirthDeathMigrationDistribution();
        densityExact.initByName(
                "method", method,
                "parameterization", parameterization,
                "startTypePriorProbs", new SimplexParam(new double[] {1.0}),
                "tree", new TreeParser(newick, false, false, true, 0),
                "conditionOnSurvival", false,
                "typeLabel", "type",
                "initialMatrixStrategy", initialStateStrategy,
                "useInverseFlow", useInverseFlow,
                "parallelize", parallelize
        );

        double logLExact = densityExact.calculateLogP();

        // Reference BDMM (version 0.2.0) 22/06/2017
        assertEquals(-21.25413884159791 + labeledTreeConversionFactor(density),
                logLExact, 1e-5);
    }

    /**
     * Direct comparison between numerical and analytical solutions for a tiny example with no rate changes.
     */
    @ParameterizedTest
    @MethodSource("data")
    public void tinyAnalyticalTest(Engine method, InitialMatrixStrategy initialStateStrategy, boolean useInverseFlow, boolean parallelize) {
        String newick = "(1[&type=0]: 1.0, 2[&type=0]: 1.0): 1.0;";

        Parameterization parameterization = new CanonicalParameterization();
        parameterization.initByName(
                "typeSet", new TypeSet(1),
                "processLength", new RealScalarParam<>(2.0, NonNegativeReal.INSTANCE),
                "birthRate", new SkylineVectorParameter(
                        null,
                        new RealVectorParam<>(new double[] {2.0}, NonNegativeReal.INSTANCE)),
                "deathRate", new SkylineVectorParameter(
                        null,
                        new RealVectorParam<>(new double[] {0.5}, NonNegativeReal.INSTANCE)),
                "birthRateAmongDemes", new SkylineMatrixParameter(null, null),
                "migrationRate", new SkylineMatrixParameter(null, null),
                "samplingRate", new SkylineVectorParameter(
                        null,
                        new RealVectorParam<>(new double[] {0.5}, NonNegativeReal.INSTANCE)),
                "rhoSampling", new TimedParameter(
                        new RealVectorParam<>(new double[] {2.0}, NonNegativeReal.INSTANCE),
                        new RealVectorParam<>(new double[] {0.5}, UnitInterval.INSTANCE)),
                "removalProb", new SkylineVectorParameter(
                        null,
                        new RealVectorParam<>(new double[] {0.5}, UnitInterval.INSTANCE)));

        BirthDeathMigrationDistribution density = new BirthDeathMigrationDistribution();
        density.initByName(
                "method", method,
                "parameterization", parameterization,
                "startTypePriorProbs", new SimplexParam(new double[] {1.0}),
                "tree", new TreeParser(newick, false, false, true, 0),
                "conditionOnSurvival", false,
                "typeLabel", "type",
                "initialMatrixStrategy", initialStateStrategy,
                "useInverseFlow", useInverseFlow,
                "parallelize", parallelize
                );

        double logLnumerical = density.calculateLogP();

        BirthDeathMigrationDistribution densityExact = new BirthDeathMigrationDistribution();
        densityExact.initByName(
                "method", method,
                "parameterization", parameterization,
                "startTypePriorProbs", new SimplexParam(new double[] {1.0}),
                "tree", new TreeParser(newick, false, false, true, 0),
                "conditionOnSurvival", false,
                "typeLabel", "type",
                "initialMatrixStrategy", initialStateStrategy,
                "useInverseFlow", useInverseFlow,
                "parallelize", parallelize
                );

        double logLanalytical = densityExact.calculateLogP();

        assertEquals(logLnumerical, logLanalytical, 1e-5);
    }

    /**
     * Two-state test for removal-probability rate change
     * Reference from BDMM itself
     */
    @ParameterizedTest
    @MethodSource("data")
    public void testLikelihoodRemovalProbChangeTwoState(Engine method, InitialMatrixStrategy initialStateStrategy, boolean useInverseFlow, boolean parallelize) {

        String newick = "((1[&type=0]: 1.5, 2[&type=1]: 0)3[&type=0]: 3.5, 4[&type=1]: 4) ;";

        Parameterization parameterization = new EpiParameterization();
        parameterization.initByName(
                "typeSet", new TypeSet(2),
                "processLength", new RealScalarParam<>(6.0, NonNegativeReal.INSTANCE),
                "Re", new SkylineVectorParameter(
                        new RealVectorParam<>(new double[] {1.0}, NonNegativeReal.INSTANCE),
                        new RealVectorParam<>(new double[] {4.0/3.0, 1.1}, NonNegativeReal.INSTANCE),
                        2),
                "becomeUninfectiousRate", new SkylineVectorParameter(
                        new RealVectorParam<>(new double[] {1.0}, NonNegativeReal.INSTANCE),
                        new RealVectorParam<>(new double[] {1.5, 1.4}, NonNegativeReal.INSTANCE),
                        2),
                "ReAmongDemes", new SkylineMatrixParameter(
                        null,
                        new RealVectorParam<>(new double[] {0.0}, NonNegativeReal.INSTANCE),
                        2),
                "migrationRate", new SkylineMatrixParameter(
                        new RealVectorParam<>(new double[] {1.0}, NonNegativeReal.INSTANCE),
                        new RealVectorParam<>(new double[] {0.2, 0.3}, NonNegativeReal.INSTANCE),
                        2),
                "samplingProportion", new SkylineVectorParameter(
                        null,
                        new RealVectorParam<>(new double[] {0.33}, UnitInterval.INSTANCE),
                        2),
                "removalProb", new SkylineVectorParameter(
                        new RealVectorParam<>(new double[] {1.0}, NonNegativeReal.INSTANCE),
                        new RealVectorParam<>(new double[] {0.3, 0.4, 0.7, 0.6}, UnitInterval.INSTANCE)));

        BirthDeathMigrationDistribution density = new BirthDeathMigrationDistribution();
        density.initByName(
                "method", method,
                "parameterization", parameterization,
                "startTypePriorProbs", new SimplexParam(new double[] {0.5, 0.5}),
                "tree", new TreeParser(newick, false, false, true, 0),
                "conditionOnSurvival", false,
                "typeLabel", "type",
                "initialMatrixStrategy", initialStateStrategy,
                "useInverseFlow", useInverseFlow,
                "parallelize", parallelize);

        double logL = density.calculateLogP();

        // Reference BDMM (version 0.2.0) 29/03/2018
        assertEquals(-22.82747259570373, logL, 1e-5);
    }

    /**
     * Basic 1-dim test
     * No rate change, 1 state, no rho-sampling
     * Reference from BDSKY
     */
    @ParameterizedTest
    @MethodSource("data")
    public void testLikelihood1dim(Engine method, InitialMatrixStrategy initialStateStrategy, boolean useInverseFlow, boolean parallelize) {

        Tree tree = new TreeParser( "((3[&state=0] : 1.5, 4[&state=0] : 0.5)[&state=0] : 1 , (1[&state=0] : 2, 2[&state=0] : 1)[&state=0] : 3)[&state=0];",
                false);

        Parameterization parameterization = new EpiParameterization();
        parameterization.initByName(
                "processLength", new RealScalarParam<>(6.0, NonNegativeReal.INSTANCE),
                "typeSet", new TypeSet(1),
                "Re", new SkylineVectorParameter(
                        null,
                        new RealVectorParam<>(new double[] {1.3333333334}, NonNegativeReal.INSTANCE)),
                "becomeUninfectiousRate", new SkylineVectorParameter(
                        null,
                        new RealVectorParam<>(new double[] {1.5}, NonNegativeReal.INSTANCE)),
                "samplingProportion", new SkylineVectorParameter(
                        null,
                        new RealVectorParam<>(new double[] {0.33333333333}, UnitInterval.INSTANCE)),
                "removalProb", new SkylineVectorParameter(
                        null,
                        new RealVectorParam<>(new double[] {1.0}, UnitInterval.INSTANCE)));

        BirthDeathMigrationDistribution density = new BirthDeathMigrationDistribution();
        density.initByName("method", method,
                "parameterization", parameterization,
                "startTypePriorProbs", new SimplexParam(new double[] {1.0}),
                "tree", tree,
                "conditionOnSurvival", false,
                "typeLabel", "state",
                "initialMatrixStrategy", initialStateStrategy,
                "useInverseFlow", useInverseFlow,
                "parallelize", parallelize
                );

        assertEquals(-19.019796073623493 + labeledTreeConversionFactor(density), density.calculateLogP(), 1e-5);   // Reference BDSKY (version 1.3.3)
    }

    /**
     * 1-dim and 1 rate-change test
     * reference from BDSKY
     */
    @ParameterizedTest
    @MethodSource("data")
    public void testLikelihoodRateChange1dim(Engine method, InitialMatrixStrategy initialStateStrategy, boolean useInverseFlow, boolean parallelize) {

        Tree tree = new TreeParser("((3[&state=0] : 1.5, 4[&state=0] : 0.5)[&state=0] : 1 , (1[&state=0] : 2, 2[&state=0] : 1)[&state=0] : 3)[&state=0];",
                false);

        Parameterization parameterization = new EpiParameterization();
        parameterization.initByName(
                "processLength", new RealScalarParam<>(6.0, NonNegativeReal.INSTANCE),
                "typeSet", new TypeSet(1),
                "Re", new SkylineVectorParameter(
                        new RealVectorParam<>(new double[] {3.0}, NonNegativeReal.INSTANCE),
                        new RealVectorParam<>(new double[] {0.6666666667, 1.3333333334}, NonNegativeReal.INSTANCE)),
                "becomeUninfectiousRate", new SkylineVectorParameter(
                        new RealVectorParam<>(new double[] {3.0}, NonNegativeReal.INSTANCE),
                        new RealVectorParam<>(new double[] {4.5, 1.5}, NonNegativeReal.INSTANCE)),
                "samplingProportion", new SkylineVectorParameter(
                        new RealVectorParam<>(new double[] {3.0}, NonNegativeReal.INSTANCE),
                        new RealVectorParam<>(new double[] {0.4444444444, 0.33333333333}, UnitInterval.INSTANCE)),
                "removalProb", new SkylineVectorParameter(
                        null,
                        new RealVectorParam<>(new double[] {1.0}, UnitInterval.INSTANCE)));

        BirthDeathMigrationDistribution density = new BirthDeathMigrationDistribution();
        density.initByName("method", method,
                "parameterization", parameterization,
                "startTypePriorProbs", new SimplexParam(new double[] {1.0}),
                "tree", tree,
                "conditionOnSurvival", false,
                "typeLabel", "state",
                "initialMatrixStrategy", initialStateStrategy,
                "useInverseFlow", useInverseFlow,
                "parallelize", parallelize);

        assertEquals(-33.7573 + labeledTreeConversionFactor(density), density.calculateLogP(), 1e-4); // Reference BDSKY
    }

    /**
     * Basic tests on 2 types situations with migration or birth among demes
     * reference from R
     * @throws Exception
     */
    @ParameterizedTest
    @MethodSource("data")
    public void testLikelihoodCalculationMigTiny(Engine method, InitialMatrixStrategy initialStateStrategy, boolean useInverseFlow, boolean parallelize) throws Exception {

        // migration and no birth among demes

        Tree tree = new TreeParser("(1[&state=0] : 1.5, 2[&state=1] : 0.5)[&state=0];", false);
        Parameterization parameterization = new EpiParameterization();
        parameterization.initByName(
                "processLength", new RealScalarParam<>(2.5, NonNegativeReal.INSTANCE),
                "typeSet", new TypeSet(2),
                "Re", new SkylineVectorParameter(
                        null,
                        new RealVectorParam<>(new double[] {4.0/3.0, 4.0/3.0}, NonNegativeReal.INSTANCE), 2),
                "becomeUninfectiousRate", new SkylineVectorParameter(
                        null,
                        new RealVectorParam<>(new double[] {1.5, 1.5}, NonNegativeReal.INSTANCE), 2),
                "samplingProportion", new SkylineVectorParameter(
                        null,
                        new RealVectorParam<>(new double[] {1.0/3.0, 1.0/3.0}, UnitInterval.INSTANCE), 2),
                "migrationRate", new SkylineMatrixParameter(
                        null,
                        new RealVectorParam<>(new double[] {0.1}, NonNegativeReal.INSTANCE), 2),

                "removalProb", new SkylineVectorParameter(
                        null,
                        new RealVectorParam<>(new double[] {1.0, 1.0}, UnitInterval.INSTANCE), 2));

        BirthDeathMigrationDistribution density = new BirthDeathMigrationDistribution();
        density.initByName("method", method,
                "parameterization", parameterization,
                "startTypePriorProbs", new SimplexParam(new double[] {0.5, 0.5}),
                "tree", tree,
                "conditionOnSurvival", false,
                "typeLabel", "state",
                "initialMatrixStrategy", initialStateStrategy,
                "useInverseFlow", useInverseFlow,
                "parallelize", parallelize
                );

        assertEquals(-7.215222 + labeledTreeConversionFactor(density), density.calculateLogP(), 1e-6); // result from R

        // no migration, symmetric birth among demes

        parameterization.setInputValue("migrationRate", null);
        parameterization.setInputValue("ReAmongDemes", new SkylineMatrixParameter(
                null,
                new RealVectorParam<>(new double[] {0.0666667}, NonNegativeReal.INSTANCE), 2));
        parameterization.initAndValidate();
        density.initAndValidate();

        assertEquals(-7.404888 + labeledTreeConversionFactor(density), density.calculateLogP(), 1e-6); // result from R

        // no migration, asymmetric birth among demes

        parameterization.setInputValue("ReAmongDemes", new SkylineMatrixParameter(
                null,
                new RealVectorParam<>(new double[] {0.0666667, 0.1}, NonNegativeReal.INSTANCE), 2));
        parameterization.initAndValidate();
        density.initAndValidate();

        assertEquals(-7.18723 + labeledTreeConversionFactor(density), density.calculateLogP(), 1e-6); // result from R


        // no migration, asymmetric R0, asymmetric birth among demes

        parameterization.setInputValue("Re", new SkylineVectorParameter(
                null,
                new RealVectorParam<>(new double[] {2.0, 1.3333333}, NonNegativeReal.INSTANCE)));
        parameterization.initAndValidate();
        density.initAndValidate();

        assertEquals(-7.350649 + labeledTreeConversionFactor(density), density.calculateLogP(), 1e-6); // result from R

        // no migration, asymmetric R0, birth among demes, BU rate, samp proportion

        parameterization.setInputValue("Re", new SkylineVectorParameter(
                null,
                new RealVectorParam<>(new double[] {2.0, 1.5}, NonNegativeReal.INSTANCE)));
        parameterization.setInputValue("becomeUninfectiousRate", new SkylineVectorParameter(
                null,
                new RealVectorParam<>(new double[] {2.0, 1.0}, NonNegativeReal.INSTANCE)));
        parameterization.setInputValue("samplingProportion", new SkylineVectorParameter(
                null,
                new RealVectorParam<>(new double[] {0.5, 0.3}, UnitInterval.INSTANCE)));
        parameterization.setInputValue("ReAmongDemes", new SkylineMatrixParameter(
                null,
                new RealVectorParam<>(new double[] {0.1, 0.5}, NonNegativeReal.INSTANCE)));
        parameterization.initAndValidate();
        density.initAndValidate();

        assertEquals(-6.504139 + labeledTreeConversionFactor(density), density.calculateLogP(), 1e-6); // result from R

        // Same params as last test, swapped leaf states

        tree = new TreeParser("(1[&state=1] : 1.5, 2[&state=0] : 0.5)[&state=0];", false);
        density.setInputValue("tree", tree);
        density.initAndValidate();

        assertEquals(-7.700916 + labeledTreeConversionFactor(density), density.calculateLogP(), 1e-6); // result from R
    }

    /**
     * Test migration
     * 2 types, migration, no birth among demes
     * Adapted from BDSKY
     */
    @ParameterizedTest
    @MethodSource("data")
    public void testLikelihoodCalculationMig(Engine method, InitialMatrixStrategy initialStateStrategy, boolean useInverseFlow, boolean parallelize) {

        // uncoloured tree, asymmetric types
        Tree tree = new TreeParser(
                "((3[&type=0] : 1.5, 4[&type=1] : 0.5) : 1 , (1[&type=1] : 2, 2[&type=0] : 1) : 3);",
                false);

        Parameterization parameterization = new EpiParameterization();
        parameterization.initByName(
                "processLength", new RealScalarParam<>(6.0, NonNegativeReal.INSTANCE),
                "typeSet", new TypeSet(2),
                "Re", new SkylineVectorParameter(
                        null,
                        new RealVectorParam<>(new double[] {4.0 / 3.0, 5.0}, NonNegativeReal.INSTANCE)),
                "becomeUninfectiousRate", new SkylineVectorParameter(
                        null,
                        new RealVectorParam<>(new double[] {1.5, 1.25}, NonNegativeReal.INSTANCE)),
                "samplingProportion", new SkylineVectorParameter(
                        null,
                        new RealVectorParam<>(new double[] {1.0 / 3.0, 1.0/2.0}, UnitInterval.INSTANCE)),
                "migrationRate", new SkylineMatrixParameter(
                        null,
                        new RealVectorParam<>(new double[] {0.2, 0.1}, NonNegativeReal.INSTANCE)),

                "removalProb", new SkylineVectorParameter(
                        null,
                        new RealVectorParam<>(new double[] {1.0, 1.0}, UnitInterval.INSTANCE), 2));

        BirthDeathMigrationDistribution density = new BirthDeathMigrationDistribution();
        density.initByName("method", method,
                "parameterization", parameterization,
                "startTypePriorProbs", new SimplexParam(new double[] {0.5, 0.5}),
                "tree", tree,
                "conditionOnSurvival", false,
                "typeLabel", "type",
                "initialMatrixStrategy", initialStateStrategy,
                "useInverseFlow", useInverseFlow,
                "parallelize", parallelize
        );

        assertEquals(-26.53293 + labeledTreeConversionFactor(density), density.calculateLogP(), 1e-4);
    }

    /**
     * Test of migration and infection among demes with rate changes
     * 2 types, no SA
     * Uncoloured tree
     * Reference from BDMM itself (version 0.2.0 28/06/2017)
     * @throws Exception
     */
    @ParameterizedTest
    @MethodSource("data")
    public void testAmongRateChange(Engine method, InitialMatrixStrategy initialStateStrategy, boolean useInverseFlow, boolean parallelize) throws Exception {

        Tree tree = new TreeParser("((3[&type=0]:1.5,4[&type=1]:0.5):1,(1[&type=1]:1,2[&type=0]:1):3);",
                false);

        Parameterization parameterization = new EpiParameterization();
        parameterization.initByName(
                "processLength", new RealScalarParam<>(4.1, NonNegativeReal.INSTANCE),
                "typeSet", new TypeSet(2),
                "Re", new SkylineVectorParameter(
                        new RealVectorParam<>(new double[] {1.0}, NonNegativeReal.INSTANCE),
                        new RealVectorParam<>(new double[] {6.0, 5.0, 2.0, 2.5}, NonNegativeReal.INSTANCE), 2),
                "becomeUninfectiousRate", new SkylineVectorParameter(
                        new RealVectorParam<>(new double[] {1.0}, NonNegativeReal.INSTANCE),
                        new RealVectorParam<>(new double[] {0.5, 0.55, 0.45, 0.6}, NonNegativeReal.INSTANCE), 2),
                "samplingProportion", new SkylineVectorParameter(
                        new RealVectorParam<>(new double[] {1.0}, NonNegativeReal.INSTANCE),
                        new RealVectorParam<>(new double[] {0.5, 0.45, 0.333333, 0.35}, UnitInterval.INSTANCE), 2),
                "ReAmongDemes", new SkylineMatrixParameter(
                        new RealVectorParam<>(new double[] {1.0}, NonNegativeReal.INSTANCE),
                        new RealVectorParam<>(new double[] {1.1, 1.3, 1.2, 1.15}, NonNegativeReal.INSTANCE), 2),
                "migrationRate", new SkylineMatrixParameter(
                        new RealVectorParam<>(new double[] {1.0}, NonNegativeReal.INSTANCE),
                        new RealVectorParam<>(new double[] {0.1, 0.15, 0.2, 0.25}, NonNegativeReal.INSTANCE)),
                "removalProb", new SkylineVectorParameter(
                        null,
                        new RealVectorParam<>(new double[] {1.0, 1.0}, UnitInterval.INSTANCE), 2));

        BirthDeathMigrationDistribution density = new BirthDeathMigrationDistribution();
        density.initByName("method", method,
                "parameterization", parameterization,
                "startTypePriorProbs", new SimplexParam(new double[] {0.5, 0.5}),
                "tree", tree,
                "conditionOnSurvival", false,
                "typeLabel", "type",
                "initialMatrixStrategy", initialStateStrategy,
                "useInverseFlow", useInverseFlow,
                "parallelize", parallelize
        );

        assertEquals(-16.466832439520886 + labeledTreeConversionFactor(density), density.calculateLogP(), 1e-4); // result from BDMM, 28/06/2017
    }

    /**
     * Test of migration and infection among demes with rate changes
     * 2 types, no SA
     * Uncoloured tree
     * Reference from BDMM itself (version 0.2.0 28/06/2017)
     * @throws Exception
     */
    @ParameterizedTest
    @MethodSource("data")
    public void testAmongNoRateChange(Engine method, InitialMatrixStrategy initialStateStrategy, boolean useInverseFlow, boolean parallelize) throws Exception {

        Tree tree = new TreeParser("((3[&type=1]:1.5,4[&type=1]:0.5):1,(1[&type=1]:2,2[&type=1]:1):3);",
                false);

        Parameterization parameterization = new EpiParameterization();
        parameterization.initByName(
                "processLength", new RealScalarParam<>(6.0, NonNegativeReal.INSTANCE),
                "typeSet", new TypeSet(2),
                "Re", new SkylineVectorParameter(
                        null,
                        new RealVectorParam<>(new double[] {0.0, 0.0}, NonNegativeReal.INSTANCE), 2),
                "becomeUninfectiousRate", new SkylineVectorParameter(
                        null,
                        new RealVectorParam<>(new double[] {0.0, 0.75}, NonNegativeReal.INSTANCE), 2),
                "samplingProportion", new SkylineVectorParameter(
                        null,
                        new RealVectorParam<>(new double[] {0.0, 0.7}, UnitInterval.INSTANCE), 2),
                "ReAmongDemes", new SkylineMatrixParameter(
                        null,
                        new RealVectorParam<>(new double[] {0.0, 2.0}, NonNegativeReal.INSTANCE), 2),
                "migrationRate", new SkylineMatrixParameter(
                        null,
                        new RealVectorParam<>(new double[] {0.5, 0.0}, NonNegativeReal.INSTANCE)),
                "removalProb", new SkylineVectorParameter(
                        null,
                        new RealVectorParam<>(new double[] {1.0, 1.0}, UnitInterval.INSTANCE), 2));

        BirthDeathMigrationDistribution density = new BirthDeathMigrationDistribution();
        density.initByName("method", method,
                "parameterization", parameterization,
                "startTypePriorProbs", new SimplexParam(new double[] {1.0, 0.0}),
                "tree", tree,
                "conditionOnSurvival", false,
                "typeLabel", "type",
                "initialMatrixStrategy", initialStateStrategy,
                "useInverseFlow", useInverseFlow,
                "parallelize", parallelize
        );

        assertEquals(-12.1441 + labeledTreeConversionFactor(density), density.calculateLogP(), 1e-4); // tanja's result from R
    }

    /**
     * Test of migration with 3 types
     * No rate change, no SA
     * Reference from BDMM itself
     * @throws Exception
     */
    @ParameterizedTest
    @MethodSource("data")
    public void testMig3types(Engine method, InitialMatrixStrategy initialStateStrategy, boolean useInverseFlow, boolean parallelize) throws Exception {

        Tree tree = new TreeParser("((3[&type=2]:1.5,4[&type=1]:0.5):1,(1[&type=1]:1,2[&type=0]:1):3);",
                false);

        Parameterization parameterization = new EpiParameterization();
        parameterization.initByName(
                "processLength", new RealScalarParam<>(4.1, NonNegativeReal.INSTANCE),
                "typeSet", new TypeSet(3),
                "Re", new SkylineVectorParameter(
                        null,
                        new RealVectorParam<>(new double[] {6.0, 2.0, 5.0}, NonNegativeReal.INSTANCE)),
                "becomeUninfectiousRate", new SkylineVectorParameter(
                        null,
                        new RealVectorParam<>(new double[] {0.5, 0.45, 0.55}, NonNegativeReal.INSTANCE)),
                "samplingProportion", new SkylineVectorParameter(
                        null,
                        new RealVectorParam<>(new double[] {0.5, 0.333333, 0.45}, UnitInterval.INSTANCE)),
                "migrationRate", new SkylineMatrixParameter(
                        null,
                        new RealVectorParam<>(new double[] {0.1, 0.2, 0.15, 0.12, 0.12, 0.15}, NonNegativeReal.INSTANCE)),
                "removalProb", new SkylineVectorParameter(
                        null,
                        new RealVectorParam<>(new double[] {1.0, 1.0, 1.0}, UnitInterval.INSTANCE), 3));

        BirthDeathMigrationDistribution density = new BirthDeathMigrationDistribution();
        density.initByName("method", method,
                "parameterization", parameterization,
                "startTypePriorProbs", new SimplexParam(new double[] {1.0/3.0, 1.0/3.0, 1.0/3.0}),
                "tree", tree,
                "conditionOnSurvival", false,
                "typeLabel", "type",
                "initialMatrixStrategy", initialStateStrategy,
                "useInverseFlow", useInverseFlow,
                "parallelize", parallelize
        );

        assertEquals(-16.88601100061662 + labeledTreeConversionFactor(density), density.calculateLogP(), 1e-4); // result from BDMM, version 0.2.0, 06/07/2017
    }

    /**
     * Likelihood test from the Sasha's SA package.
     */
    @ParameterizedTest
    @MethodSource("data")
    public void testSALikelihoodMini3(Engine method, InitialMatrixStrategy initialStateStrategy, boolean useInverseFlow, boolean parallelize) {
        String newick = "((1:1.0,2:0.0):1.0,3:0):0.0";

        Parameterization parameterization = new CanonicalParameterization();
        parameterization.initByName(
                "processLength", new RealScalarParam<>(10.0, NonNegativeReal.INSTANCE),
                "birthRate", new SkylineVectorParameter(
                        null,
                        new RealVectorParam<>(new double[] {2.0}, NonNegativeReal.INSTANCE)),
                "deathRate", new SkylineVectorParameter(
                        null,
                        new RealVectorParam<>(new double[] {0.99}, NonNegativeReal.INSTANCE)),
                "samplingRate", new SkylineVectorParameter(
                        null,
                        new RealVectorParam<>(new double[] {0.5}, NonNegativeReal.INSTANCE)),
                "removalProb", new SkylineVectorParameter(
                        null,
                        new RealVectorParam<>(new double[] {0.9}, UnitInterval.INSTANCE)));

        BirthDeathMigrationDistribution density = new BirthDeathMigrationDistribution();
        density.initByName(
                "method", method,
                "parameterization", parameterization,
                "tree", new TreeParser(newick, false, false, true,0),
                "conditionOnSurvival", false,
                "parallelize", parallelize
        );

        // this value was calculated by Sasha with Mathematica
        assertEquals(-25.3707 + labeledTreeConversionFactor(density),
                density.calculateLogP(), 1e-5); // likelihood conditioning on at least one sampled individual
    }

    /**
     * 1-dim and 1 rate-change test
     * reference from BDSKY
     */
    @ParameterizedTest
    @MethodSource("data")
    public void testLikelihoodRateChangeCondOnSampling1dim(Engine method, InitialMatrixStrategy initialStateStrategy, boolean useInverseFlow, boolean parallelize) {

        Tree tree = new TreeParser("((3[&state=0] : 1.5, 4[&state=0] : 0.5)[&state=0] : 1 , (1[&state=0] : 2, 2[&state=0] : 1)[&state=0] : 3)[&state=0];",
                false);

        Parameterization parameterization = new EpiParameterization();
        parameterization.initByName(
                "processLength", new RealScalarParam<>(6.0, NonNegativeReal.INSTANCE),
                "typeSet", new TypeSet(1),
                "Re", new SkylineVectorParameter(
                        new RealVectorParam<>(new double[] {3.0}, NonNegativeReal.INSTANCE),
                        new RealVectorParam<>(new double[] {0.6666666667, 1.3333333334}, NonNegativeReal.INSTANCE)),
                "becomeUninfectiousRate", new SkylineVectorParameter(
                        new RealVectorParam<>(new double[] {3.0}, NonNegativeReal.INSTANCE),
                        new RealVectorParam<>(new double[] {4.5, 1.5}, NonNegativeReal.INSTANCE)),
                "samplingProportion", new SkylineVectorParameter(
                        new RealVectorParam<>(new double[] {2.4}, NonNegativeReal.INSTANCE),
                        new RealVectorParam<>(new double[] {0.0, 0.33333333333}, UnitInterval.INSTANCE)),
                "removalProb", new SkylineVectorParameter(
                        null,
                        new RealVectorParam<>(new double[] {1.0}, UnitInterval.INSTANCE)));

        BirthDeathMigrationDistribution density = new BirthDeathMigrationDistribution();
        density.initByName("method", method,
                "parameterization", parameterization,
                "startTypePriorProbs", new SimplexParam(new double[] {1.0}),
                "conditionOnSurvival", true,
                "tree", tree,
                "typeLabel", "state",
                "initialMatrixStrategy", initialStateStrategy,
                "useInverseFlow", useInverseFlow,
                "parallelize", parallelize
        );

        double logPnumeric = density.calculateLogP();
        System.out.println("Numerical solution: " + logPnumeric);

        double logPanalytic = density.calculateLogP();
        System.out.println("Analytical solution: " + logPnumeric);

        assertEquals(logPnumeric, logPanalytic, 1e-4);
    }

    /**
     * Test infection among demes
     * No rate changes
     * Symmetric configuration
     * reference from R
     */
    @ParameterizedTest
    @MethodSource("data")
    public void testLikelihoodCalculationInfAmongDemesSymmetric(Engine method, InitialMatrixStrategy initialStateStrategy, boolean useInverseFlow, boolean parallelize) {

        // uncoloured, symmetric tree

        Tree tree = new TreeParser("((t3[&type=1]:0.004214277605,t4[&type=1]:0.02157681391):0.229186993,(t2[&type=0]:0.624713651,t1[&type=1]:1.347400211):0.06231047755);",
                false);

        Parameterization parameterization = new CanonicalParameterization();
        parameterization.initByName(
                "processLength", new RealScalarParam<>(tree.getRoot().getHeight() + 0.02686563367, NonNegativeReal.INSTANCE),
                "typeSet", new TypeSet(2),
                "birthRate", new SkylineVectorParameter(
                        null,
                        new RealVectorParam<>(new double[] {2.0}, NonNegativeReal.INSTANCE), 2),
                "deathRate", new SkylineVectorParameter(
                        null,
                        new RealVectorParam<>(new double[] {0.5}, NonNegativeReal.INSTANCE), 2),
                "samplingRate", new SkylineVectorParameter(
                        null,
                        new RealVectorParam<>(new double[] {0.5}, NonNegativeReal.INSTANCE), 2),
                "birthRateAmongDemes", new SkylineMatrixParameter(
                        null,
                        new RealVectorParam<>(new double[] {1.0}, NonNegativeReal.INSTANCE), 2),
                "removalProb", new SkylineVectorParameter(
                        null,
                        new RealVectorParam<>(new double[] {1.0, 1.0}, UnitInterval.INSTANCE), 2));

        BirthDeathMigrationDistribution density = new BirthDeathMigrationDistribution();
        density.initByName("method", method,
                "parameterization", parameterization,
                "startTypePriorProbs", new SimplexParam(new double[] {0.5, 0.5}),
                "conditionOnSurvival", true,
                "tree", tree,
                "typeLabel", "type",
                "initialMatrixStrategy", initialStateStrategy,
                "useInverseFlow", useInverseFlow,
                "parallelize", parallelize
        );

        //System.out.println("Log-likelihood " + logL + " - testLikelihoodCalculationInfAmongDemes \t");
        assertEquals(-5.1966118470881 + labeledTreeConversionFactor(density), density.calculateLogP(), 1e-3);

    }

    /**
     * Test infection among demes
     * No rate changes
     * Asymmetric configuration
     * reference from R
     */
    @ParameterizedTest
    @MethodSource("data")
    public void testLikelihoodCalculationInfAmongDemesAsymmetric(Engine method, InitialMatrixStrategy initialStateStrategy, boolean useInverseFlow, boolean parallelize) {

        Tree tree = new TreeParser("((3[&type=1]:1.5,4[&type=0]:0.5):1,(1[&type=0]:2,2[&type=1]:1):3);",
                false);

        Parameterization parameterization = new CanonicalParameterization();
        parameterization.initByName(
                "processLength", new RealScalarParam<>(6.0, NonNegativeReal.INSTANCE),
                "typeSet", new TypeSet(2),
                "birthRate", new SkylineVectorParameter(
                        null,
                        new RealVectorParam<>(new double[] {2.0, 6.25}, NonNegativeReal.INSTANCE), 2),
                "deathRate", new SkylineVectorParameter(
                        null,
                        new RealVectorParam<>(new double[] {1.2, 0.625}, NonNegativeReal.INSTANCE), 2),
                "samplingRate", new SkylineVectorParameter(
                        null,
                        new RealVectorParam<>(new double[] {0.3, 0.625}, NonNegativeReal.INSTANCE), 2),
                "birthRateAmongDemes", new SkylineMatrixParameter(
                        null,
                        new RealVectorParam<>(new double[] {0.2, 0.1}, NonNegativeReal.INSTANCE), 2),
                "removalProb", new SkylineVectorParameter(
                        null,
                        new RealVectorParam<>(new double[] {1.0, 1.0}, UnitInterval.INSTANCE), 2));

        BirthDeathMigrationDistribution density = new BirthDeathMigrationDistribution();
        density.initByName("method", method,
                "parameterization", parameterization,
                "startTypePriorProbs", new SimplexParam(new double[] {0.5, 0.5}),
                "conditionOnSurvival", true,
                "tree", tree,
                "typeLabel", "type",
                "initialMatrixStrategy", initialStateStrategy,
                "useInverseFlow", useInverseFlow,
                "parallelize", parallelize
        );

        assertEquals(-26.7939 + labeledTreeConversionFactor(density), density.calculateLogP(), 1e-5);  //result from R
    }

    /**
     * Basic test on sampled-ancestors lik. calculation.
     * 2 leaves, 1 SA. 1 type, no rho-sampling, no rate-change
     * Reference value from BDSKY (23/03/2017)
     */
    @ParameterizedTest
    @MethodSource("data")
    public void testSALikelihoodMini(Engine method, InitialMatrixStrategy initialStateStrategy, boolean useInverseFlow, boolean parallelize) {

        Tree tree = new TreeParser("((3[&type=0]: 1.5, 6[&type=0]: 0)5[&type=0]: 3.5, 4[&type=0]: 4) ;",
                false);

        Parameterization parameterization = new EpiParameterization();
        parameterization.initByName(
                "processLength", new RealScalarParam<>(6.0, NonNegativeReal.INSTANCE),
                "typeSet", new TypeSet(1),
                "Re", new SkylineVectorParameter(
                        null,
                        new RealVectorParam<>(new double[] {1.5}, NonNegativeReal.INSTANCE)),
                "becomeUninfectiousRate", new SkylineVectorParameter(
                        null,
                        new RealVectorParam<>(new double[] {1.5}, NonNegativeReal.INSTANCE)),
                "samplingProportion", new SkylineVectorParameter(
                        null,
                        new RealVectorParam<>(new double[] {0.2}, UnitInterval.INSTANCE)),
                "removalProb", new SkylineVectorParameter(
                        null,
                        new RealVectorParam<>(new double[] {0.9}, UnitInterval.INSTANCE)));

        BirthDeathMigrationDistribution density = new BirthDeathMigrationDistribution();
        density.initByName("method", method,
                "parameterization", parameterization,
                "startTypePriorProbs", new SimplexParam(new double[] {1.0}),
                "conditionOnSurvival", true,
                "tree", tree,
                "typeLabel", "type",
                "initialMatrixStrategy", initialStateStrategy,
                "useInverseFlow", useInverseFlow,
                "parallelize", parallelize
        );

        assertEquals(-18.854438107814335 + labeledTreeConversionFactor(density), density.calculateLogP(), 1e-4); //Reference value from BDSKY (23/03/2017)
    }

    /**
     * Likelihood test from the Sasha's SA package.
     */
    @ParameterizedTest
    @MethodSource("data")
    public void testSALikelihoodMini2(Engine method, InitialMatrixStrategy initialStateStrategy, boolean useInverseFlow, boolean parallelize) {
        String newick = "((1:1.5,2:0.5):0.5,3:0.0)4:0.0;";

        Parameterization parameterization = new CanonicalParameterization();
        parameterization.initByName(
                "processLength", new RealScalarParam<>(10.0, NonNegativeReal.INSTANCE),
                "birthRate", new SkylineVectorParameter(
                        null,
                        new RealVectorParam<>(new double[] {2.0}, NonNegativeReal.INSTANCE)),
                "deathRate", new SkylineVectorParameter(
                        null,
                        new RealVectorParam<>(new double[] {0.99}, NonNegativeReal.INSTANCE)),
                "samplingRate", new SkylineVectorParameter(
                        null,
                        new RealVectorParam<>(new double[] {0.5}, NonNegativeReal.INSTANCE)),
                "removalProb", new SkylineVectorParameter(
                        null,
                        new RealVectorParam<>(new double[] {0.9}, UnitInterval.INSTANCE)));

        BirthDeathMigrationDistribution density = new BirthDeathMigrationDistribution();
        density.initByName(
                "method", method,
                "parameterization", parameterization,
                "conditionOnSurvival", true,
                "tree", new TreeParser(newick, false, false, true,0),
                "parallelize", parallelize
        );

        // this value was calculated by Sasha with Mathematica
        assertEquals(-22.08332 + labeledTreeConversionFactor(density),
                density.calculateLogP(), 1e-5); // likelihood conditioning on at least one sampled individual
    }

    /**
     * Test on sampled-ancestors lik. calculation with no sampled ancestor
     * No rate-change, one state, 4 tips
     * This state is just there in case something is broken with sampled ancestors,
     * helps for debugging if combined with testSALikelihoodMini for instance
     */
    @ParameterizedTest
    @MethodSource("data")
    public void testSALikelihoodCalculationWithoutAncestors(Engine method, InitialMatrixStrategy initialStateStrategy, boolean useInverseFlow, boolean parallelize) {

        Tree tree = new TreeParser("((3[&type=0] : 1.5, 4[&type=0] : 0.5) : 1 , (1[&type=0] : 2, 2[&type=0] : 1) : 3);",
                false);

        Parameterization parameterization = new EpiParameterization();
        parameterization.initByName(
                "processLength", new ProcessLength(tree),
                "typeSet", new TypeSet(1),
                "Re", new SkylineVectorParameter(
                        null,
                        new RealVectorParam<>(new double[] {1.5}, NonNegativeReal.INSTANCE)),
                "becomeUninfectiousRate", new SkylineVectorParameter(
                        null,
                        new RealVectorParam<>(new double[] {1.5}, NonNegativeReal.INSTANCE)),
                "samplingProportion", new SkylineVectorParameter(
                        null,
                        new RealVectorParam<>(new double[] {0.3}, UnitInterval.INSTANCE)),
                "removalProb", new SkylineVectorParameter(
                        null,
                        new RealVectorParam<>(new double[] {0.9}, UnitInterval.INSTANCE)));

        BirthDeathMigrationDistribution density = new BirthDeathMigrationDistribution();
        density.initByName("method", method,
                "parameterization", parameterization,
                "startTypePriorProbs", new SimplexParam(new double[] {1.0}),
                "conditionOnSurvival", true,
                "conditionOnRoot", true,
                "tree", tree,
                "typeLabel", "type",
                "initialMatrixStrategy", initialStateStrategy,
                "useInverseFlow", useInverseFlow,
                "parallelize", parallelize
        );

        // Conditioned on root:

        assertEquals(-15.545323363405362 + labeledTreeConversionFactor(density), density.calculateLogP(), 1e-4);

        // Conditioned on origin:

        parameterization.setInputValue("processLength", new RealScalarParam<>(10.0, NonNegativeReal.INSTANCE));
        parameterization.initAndValidate();
        density.setInputValue("conditionOnRoot", false);
        density.initAndValidate();

        assertEquals(-25.991511346557598 + labeledTreeConversionFactor(density), density.calculateLogP(), 1e-4);
    }

    /**
     * Tests the case where we have direct ancestors (SA nodes).
     * Tests if two identical trees but with different newick representations have the same likelihood.
     */
    @ParameterizedTest
    @MethodSource("data")
    public void testDirectAncestor(Engine method, InitialMatrixStrategy initialStateStrategy, boolean useInverseFlow, boolean parallelize) {

        // two identical trees up to rotation (the two root children are rotated)
        String newick1 = "((1[&type=0]: 1.5, 2[&type=1]: 0.0)3[&type=0]: 3.5, (4[&type=0]: 1.5, 5[&type=1]: 1.5)6[&type=0]: 3.5) ;";
        String newick2 = "((1[&type=0]: 1.5, 2[&type=1]: 1.5)3[&type=0]: 3.5, (4[&type=0]: 1.5, 5[&type=1]: 0.0)6[&type=0]: 3.5) ;";

        Parameterization parameterization = new EpiParameterization();
        parameterization.initByName(
                "typeSet", new TypeSet(2),
                "processLength", new RealScalarParam<>(6.0, NonNegativeReal.INSTANCE),
                "Re", new SkylineVectorParameter(
                        null,
                        new RealVectorParam<>(new double[] {4.0/3.0, 1.1}, NonNegativeReal.INSTANCE),
                        2),
                "becomeUninfectiousRate", new SkylineVectorParameter(
                        null,
                        new RealVectorParam<>(new double[] {1.5, 1.4}, NonNegativeReal.INSTANCE),
                        2),
                "ReAmongDemes", new SkylineMatrixParameter(
                        null,
                        new RealVectorParam<>(new double[] {0.0}, NonNegativeReal.INSTANCE),
                        2),
                "migrationRate", new SkylineMatrixParameter(
                        null,
                        new RealVectorParam<>(new double[] {0.2, 0.3}, NonNegativeReal.INSTANCE),
                        2),
                "samplingProportion", new SkylineVectorParameter(
                        null,
                        new RealVectorParam<>(new double[] {0.33}, UnitInterval.INSTANCE),
                        2),
                "removalProb", new SkylineVectorParameter(
                        null,
                        new RealVectorParam<>(new double[] {0.3, 0.4}, UnitInterval.INSTANCE)));

        BirthDeathMigrationDistribution density = new BirthDeathMigrationDistribution();
        density.initByName(
                "method", method,
                "parameterization", parameterization,
                "tree", new TreeParser(newick1, false, false, true,0),
                "startTypePriorProbs", new SimplexParam(new double[] {0.5, 0.5}),
                "conditionOnSurvival", false,
                "typeLabel", "type",
                "initialMatrixStrategy", initialStateStrategy,
                "useInverseFlow", useInverseFlow,
                "parallelize", parallelize
        );

        density.setInputValue("tree", new TreeParser(newick1, false, false, true,0));
        density.initAndValidate();
        double logL1 = density.calculateLogP();

        density.setInputValue("tree", new TreeParser(newick2, false, false, true,0));
        density.initAndValidate();
        double logL2 = density.calculateLogP();

        assertEquals(logL1, logL2, 1e-5);
    }

    /**
     * Test simple configuration with one rho-sampling event.
     * One type, no psi-sampling, no sampled-ancestor
     * No rate-changes
     * Reference: R
     */
    @ParameterizedTest
    @MethodSource("data")
    public void testSingleRho(Engine method, InitialMatrixStrategy initialStateStrategy, boolean useInverseFlow, boolean parallelize) {

        Tree tree = new TreeParser("((1[&type=0]: 4.5, 2[&type=0]: 4.5):1,3[&type=0]:5.5);",false);

        Parameterization parameterization = new EpiParameterization();
        parameterization.initByName(
                "typeSet", new TypeSet(1),
                "processLength", new ProcessLength(tree),
                "Re", new SkylineVectorParameter(
                        null,
                        new RealVectorParam<>(new double[] {1.5}, NonNegativeReal.INSTANCE)),
                "becomeUninfectiousRate", new SkylineVectorParameter(
                        null,
                        new RealVectorParam<>(new double[] {1.5}, NonNegativeReal.INSTANCE)),
                "ReAmongDemes", new SkylineMatrixParameter(
                        null,
                        null),
                "migrationRate", new SkylineMatrixParameter(
                        null,
                        null),
                "samplingProportion", new SkylineVectorParameter(
                        null,
                        new RealVectorParam<>(new double[] {0.0}, UnitInterval.INSTANCE)),
                "removalProb", new SkylineVectorParameter(
                        null,
                        new RealVectorParam<>(new double[] {1.0}, UnitInterval.INSTANCE)),
                "rhoSampling", new TimedParameter(
                        new RealVectorParam<>(new double[] {5.5}, NonNegativeReal.INSTANCE),
                        new RealVectorParam<>(new double[] {0.01}, UnitInterval.INSTANCE)));

        BirthDeathMigrationDistribution density = new BirthDeathMigrationDistribution();
        density.initByName(
                "method", method,
                "parameterization", parameterization,
                "startTypePriorProbs", new SimplexParam(new double[] {1.0}),
                "conditionOnRoot", true,
                "tree", tree,
                "typeLabel", "type",
                "initialMatrixStrategy", initialStateStrategy,
                "useInverseFlow", useInverseFlow,
                "parallelize", parallelize
        );

        double logL = density.calculateLogP();

        // this result is from the R package, TreePar:
        // LikConstant(2.25,1.5,0.01,c(4.5,5.5),root=1,survival=1)
        assertEquals(-3.72382 + labeledTreeConversionFactor(density), logL, 1e-4);

        // test with conditioned-on-survival tree
        parameterization.setInputValue("processLength", new RealScalarParam<>(10.0, NonNegativeReal.INSTANCE));
        parameterization.setInputValue("rhoSampling",
                new TimedParameter(new RealVectorParam<>(new double[] {10.0}, NonNegativeReal.INSTANCE),
                        new RealVectorParam<>(new double[] {0.01}, UnitInterval.INSTANCE)));
        parameterization.initAndValidate();

        density.setInputValue("conditionOnSurvival", true);
        density.setInputValue("conditionOnRoot", false);
        density.initAndValidate();

        double logL2 = density.calculateLogP();

        // this result is from R: LikConstant(2.25,1.5,0.01,c(4.5,5.5,5.5+1e-100),root=0,survival=1)
        assertEquals(-7.404227 + labeledTreeConversionFactor(density), logL2, 1e-4);
    }

    /**
     * Basic test for rho-sampling in the past
     * One type, 2 tips, one state
     * No rate changes
     * @throws Exception
     */
    @ParameterizedTest
    @MethodSource("data")
    public void testMultiRho2tips(Engine method, InitialMatrixStrategy initialStateStrategy, boolean useInverseFlow, boolean parallelize) throws Exception {

        // two tips sampled at the same time
        Tree tree = new TreeParser("(3[&type=0]: 4, 4[&type=0]: 4) ;",false);

        RealScalarParam<NonNegativeReal> originParam = new RealScalarParam<>(5.0, NonNegativeReal.INSTANCE);

        Parameterization parameterization = new EpiParameterization();
        parameterization.initByName(
                "typeSet", new TypeSet(1),
                "processLength", originParam,
                "Re", new SkylineVectorParameter(
                        null,
                        new RealVectorParam<>(new double[] {1.5}, NonNegativeReal.INSTANCE)),
                "becomeUninfectiousRate", new SkylineVectorParameter(
                        null,
                        new RealVectorParam<>(new double[] {1.5}, NonNegativeReal.INSTANCE)),
                "ReAmongDemes", new SkylineMatrixParameter(
                        null,
                        null),
                "migrationRate", new SkylineMatrixParameter(
                        null,
                        null),
                "samplingProportion", new SkylineVectorParameter(
                        null,
                        new RealVectorParam<>(new double[] {0.0}, UnitInterval.INSTANCE)),
                "removalProb", new SkylineVectorParameter(
                        null,
                        new RealVectorParam<>(new double[] {1.0}, UnitInterval.INSTANCE)),
                "rhoSampling", new TimedParameter(
                        new RealVectorParam<>(new double[] {0.0, 2.5}, NonNegativeReal.INSTANCE),
                        new RealVectorParam<>(new double[] {1.0, 0.2}, UnitInterval.INSTANCE),
                        originParam));


        BirthDeathMigrationDistribution density = new BirthDeathMigrationDistribution();
        density.initByName(
                "method", method,
                "parameterization", parameterization,
                "startTypePriorProbs", new SimplexParam(new double[] {1.0}),
                "conditionOnSurvival", true,
                "tree", tree,
                "typeLabel", "type",
                "initialMatrixStrategy", initialStateStrategy,
                "useInverseFlow", useInverseFlow,
                "parallelize", parallelize);

        double logL = density.calculateLogP();

        // this result is from BEAST: BDSKY, not double checked in R
        assertEquals(-10.569863754307026, logL, 1e-4);

        // tips sampled at two different times
        tree = new TreeParser("(3[&type=0]: 1.5, 4[&type=0]: 4) ;",false);
        density.setInputValue("tree", tree);
        density.initAndValidate();

        double logL2 = density.calculateLogP();

        // this result is from BEAST: BDSKY, not double checked in R
        assertEquals(-8.099631076932816, logL2, 1e-4);
    }

    /**
     * Test with combined multiple-rho-sampling events in the past and psi-sampling
     * was "TestRhoSasha"
     * One type, no rate-changes, no sampled-ancestors
     * 26 tips
     * @throws Exception
     */
    @ParameterizedTest
    @MethodSource("data")
    public void testMultiRhoSampling(Engine method, InitialMatrixStrategy initialStateStrategy, boolean useInverseFlow, boolean parallelize) throws Exception {
        // Uncoloured tree
        Tree tree = new TreeParser("(((((t1[&type=0]:0.4595008531,t25[&type=0]:0.4595008531)[&type=0]:0.3373053072,t23[&type=0]:0.3567584538)[&type=0]:0.007310819036,t16[&type=0]:0.3489190732)[&type=0]:0.331009529,((t18[&type=0]:0.03315384045,t14[&type=0]:0.03315384045)[&type=0]:0.5063451374,(t10[&type=0]:0.4211543131,t15[&type=0]:0.4211543131)[&type=0]:0.1183446648)[&type=0]:0.5956275305)[&type=0]:0.1158090878,((t19[&type=0]:0.9429393194,((t6[&type=0]:0.363527235,t11[&type=0]:0.4417423167)[&type=0]:0.01881829549,((((t3[&type=0]:0.3071904376,(((t24[&type=0]:0.01065209364,t13[&type=0]:0.01065209364)[&type=0]:0.06076485145,t8[&type=0]:0.07141694509)[&type=0]:0.123620245,(t22[&type=0]:0.1616119808,t2[&type=0]:0.1616119808)[&type=0]:0.03342520927)[&type=0]:0.1121532475)[&type=0]:0.24520579,t9[&type=0]:0.5523962276)[&type=0]:0.3852615426,(((t20[&type=0]:0.2935970782,(t17[&type=0]:0.06569090089,t4[&type=0]:0.06569090089)[&type=0]:0.2279061773)[&type=0]:0.08350780408,(t21[&type=0]:0.05109047139,t5[&type=0]:0.05109047139)[&type=0]:0.3260144109)[&type=0]:0.2298344132,t7[&type=0]:0.6069392955)[&type=0]:0.3307184747)[&type=0]:0.01206284377,t26[&type=0]:0.9497206139)[&type=0]:0.05755333197)[&type=0]:0.03290891884)[&type=0]:0.07263755325,t12[&type=0]:1.112820418)[&type=0]:0.1381151782);",false);

        RealScalarParam<NonNegativeReal> originParam = new RealScalarParam<>(2.0, NonNegativeReal.INSTANCE);

        Parameterization parameterization = new EpiParameterization();
        parameterization.initByName(
                "typeSet", new TypeSet(1),
                "processLength", originParam,
                "Re", new SkylineVectorParameter(
                        null,
                        new RealVectorParam<>(new double[] {3.0/4.5}, NonNegativeReal.INSTANCE)),
                "becomeUninfectiousRate", new SkylineVectorParameter(
                        null,
                        new RealVectorParam<>(new double[] {4.5}, NonNegativeReal.INSTANCE)),
                "ReAmongDemes", new SkylineMatrixParameter(
                        null,
                        null),
                "migrationRate", new SkylineMatrixParameter(
                        null,
                        null),
                "samplingProportion", new SkylineVectorParameter(
                        null,
                        new RealVectorParam<>(new double[] {2.0/4.5}, UnitInterval.INSTANCE)),
                "removalProb", new SkylineVectorParameter(
                        null,
                        new RealVectorParam<>(new double[] {1.0}, UnitInterval.INSTANCE)),
                "rhoSampling", new TimedParameter(
                        new RealVectorParam<>(new double[] {1.0, 1.5, 2.0}, NonNegativeReal.INSTANCE),
                        new RealVectorParam<>(new double[] {0.0, 0.05, 0.01}, UnitInterval.INSTANCE)));


        BirthDeathMigrationDistribution density = new BirthDeathMigrationDistribution();
        density.initByName(
                "method", method,
                "parameterization", parameterization,
                "startTypePriorProbs", new SimplexParam(new double[] {1.0}),
                "conditionOnSurvival", false,
                "tree", tree,
                "typeLabel", "type",
                "initialMatrixStrategy", initialStateStrategy,
                "useInverseFlow", useInverseFlow,
                "parallelize", parallelize
        );

        assertEquals(-124.96086690757612 + labeledTreeConversionFactor(density),
                density.calculateLogP(), 1e-2);     // this result is from BEAST, not double checked in R

        parameterization.setInputValue("rhoSampling",
                new TimedParameter(new RealVectorParam<>(new double[] {0.0, 0.5, 1.0}, NonNegativeReal.INSTANCE),
                        new RealVectorParam<>(new double[] {0.01, 0.05, 0.0}, UnitInterval.INSTANCE),
                        originParam));
        parameterization.initAndValidate();
        density.initAndValidate();

        assertEquals(-124.96086690757612 + labeledTreeConversionFactor(density),
                density.calculateLogP(), 1e-2);     // this result is from BEAST, not double checked in R
    }

    /**
     * Test with multiple cases for rho-sampling in the past combined with rate changes
     * 1 state, no sampled ancestors
     * 26 tips
     */
    @ParameterizedTest
    @MethodSource("data")
    public void testMultiRhoWithRateChanges1(Engine method, InitialMatrixStrategy initialStateStrategy, boolean useInverseFlow, boolean parallelize) {

        Tree tree = new TreeParser("(((((t1[&type=0]:0.4595008531,t25[&type=0]:0.4595008531)[&type=0]:0.3373053072,t23[&type=0]:0.3567584538)[&type=0]:0.007310819036,t16[&type=0]:0.3489190732)[&type=0]:0.331009529,((t18[&type=0]:0.03315384045,t14[&type=0]:0.03315384045)[&type=0]:0.5063451374,(t10[&type=0]:0.4211543131,t15[&type=0]:0.4211543131)[&type=0]:0.1183446648)[&type=0]:0.5956275305)[&type=0]:0.1158090878,((t19[&type=0]:0.9429393194,((t6[&type=0]:0.363527235,t11[&type=0]:0.4417423167)[&type=0]:0.01881829549,((((t3[&type=0]:0.3071904376,(((t24[&type=0]:0.01065209364,t13[&type=0]:0.01065209364)[&type=0]:0.06076485145,t8[&type=0]:0.07141694509)[&type=0]:0.123620245,(t22[&type=0]:0.1616119808,t2[&type=0]:0.1616119808)[&type=0]:0.03342520927)[&type=0]:0.1121532475)[&type=0]:0.24520579,t9[&type=0]:0.5523962276)[&type=0]:0.3852615426,(((t20[&type=0]:0.2935970782,(t17[&type=0]:0.06569090089,t4[&type=0]:0.06569090089)[&type=0]:0.2279061773)[&type=0]:0.08350780408,(t21[&type=0]:0.05109047139,t5[&type=0]:0.05109047139)[&type=0]:0.3260144109)[&type=0]:0.2298344132,t7[&type=0]:0.6069392955)[&type=0]:0.3307184747)[&type=0]:0.01206284377,t26[&type=0]:0.9497206139)[&type=0]:0.05755333197)[&type=0]:0.03290891884)[&type=0]:0.07263755325,t12[&type=0]:1.112820418)[&type=0]:0.1381151782);", false);

        // no rate-change, rho-sampling at present
        RealScalarParam<NonNegativeReal> originParam = new RealScalarParam<>(2.0, NonNegativeReal.INSTANCE);
        Parameterization parameterization = new CanonicalParameterization();
        parameterization.initByName(
                "typeSet", new TypeSet(1),
                "processLength", originParam,
                "birthRate", new SkylineVectorParameter(
                        null,
                        new RealVectorParam<>(new double[] {2.0}, NonNegativeReal.INSTANCE)),
                "deathRate", new SkylineVectorParameter(
                        null,
                        new RealVectorParam<>(new double[] {0.5}, NonNegativeReal.INSTANCE)),
                "samplingRate", new SkylineVectorParameter(
                        null,
                        new RealVectorParam<>(new double[] {0.5}, NonNegativeReal.INSTANCE)),
                "removalProb", new SkylineVectorParameter(
                        null,
                        new RealVectorParam<>(new double[] {1.0}, UnitInterval.INSTANCE)),
                "birthRateAmongDemes", new SkylineMatrixParameter(
                        null,
                        null),
                "migrationRate", new SkylineMatrixParameter(
                        null,
                        null),
                "rhoSampling", new TimedParameter(
                        new RealVectorParam<>(new double[] {0.0}, NonNegativeReal.INSTANCE),
                        new RealVectorParam<>(new double[] {1.0}, UnitInterval.INSTANCE),
                        originParam));

        BirthDeathMigrationDistribution density = new BirthDeathMigrationDistribution();
        density.initByName("method", method,
                "parameterization", parameterization,
                "startTypePriorProbs", new SimplexParam(new double[] {1.0}),
                "conditionOnSurvival", true,
                "tree", tree,
                "typeLabel", "type",
                "initialMatrixStrategy", initialStateStrategy,
                "useInverseFlow", useInverseFlow,
                "parallelize", parallelize
        );

        assertEquals(-21.42666177086957 + labeledTreeConversionFactor(density), density.calculateLogP(), 1e-5);
    }


    @ParameterizedTest
    @MethodSource("data")
    public void testMultiRhoWithRateChanges2(Engine method, InitialMatrixStrategy initialStateStrategy, boolean useInverseFlow, boolean parallelize) {

        Tree tree = new TreeParser("(((((t1[&type=0]:0.4595008531,t25[&type=0]:0.4595008531)[&type=0]:0.3373053072,t23[&type=0]:0.3567584538)[&type=0]:0.007310819036,t16[&type=0]:0.3489190732)[&type=0]:0.331009529,((t18[&type=0]:0.03315384045,t14[&type=0]:0.03315384045)[&type=0]:0.5063451374,(t10[&type=0]:0.4211543131,t15[&type=0]:0.4211543131)[&type=0]:0.1183446648)[&type=0]:0.5956275305)[&type=0]:0.1158090878,((t19[&type=0]:0.9429393194,((t6[&type=0]:0.363527235,t11[&type=0]:0.4417423167)[&type=0]:0.01881829549,((((t3[&type=0]:0.3071904376,(((t24[&type=0]:0.01065209364,t13[&type=0]:0.01065209364)[&type=0]:0.06076485145,t8[&type=0]:0.07141694509)[&type=0]:0.123620245,(t22[&type=0]:0.1616119808,t2[&type=0]:0.1616119808)[&type=0]:0.03342520927)[&type=0]:0.1121532475)[&type=0]:0.24520579,t9[&type=0]:0.5523962276)[&type=0]:0.3852615426,(((t20[&type=0]:0.2935970782,(t17[&type=0]:0.06569090089,t4[&type=0]:0.06569090089)[&type=0]:0.2279061773)[&type=0]:0.08350780408,(t21[&type=0]:0.05109047139,t5[&type=0]:0.05109047139)[&type=0]:0.3260144109)[&type=0]:0.2298344132,t7[&type=0]:0.6069392955)[&type=0]:0.3307184747)[&type=0]:0.01206284377,t26[&type=0]:0.9497206139)[&type=0]:0.05755333197)[&type=0]:0.03290891884)[&type=0]:0.07263755325,t12[&type=0]:1.112820418)[&type=0]:0.1381151782);", false);

        // rate-changes, rho-sampling in the past
        RealScalarParam<NonNegativeReal> originParam = new RealScalarParam<>(2.0, NonNegativeReal.INSTANCE);
        Parameterization parameterization = new EpiParameterization();
        parameterization.initByName(
                "typeSet", new TypeSet(1),
                "processLength", originParam,
                "Re", new SkylineVectorParameter(
                        new RealVectorParam<>(new double[] {1.0, 1.5}, NonNegativeReal.INSTANCE),
                        new RealVectorParam<>(new double[] {3.0/4.5, 2.0/1.5, 4.0/1.5}, NonNegativeReal.INSTANCE)),
                "becomeUninfectiousRate", new SkylineVectorParameter(
                        new RealVectorParam<>(new double[] {1.0, 1.5}, NonNegativeReal.INSTANCE),
                        new RealVectorParam<>(new double[] {4.5, 1.5, 1.5}, NonNegativeReal.INSTANCE)),
                "samplingProportion", new SkylineVectorParameter(
                        new RealVectorParam<>(new double[] {1.0, 1.5}, NonNegativeReal.INSTANCE),
                        new RealVectorParam<>(new double[] {2.0/4.5, 0.5/1.5, 1.0/1.5}, UnitInterval.INSTANCE)),
                "removalProb", new SkylineVectorParameter(
                        null,
                        new RealVectorParam<>(new double[] {1.0}, UnitInterval.INSTANCE)),
                "ReAmongDemes", new SkylineMatrixParameter(
                        null,
                        null),
                "migrationRate", new SkylineMatrixParameter(
                        null,
                        null),
                "rhoSampling", new TimedParameter(
                        new RealVectorParam<>(new double[] {2.0}, NonNegativeReal.INSTANCE),
                        new RealVectorParam<>(new double[] {0.01}, UnitInterval.INSTANCE)));

        BirthDeathMigrationDistribution density = new BirthDeathMigrationDistribution();
        density.initByName("method", method,
                "parameterization", parameterization,
                "startTypePriorProbs", new SimplexParam(new double[] {1.0}),
                "conditionOnSurvival", false,
                "tree", tree,
                "typeLabel", "type",
                "initialMatrixStrategy", initialStateStrategy,
                "useInverseFlow", useInverseFlow,
                "parallelize", parallelize
        );

        assertEquals(-87.59718586549747 + labeledTreeConversionFactor(density), density.calculateLogP(), 1e-4);
    }

    @ParameterizedTest
    @MethodSource("data")
    public void testMultiRhoWithRateChanges3(Engine method, InitialMatrixStrategy initialStateStrategy, boolean useInverseFlow, boolean parallelize) {

        Tree tree = new TreeParser("(((((t1[&type=0]:0.4595008531,t25[&type=0]:0.4595008531)[&type=0]:0.3373053072,t23[&type=0]:0.3567584538)[&type=0]:0.007310819036,t16[&type=0]:0.3489190732)[&type=0]:0.331009529,((t18[&type=0]:0.03315384045,t14[&type=0]:0.03315384045)[&type=0]:0.5063451374,(t10[&type=0]:0.4211543131,t15[&type=0]:0.4211543131)[&type=0]:0.1183446648)[&type=0]:0.5956275305)[&type=0]:0.1158090878,((t19[&type=0]:0.9429393194,((t6[&type=0]:0.363527235,t11[&type=0]:0.4417423167)[&type=0]:0.01881829549,((((t3[&type=0]:0.3071904376,(((t24[&type=0]:0.01065209364,t13[&type=0]:0.01065209364)[&type=0]:0.06076485145,t8[&type=0]:0.07141694509)[&type=0]:0.123620245,(t22[&type=0]:0.1616119808,t2[&type=0]:0.1616119808)[&type=0]:0.03342520927)[&type=0]:0.1121532475)[&type=0]:0.24520579,t9[&type=0]:0.5523962276)[&type=0]:0.3852615426,(((t20[&type=0]:0.2935970782,(t17[&type=0]:0.06569090089,t4[&type=0]:0.06569090089)[&type=0]:0.2279061773)[&type=0]:0.08350780408,(t21[&type=0]:0.05109047139,t5[&type=0]:0.05109047139)[&type=0]:0.3260144109)[&type=0]:0.2298344132,t7[&type=0]:0.6069392955)[&type=0]:0.3307184747)[&type=0]:0.01206284377,t26[&type=0]:0.9497206139)[&type=0]:0.05755333197)[&type=0]:0.03290891884)[&type=0]:0.07263755325,t12[&type=0]:1.112820418)[&type=0]:0.1381151782);", false);

        // rate-changes, rho-sampling in the past and present
        RealScalarParam<NonNegativeReal> originParam = new RealScalarParam<>(2.0, NonNegativeReal.INSTANCE);
        Parameterization parameterization = new EpiParameterization();
        parameterization.initByName(
                "typeSet", new TypeSet(1),
                "processLength", originParam,
                "Re", new SkylineVectorParameter(
                        new RealVectorParam<>(new double[] {1.0, 1.5}, NonNegativeReal.INSTANCE),
                        new RealVectorParam<>(new double[] {3.0/4.5, 2.0/1.5, 4.0/1.5}, NonNegativeReal.INSTANCE)),
                "becomeUninfectiousRate", new SkylineVectorParameter(
                        new RealVectorParam<>(new double[] {1.0, 1.5}, NonNegativeReal.INSTANCE),
                        new RealVectorParam<>(new double[] {4.5, 1.5, 1.5}, NonNegativeReal.INSTANCE)),
                "samplingProportion", new SkylineVectorParameter(
                        new RealVectorParam<>(new double[] {1.0, 1.5}, NonNegativeReal.INSTANCE),
                        new RealVectorParam<>(new double[] {2.0/4.5, 0.5/1.5, 1.0/1.5}, UnitInterval.INSTANCE)),
                "removalProb", new SkylineVectorParameter(
                        null,
                        new RealVectorParam<>(new double[] {1.0}, UnitInterval.INSTANCE)),
                "ReAmongDemes", new SkylineMatrixParameter(
                        null,
                        null),
                "migrationRate", new SkylineMatrixParameter(
                        null,
                        null),
                "rhoSampling", new TimedParameter(
                        new RealVectorParam<>(new double[] {1.0, 2.0}, NonNegativeReal.INSTANCE),
                        new RealVectorParam<>(new double[] {0.05, 0.01}, UnitInterval.INSTANCE)));

        BirthDeathMigrationDistribution density = new BirthDeathMigrationDistribution();
        density.initByName("method", method,
                "parameterization", parameterization,
                "startTypePriorProbs", new SimplexParam(new double[] {1.0}),
                "conditionOnSurvival", false,
                "tree", tree,
                "typeLabel", "type",
                "initialMatrixStrategy", initialStateStrategy,
                "useInverseFlow", useInverseFlow,
                "parallelize", parallelize
        );

        assertEquals(-87.96488 + labeledTreeConversionFactor(density), density.calculateLogP(), 1e-1);
    }

    @ParameterizedTest
    @MethodSource("data")
    public void testMultiRhoWithRateChanges4(Engine method, InitialMatrixStrategy initialStateStrategy, boolean useInverseFlow, boolean parallelize) {

        Tree tree = new TreeParser("(((((t1[&type=0]:0.4595008531,t25[&type=0]:0.4595008531)[&type=0]:0.3373053072,t23[&type=0]:0.3567584538)[&type=0]:0.007310819036,t16[&type=0]:0.3489190732)[&type=0]:0.331009529,((t18[&type=0]:0.03315384045,t14[&type=0]:0.03315384045)[&type=0]:0.5063451374,(t10[&type=0]:0.4211543131,t15[&type=0]:0.4211543131)[&type=0]:0.1183446648)[&type=0]:0.5956275305)[&type=0]:0.1158090878,((t19[&type=0]:0.9429393194,((t6[&type=0]:0.363527235,t11[&type=0]:0.4417423167)[&type=0]:0.01881829549,((((t3[&type=0]:0.3071904376,(((t24[&type=0]:0.01065209364,t13[&type=0]:0.01065209364)[&type=0]:0.06076485145,t8[&type=0]:0.07141694509)[&type=0]:0.123620245,(t22[&type=0]:0.1616119808,t2[&type=0]:0.1616119808)[&type=0]:0.03342520927)[&type=0]:0.1121532475)[&type=0]:0.24520579,t9[&type=0]:0.5523962276)[&type=0]:0.3852615426,(((t20[&type=0]:0.2935970782,(t17[&type=0]:0.06569090089,t4[&type=0]:0.06569090089)[&type=0]:0.2279061773)[&type=0]:0.08350780408,(t21[&type=0]:0.05109047139,t5[&type=0]:0.05109047139)[&type=0]:0.3260144109)[&type=0]:0.2298344132,t7[&type=0]:0.6069392955)[&type=0]:0.3307184747)[&type=0]:0.01206284377,t26[&type=0]:0.9497206139)[&type=0]:0.05755333197)[&type=0]:0.03290891884)[&type=0]:0.07263755325,t12[&type=0]:1.112820418)[&type=0]:0.1381151782);", false);

        // rate-changes, rho-sampling in the past and present, with reversed times
        Parameterization parameterization = new EpiParameterization();
        parameterization.initByName(
                "processLength", new ProcessLength(tree),
                "typeSet", new TypeSet(1),
                "Re", new SkylineVectorParameter(
                        new RealVectorParam<>(new double[] {0.5, 1.0, 1.1}, NonNegativeReal.INSTANCE),
                        new RealVectorParam<>(new double[] {3.0/4.5, 2.0/1.5, 4.0/1.5, 4.0/2.5}, NonNegativeReal.INSTANCE)),
                "becomeUninfectiousRate", new SkylineVectorParameter(
                        new RealVectorParam<>(new double[] {0.5, 1.0, 1.1}, NonNegativeReal.INSTANCE),
                        new RealVectorParam<>(new double[] {4.5, 1.5, 1.5, 2.5}, NonNegativeReal.INSTANCE)),
                "samplingProportion", new SkylineVectorParameter(
                        new RealVectorParam<>(new double[] {0.5, 1.0, 1.1}, NonNegativeReal.INSTANCE),
                        new RealVectorParam<>(new double[] {2.0/4.5, 0.5/1.5, 1.0/1.5, 2.0/2.5}, UnitInterval.INSTANCE)),
                "removalProb", new SkylineVectorParameter(
                        null,
                        new RealVectorParam<>(new double[] {1.0}, UnitInterval.INSTANCE)),
                "ReAmongDemes", new SkylineMatrixParameter(
                        null,
                        null),
                "migrationRate", new SkylineMatrixParameter(
                        null,
                        null),
                "rhoSampling", new TimedParameter(
                        new RealVectorParam<>(new double[] {1.0, tree.getRoot().getHeight()}, NonNegativeReal.INSTANCE),
                        new RealVectorParam<>(new double[] {0.05, 0.01}, UnitInterval.INSTANCE)));

        BirthDeathMigrationDistribution density = new BirthDeathMigrationDistribution();
        density.initByName("method", method,
                "parameterization", parameterization,
                "startTypePriorProbs", new SimplexParam(new double[] {1.0}),
                "conditionOnSurvival", false,
                "conditionOnRoot", true,
                "tree", tree,
                "typeLabel", "type",
                "initialMatrixStrategy", initialStateStrategy,
                "useInverseFlow", useInverseFlow,
                "parallelize", parallelize
        );

        assertEquals(-99.0428845398644 + labeledTreeConversionFactor(density), density.calculateLogP(), 1e-1);
    }

    /**
     * Test on combining migration with rho-sampling
     * Reference from BDMM
     */
    @ParameterizedTest
    @MethodSource("data")
    public void testLikelihoodMigrationRhoSampling(Engine method, InitialMatrixStrategy initialStateStrategy, boolean useInverseFlow, boolean parallelize) {

        Tree tree = new TreeParser("((1[&type=0]: 4.5, 2[&type=1]: 4.5):1,3[&type=0]:5.5);",
                false);

        Parameterization parameterization = new EpiParameterization();
        parameterization.initByName(
                "processLength", new ProcessLength(tree),
                "typeSet", new TypeSet(2),
                "Re", new SkylineVectorParameter(
                        null,
                        new RealVectorParam<>(new double[] {1.5, 1.4}, NonNegativeReal.INSTANCE)),
                "becomeUninfectiousRate", new SkylineVectorParameter(
                        null,
                        new RealVectorParam<>(new double[] {1.5, 1.3}, NonNegativeReal.INSTANCE)),
                "samplingProportion", new SkylineVectorParameter(
                        null,
                        new RealVectorParam<>(new double[] {0.0}, NonNegativeReal.INSTANCE), 2),
                "removalProb", new SkylineVectorParameter(
                        null,
                        new RealVectorParam<>(new double[] {1.0}, NonNegativeReal.INSTANCE), 2),
                "migrationRate", new SkylineMatrixParameter(
                        null,
                        new RealVectorParam<>(new double[] {0.3, 0.4}, NonNegativeReal.INSTANCE)),
                "rhoSampling", new TimedParameter(
                        new RealVectorParam<>(new double[] {0.0}, NonNegativeReal.INSTANCE),
                        new RealVectorParam<>(new double[] {0.01, 0.015}, UnitInterval.INSTANCE),
                        new ProcessLength(tree)));

        BirthDeathMigrationDistribution density = new BirthDeathMigrationDistribution();
        density.initByName("method", method,
                "parameterization", parameterization,
                "startTypePriorProbs", new SimplexParam(new double[] {0.6, 0.4}),
                "conditionOnSurvival", false,
                "conditionOnRoot", true,
                "tree", tree,
                "typeLabel", "type",
                "initialMatrixStrategy", initialStateStrategy,
                "useInverseFlow", useInverseFlow,
                "parallelize", parallelize
        );

        // Corrected value from BDMM (original was incorrectly conditioned)
        assertEquals(-5.5751511486962215 + labeledTreeConversionFactor(density), density.calculateLogP(), 1e-4);
    }

    /**
     * Basic test on sampled-ancestors lik. calculation with multi-rho sampling.
     * 2 tips, 1 SA. 1 type, no rate-change
     * Coloured and uncoloured trees
     * Reference value from BDSKY (06/04/2017)
     * @throws Exception
     */
    @ParameterizedTest
    @MethodSource("data")
    public void testSALikelihoodMultiRho(Engine method, InitialMatrixStrategy initialStateStrategy, boolean useInverseFlow, boolean parallelize) throws Exception {

        Tree tree = new TreeParser("((3[&type=0]: 1.5, 6[&type=0]: 0)5[&type=0]: 3.5, 4[&type=0]: 4) ;",false);

        RealScalarParam<NonNegativeReal> origin = new RealScalarParam<>(6.0, NonNegativeReal.INSTANCE);

        Parameterization parameterization = new EpiParameterization();
        parameterization.initByName(
                "processLength", origin,
                "typeSet", new TypeSet(1),
                "Re", new SkylineVectorParameter(
                        null,
                        new RealVectorParam<>(new double[] {1.5}, NonNegativeReal.INSTANCE)),
                "becomeUninfectiousRate", new SkylineVectorParameter(
                        null,
                        new RealVectorParam<>(new double[] {1.5}, NonNegativeReal.INSTANCE)),
                "samplingProportion", new SkylineVectorParameter(
                        null,
                        new RealVectorParam<>(new double[] {0.2}, UnitInterval.INSTANCE)),
                "removalProb", new SkylineVectorParameter(
                        null,
                        new RealVectorParam<>(new double[] {0.9}, UnitInterval.INSTANCE)),
                "rhoSampling", new TimedParameter(
                        new RealVectorParam<>(new double[] {0.0, 1.5}, NonNegativeReal.INSTANCE),
                        new RealVectorParam<>(new double[] {0.05, 0.3}, UnitInterval.INSTANCE),
                        origin));

        BirthDeathMigrationDistribution density = new BirthDeathMigrationDistribution();
        density.initByName("method", method,
                "parameterization", parameterization,
                "startTypePriorProbs", new SimplexParam(new double[] {1.0}),
                "conditionOnSurvival", true,
                "tree", tree,
                "typeLabel", "type",
                "initialMatrixStrategy", initialStateStrategy,
                "useInverseFlow", useInverseFlow,
                "parallelize", parallelize
        );

        BirthDeathMigrationDistribution primeDensity = new BirthDeathMigrationDistribution();
        primeDensity.initByName("method", method,
                "parameterization", parameterization,
                "startTypePriorProbs", new SimplexParam(new double[] {1.0}),
                "conditionOnSurvival", true,
                "tree", tree,
                "typeLabel", "type",
                "initialMatrixStrategy", initialStateStrategy,
                "useInverseFlow", useInverseFlow,
                "parallelize", parallelize
        );

        // assertEquals(-22.348462265673483 + labeledTreeConversionFactor(density), density.calculateLogP(), 1e-5); //Reference value from BDSKY (06/04/2017)
        assertEquals(density.calculateLogP(), primeDensity.calculateLogP(), 1e-5); //Reference value from BDSKY (06/04/2017)
    }
}
