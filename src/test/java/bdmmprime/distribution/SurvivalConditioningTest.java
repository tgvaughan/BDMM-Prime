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

import bdmmprime.parameterization.CanonicalParameterization;
import bdmmprime.parameterization.Parameterization;
import bdmmprime.parameterization.SkylineMatrixParameter;
import bdmmprime.parameterization.SkylineVectorParameter;
import bdmmprime.parameterization.TypeSet;
import beast.base.evolution.tree.TreeParser;
import beast.base.spec.domain.NonNegativeReal;
import beast.base.spec.inference.parameter.RealScalarParam;
import beast.base.spec.inference.parameter.RealVectorParam;
import beast.base.spec.inference.parameter.SimplexParam;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.Arguments;
import org.junit.jupiter.params.provider.MethodSource;

import java.util.ArrayList;
import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;

/**
 * Tests the conditioning on survival when the survival probability at the origin is tiny. This happens when an
 * early interval without sampling is strongly subcritical, so that the extinction probability p0 at the origin is
 * close to 1 and 1 - p0 is prone to cancellation.
 * <p>
 * With a single type, the condition density is the survival probability q0 = 1 - p0 at the origin, so the
 * conditioned and unconditioned log-likelihoods differ by exactly -log(q0). For piecewise-constant rates, q0 has
 * a closed form that serves as the reference.
 */
public class SurvivalConditioningTest {

    private static final double PROCESS_LENGTH = 10.0;
    private static final double CHANGE_TIME = 5.0;

    private static final double EARLY_BIRTH_RATE = 1.0;
    private static final double RECENT_BIRTH_RATE = 2.0;
    private static final double RECENT_DEATH_RATE = 1.0;
    private static final double RECENT_SAMPLING_RATE = 0.5;

    private static final String NEWICK = "((t1[&type=0]:1.0,t2[&type=0]:0.5):1.5,(t3[&type=0]:2.0,t4[&type=0]:1.2):0.8);";

    /**
     * Creates the test cases. Format is { engine, death rate of the early interval }. Larger death rates in the
     * early interval give smaller survival probabilities at the origin, from about 1e-5 down to about 1e-30.
     */
    public static List<Arguments> getTestCases() {
        List<Arguments> testCases = new ArrayList<>();

        for (Engine engine : new Engine[] {Engine.classic, Engine.flow}) {
            for (double earlyDeathRate : new double[] {3.0, 6.0, 9.0, 12.0, 15.0}) {
                testCases.add(Arguments.of(engine, earlyDeathRate));
            }
        }

        return testCases;
    }

    @ParameterizedTest
    @MethodSource("getTestCases")
    public void testSurvivalConditioningWithTinySurvivalProbability(Engine engine, double earlyDeathRate) {
        double logLikelihoodConditioned = this.computeLogLikelihood(engine, earlyDeathRate, true);
        double logLikelihoodUnconditioned = this.computeLogLikelihood(engine, earlyDeathRate, false);

        double logSurvivalProbability = logLikelihoodUnconditioned - logLikelihoodConditioned;
        double expectedLogSurvivalProbability = Math.log(this.computeSurvivalProbabilityAtOrigin(earlyDeathRate));

        System.out.printf("engine=%s earlyDeathRate=%.1f expected log(q0)=%.6f computed log(q0)=%.6f%n",
                engine, earlyDeathRate, expectedLogSurvivalProbability, logSurvivalProbability);

        assertEquals(expectedLogSurvivalProbability, logSurvivalProbability, 1e-4);
    }

    /**
     * Computes the log-likelihood of the test tree under the two-interval single-type model.
     */
    private double computeLogLikelihood(Engine engine, double earlyDeathRate, boolean conditionOnSurvival) {
        RealScalarParam<NonNegativeReal> processLength = new RealScalarParam<>(PROCESS_LENGTH, NonNegativeReal.INSTANCE);

        Parameterization parameterization = new CanonicalParameterization();
        parameterization.initByName(
                "typeSet", new TypeSet(1),
                "processLength", processLength,
                "birthRate", this.createSkylineVector(EARLY_BIRTH_RATE, RECENT_BIRTH_RATE),
                "deathRate", this.createSkylineVector(earlyDeathRate, RECENT_DEATH_RATE),
                "samplingRate", this.createSkylineVector(0.0, RECENT_SAMPLING_RATE),
                "removalProb", new SkylineVectorParameter(
                        null,
                        new RealVectorParam<>(new double[] {1.0}, NonNegativeReal.INSTANCE), 1),
                "birthRateAmongDemes", new SkylineMatrixParameter(
                        null,
                        new RealVectorParam<>(new double[] {0.0}, NonNegativeReal.INSTANCE), 1),
                "migrationRate", new SkylineMatrixParameter(
                        null,
                        new RealVectorParam<>(new double[] {0.0}, NonNegativeReal.INSTANCE), 1));

        BirthDeathMigrationDistribution density = new BirthDeathMigrationDistribution();
        density.initByName(
                "engine", engine,
                "parameterization", parameterization,
                "startTypePriorProbs", new SimplexParam(new double[] {1.0}),
                "conditionOnSurvival", conditionOnSurvival,
                "tree", new TreeParser(NEWICK, false, false, true, 0),
                "typeLabel", "type",
                "useAnalyticalSingleTypeSolution", false,
                "parallelize", false
        );

        return density.calculateLogP();
    }

    /**
     * Creates a single-type skyline parameter with one value before and one value after the change time.
     */
    private SkylineVectorParameter createSkylineVector(double earlyValue, double recentValue) {
        return new SkylineVectorParameter(
                new RealVectorParam<>(new double[] {CHANGE_TIME}, NonNegativeReal.INSTANCE),
                new RealVectorParam<>(new double[] {earlyValue, recentValue}, NonNegativeReal.INSTANCE),
                1);
    }

    /**
     * Computes the survival probability at the origin in closed form, integrating from the present back through
     * the recent and then the early interval.
     */
    private double computeSurvivalProbabilityAtOrigin(double earlyDeathRate) {
        double survivalAtChangeTime = this.computeSurvivalProbability(
                0.0, PROCESS_LENGTH - CHANGE_TIME,
                RECENT_BIRTH_RATE, RECENT_DEATH_RATE, RECENT_SAMPLING_RATE
        );

        return this.computeSurvivalProbability(
                survivalAtChangeTime, CHANGE_TIME,
                EARLY_BIRTH_RATE, earlyDeathRate, 0.0
        );
    }

    /**
     * Solves the survival probability equation dq/dtau = s + (b - d - s) q - b q^2 in backward time tau for
     * constant rates, starting from the initial survival probability and running for the given duration.
     * <p>
     * The right-hand side factorizes as -b (q - qPlus) (q - qMinus) with roots qPlus >= 0 >= qMinus. The solution
     * satisfies (q - qPlus) / (q - qMinus) = c exp(-a tau), with a = b (qPlus - qMinus). Each root is computed with
     * the formula that avoids cancellation, so the result keeps its relative accuracy even when q is tiny.
     */
    private double computeSurvivalProbability(double initialSurvival, double duration, double birthRate,
                                              double deathRate, double samplingRate) {
        double netGrowthRate = birthRate - deathRate - samplingRate;
        double a = Math.sqrt(netGrowthRate * netGrowthRate + 4 * birthRate * samplingRate);

        double qPlus;
        double qMinus;
        if (netGrowthRate >= 0) {
            qPlus = (netGrowthRate + a) / (2 * birthRate);
            qMinus = -samplingRate / (birthRate * qPlus);
        } else {
            qMinus = (netGrowthRate - a) / (2 * birthRate);
            qPlus = -samplingRate / (birthRate * qMinus);
        }

        double c = (initialSurvival - qPlus) / (initialSurvival - qMinus);
        double decay = c * Math.exp(-a * duration);

        return (qPlus - qMinus * decay) / (1 - decay);
    }

}
