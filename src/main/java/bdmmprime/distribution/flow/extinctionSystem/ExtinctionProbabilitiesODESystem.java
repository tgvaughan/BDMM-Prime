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

package bdmmprime.distribution.flow.extinctionSystem;

import bdmmprime.distribution.flow.intervals.Interval;
import bdmmprime.distribution.flow.intervals.IntervalODESystem;
import bdmmprime.parameterization.Parameterization;
import org.apache.commons.math3.ode.ContinuousOutputModel;

import java.util.List;

/**
 * This class represents the ODE system for extinction probabilities. Its state consists of two blocks of
 * length n: the extinction probabilities p0 followed by the survival probabilities q = 1 - p0. Both blocks are
 * integrated directly, so that the relative tolerance applies to p0 when p0 is small and to q when q is small.
 */
public class ExtinctionProbabilitiesODESystem extends IntervalODESystem {

    private final double[][] birthRates;
    private final double[][] deathRates;
    private final double[][] samplingRates;
    private final double[][][] crossBirthRates;
    private final double[][][] migrationRates;

    private int currentParameterizationInterval;

    public ExtinctionProbabilitiesODESystem(Parameterization parameterization, List<Interval> intervals, double absoluteTolerance, double relativeTolerance) {
        super(parameterization, intervals, absoluteTolerance, relativeTolerance);

        this.birthRates = this.parameterization.getBirthRates();
        this.deathRates = this.parameterization.getDeathRates();
        this.samplingRates = this.parameterization.getSamplingRates();
        this.crossBirthRates = this.parameterization.getCrossBirthRates();
        this.migrationRates = this.parameterization.getMigRates();
    }

    @Override
    public int getDimension() {
        return 2 * this.parameterization.getNTypes();
    }

    /**
     * Integrates the system along the given interval using the rates of its parameterization interval.
     */
    @Override
    protected ContinuousOutputModel integrate(double[] initialState, double start, double end, Interval interval) {
        this.currentParameterizationInterval = interval.parameterizationInterval();
        return super.integrate(initialState, start, end, interval);
    }

    @Override
    public void computeDerivatives(double t, double[] y, double[] yDot) {
        if (Double.isNaN(t)) {
            throw new IllegalStateException("NaN detected during integration.");
        }

        int interval = this.currentParameterizationInterval;
        int nTypes = this.parameterization.getNTypes();

        for (int i = 0; i < nTypes; i++) {
            double birthRate = this.birthRates[interval][i];
            double deathRate = this.deathRates[interval][i];
            double samplingRate = this.samplingRates[interval][i];

            double extinctionProbability = y[i];
            double survivalProbability = y[nTypes + i];

            // derivative of p0 with the base terms written in p0
            yDot[i] = (birthRate + deathRate + samplingRate) * extinctionProbability
                    - birthRate * extinctionProbability * extinctionProbability
                    - deathRate;

            // derivative of q = 1 - p0 with the base terms written in q
            yDot[nTypes + i] = birthRate * survivalProbability * survivalProbability
                    - (birthRate - deathRate - samplingRate) * survivalProbability
                    - samplingRate;

            for (int j = 0; j < nTypes; j++) {
                if (i == j) {
                    continue;
                }

                // the cross birth term is p0_i * q_j in both blocks, which avoids 1 - p0_j in the p0 block
                double crossBirthFlux = this.crossBirthRates[interval][i][j] * extinctionProbability * y[nTypes + j];
                double migrationRate = this.migrationRates[interval][i][j];

                yDot[i] += crossBirthFlux + migrationRate * (extinctionProbability - y[j]);
                yDot[nTypes + i] -= crossBirthFlux + migrationRate * (y[nTypes + j] - survivalProbability);
            }
        }
    }

    @Override
    protected void handleParameterizationIntervalBoundary(double boundaryTime, int oldInterval, int newInterval, double[] state) {
        super.handleParameterizationIntervalBoundary(boundaryTime, oldInterval, newInterval, state);

        // include rho sampling effects

        int nTypes = this.parameterization.getNTypes();

        for (int type = 0; type < nTypes; type++) {
            double rho = this.parameterization.getRhoValues()[newInterval][type];
            double extinctionProbability = state[type];

            state[type] = extinctionProbability * (1.0 - rho);
            state[nTypes + type] += rho * extinctionProbability;
        }
    }

}
