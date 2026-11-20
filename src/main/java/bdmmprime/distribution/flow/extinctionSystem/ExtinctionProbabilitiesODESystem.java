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

import java.util.List;

/**
 * This class represents the ODE system for extinction probabilities. Its state is the survival probability
 * q = 1 - p0 instead of the extinction probability p0, so that the relative tolerance applies to q. Otherwise,
 * the absolute error of p0 can exceed q when q is tiny (e.g. for strongly subcritical early intervals), which
 * makes 1 - p0 meaningless.
 */
public class ExtinctionProbabilitiesODESystem extends IntervalODESystem {

    private final double[][] birthRates;
    private final double[][] deathRates;
    private final double[][] samplingRates;
    private final double[][][] crossBirthRates;
    private final double[][][] migrationRates;

    private int currentInterval;

    public ExtinctionProbabilitiesODESystem(Parameterization parameterization, List<Interval> intervals, double absoluteTolerance, double relativeTolerance) {
        super(parameterization, intervals, absoluteTolerance, relativeTolerance);

        this.birthRates = this.parameterization.getBirthRates();
        this.deathRates = this.parameterization.getDeathRates();
        this.samplingRates = this.parameterization.getSamplingRates();
        this.crossBirthRates = this.parameterization.getCrossBirthRates();
        this.migrationRates = this.parameterization.getMigRates();

        this.currentInterval = intervals.size() - 1;
    }

    @Override
    public int getDimension() {
        return this.parameterization.getNTypes();
    }

    @Override
    public void computeDerivatives(double t, double[] y, double[] yDot) {
        if (Double.isNaN(t)) {
            throw new IllegalStateException("NaN detected during integration.");
        }

        int interval = this.currentInterval;

        for (int i = 0; i < this.parameterization.getNTypes(); i++) {
            double birthRate = this.birthRates[interval][i];
            double deathRate = this.deathRates[interval][i];
            double samplingRate = this.samplingRates[interval][i];

            // derivative of q = 1 - p0, written without the cancellation in 1 - p0
            yDot[i] = birthRate * y[i] * y[i] - (birthRate - deathRate - samplingRate) * y[i] - samplingRate;

            for (int j = 0; j < this.parameterization.getNTypes(); j++) {
                if (i == j) {
                    continue;
                }

                yDot[i] -= (
                        this.crossBirthRates[interval][i][j] * (1 - y[i]) * y[j]
                                + this.migrationRates[interval][i][j] * (y[j] - y[i])
                );
            }
        }
    }

    @Override
    protected void handleParameterizationIntervalBoundary(double boundaryTime, int oldInterval, int newInterval, double[] state) {
        super.handleParameterizationIntervalBoundary(boundaryTime, oldInterval, newInterval, state);

        // include rho sampling effects

        for (int type = 0; type < this.parameterization.getNTypes(); type++) {
            double rho = this.parameterization.getRhoValues()[newInterval][type];
            state[type] += rho * (1.0 - state[type]);
        }

        this.currentInterval--;
    }

}
