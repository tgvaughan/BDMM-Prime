/*
 * Copyright (c) 2012-2016 Denise Kühnert
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

package bdmmprime.distribution.classic;

import bdmmprime.parameterization.Parameterization;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
import org.apache.commons.math3.ode.FirstOrderIntegrator;
import org.apache.commons.math3.ode.nonstiff.DormandPrince54Integrator;


/**
 * ODE system for the extinction probabilities p0. Its state is the survival probability q = 1 - p0 instead of p0,
 * so that the relative tolerance applies to q. Otherwise, the absolute error of p0 can exceed q when q is tiny,
 * which makes 1 - p0 meaningless.
 *
 * @author dkuh004
 *         Date: May 24, 2012
 *         Time: 6:42:00 PM
 */

public class P0System implements FirstOrderDifferentialEquations {

	public double[][] b, d, s, r,rho;
	public double[][][] M, b_ij;

    public double totalProcessLength;

	public int nTypes;
	public int nIntervals;
    public double[] intervalEndTimes;

    protected int interval;

    protected FirstOrderIntegrator p0Integrator;

    protected double integrationMinStep, integrationMaxStep;


	public P0System(Parameterization parameterization,
                    double absoluteTolerance,
                    double relativeTolerance) {

		this.b = parameterization.getBirthRates();
		this.d = parameterization.getDeathRates();
		this.s = parameterization.getSamplingRates();
		this.r = parameterization.getRemovalProbs();
		this.rho = parameterization.getRhoValues();

		this.M = parameterization.getMigRates();
        this.b_ij = parameterization.getCrossBirthRates();

        this.totalProcessLength = parameterization.getTotalProcessLength();

		this.nTypes = parameterization.getNTypes();
		this.nIntervals = parameterization.getTotalIntervalCount();

        this.intervalEndTimes = parameterization.getIntervalEndTimes();

        integrationMinStep = parameterization.getTotalProcessLength() * 1e-100;
        integrationMaxStep= parameterization.getTotalProcessLength() / 10;

        this.p0Integrator = new DormandPrince54Integrator(
                integrationMinStep, integrationMaxStep,
                absoluteTolerance, relativeTolerance);
	}

	public void setInterval(int interval) {
	    this.interval = interval;
    }

	public int getDimension() {
		return this.nTypes;
	}

	public void computeDerivatives(double t, double[] y, double[] yDot) {

		for (int i = 0; i < this.nTypes; i++) {
			double birthRate = this.b[this.interval][i];
			double deathRate = this.d[this.interval][i];
			double samplingRate = this.s[this.interval][i];

			// derivative of q = 1 - p0, written without the cancellation in 1 - p0
			yDot[i] = birthRate * y[i] * y[i] - (birthRate - deathRate - samplingRate) * y[i] - samplingRate;

			for (int j = 0; j < this.nTypes; j++) {
			    if (j == i)
			        continue;

                yDot[i] -= this.b_ij[this.interval][i][j] * (1 - y[i]) * y[j];
                yDot[i] -= this.M[this.interval][i][j] * (y[j] - y[i]);
			}
		}

	}

	/**
	 * Integrates the survival probabilities in place from tStart to tEnd.
	 */
	public void integrate(double[] survivalProbability, double tStart, double tEnd) {
        this.p0Integrator.integrate(this, tStart, survivalProbability, tEnd, survivalProbability);
    }
}


