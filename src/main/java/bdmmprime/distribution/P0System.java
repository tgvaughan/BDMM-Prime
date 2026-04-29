/*
 * Copyright (C) 2016-2025 ETH Zurich
 * Copyright (C) 2013-2018 Denise Kühnert
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

import bdmmprime.parameterization.Parameterization;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
import org.apache.commons.math3.ode.FirstOrderIntegrator;
import org.apache.commons.math3.ode.nonstiff.DormandPrince54Integrator;


/**
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

		for (int i = 0; i< nTypes; i++){

			yDot[i] = + (b[interval][i]+d[interval][i]+s[interval][i])*y[i] - d[interval][i] - b[interval][i]*y[i]*y[i] ;

			for (int j = 0; j< nTypes; j++){

			    if (j==i)
			        continue;

                yDot[i] += b_ij[interval][i][j]*y[i];
                yDot[i] -= b_ij[interval][i][j]*y[i]*y[j];

                yDot[i] += M[interval][i][j] * y[i];
                yDot[i] -= M[interval][i][j] * y[j];
			}
		}

	}

	public void integrate(P0State state, double tStart, double tEnd) {
        p0Integrator.integrate(this, tStart, state.p0, tEnd, state.p0);
    }
}


