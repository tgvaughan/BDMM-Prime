package bdmm.distributions;

import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;

import bdmm.util.Utils;


/**
 * @author dkuh004
 *         Date: May 24, 2012
 *         Time: 6:42:00 PM
 */

public class p0_ODE implements FirstOrderDifferentialEquations {

	public double[][] b, d, s, r,rho;
	public double[][][] M, b_ij;

	public int nTypes;
	public int nIntervals;

	public double[] times;

	public p0_ODE(Parameterization parameterization) {

		this.b = parameterization.getBirthRates();
		this.d = parameterization.getDeathRates();
		this.s = parameterization.getSamplingRates();
		this.r = parameterization.getRemovalProbs();
		this.rho = parameterization.getRhoValues();

		this.M = parameterization.getMigRates();
        this.b_ij = parameterization.getCrossBirthRates();

		this.nTypes = parameterization.getNTypes();
		this.nIntervals = parameterization.getTotalIntervalCount();

		this.times = parameterization.getIntervalStartTimes();

	}

	public int getDimension() {
		return this.nTypes;
	}

	public void computeDerivatives(double t, double[] y, double[] yDot) {

		int interval = Utils.index(t, times, nIntervals); //finds the indexTimeInterval of the time interval t lies in

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
}


