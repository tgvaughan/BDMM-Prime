package bdmm.distributions;

import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
import org.apache.commons.math3.ode.FirstOrderIntegrator;
import org.apache.commons.math3.ode.nonstiff.DormandPrince853Integrator;

import bdmm.util.Utils;


/**
 * @author dkuh004
 *         Date: May 24, 2012
 *         Time: 6:42:00 PM
 */

public class p0_ODE implements FirstOrderDifferentialEquations {

	double[][] b, d, s;
	double[][] M, b_ij;

	int nTypes;
	int intervalCount;

	double[] times;

	public p0_ODE(Parameterization parameterization) {

		this.b = parameterization.getBirthRates();
		this.b_ij = parameterization.getCrossBirthRates();
		this.d = parameterization.getDeathRates();
		this.s = parameterization.getSamplingRates();
		this.M = parameterization.getMigRates();

		this.nTypes = parameterization.getNTypes();
		this.intervalCount = parameterization.getTotalIntervalCount();

		this.times = parameterization.getIntervalStartTimes();

	}

	public int getDimension() {
		return this.nTypes;
	}

	public void computeDerivatives(double t, double[] y, double[] yDot) {

		int interval = Utils.index(t, times, intervalCount); //finds the indexTimeInterval of the time interval t lies in

		for (int i = 0; i< nTypes; i++){

			yDot[i] = + (b[interval][i]+d[interval][i]+s[interval][i])*y[i] - d[interval][i] - b[interval][i]*y[i]*y[i] ;

			for (int j = 0; j< nTypes; j++){

			    if (j==i)
			        continue;

				int ijOffset = Utils.getArrayOffset(i, j, nTypes);

                yDot[i] += b_ij[interval][ijOffset]*y[i];
                yDot[i] -= b_ij[interval][ijOffset]*y[i]*y[j];

                yDot[i] += M[interval][ijOffset] * y[i];
                yDot[i] -= M[interval][ijOffset] * y[j];
			}
		}

	}
}


