package bdmm.distributions;


import bdmm.parameterization.Parameterization;
import bdmm.util.Utils;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
import org.apache.commons.math3.ode.FirstOrderIntegrator;

/**
 * User: Denise
 * Date: Jul 11, 2013
 * Time: 5:56:21 PM
 */


public class P0GeSystem implements FirstOrderDifferentialEquations {

	public double[][] b, d, s, r, rho;
	public double[][][] b_ij, M;

	public double origin;

	public int nTypes; /* ODE numberOfDemes = stateNumber */
	public int nIntervals;
	public double[] intervalStartTimes;

	int maxEvals;
	public int maxEvalsUsed;
	public static double globalPrecisionThreshold;


	public P0GeSystem(Parameterization parameterization, int maxEvals){


		this.b = parameterization.getBirthRates();
		this.d = parameterization.getDeathRates();
		this.s = parameterization.getSamplingRates();
		this.r = parameterization.getRemovalProbs();
		this.rho = parameterization.getRhoValues();

        this.b_ij = parameterization.getCrossBirthRates();
		this.M = parameterization.getMigRates();

		this.nTypes = parameterization.getNTypes();
		this.nIntervals = parameterization.getTotalIntervalCount();

		this.maxEvals = maxEvals;
		maxEvalsUsed = 0;

		this.origin = parameterization.getOrigin();
		this.intervalStartTimes = parameterization.getIntervalStartTimes();
	}

	public int getDimension() {
		return 2*this.nTypes;
	}

	public void computeDerivatives(double t, double[] g, double[] gDot) {

		int interval = Utils.getIntervalIndex(t, intervalStartTimes);

		for (int i = 0; i<nTypes; i++){

			/*  p0 equations (0 .. dim-1) */

			gDot[i] = + (b[interval][i]+d[interval][i]+s[interval][i]
					- b[interval][i] * g[i]) * g[i]
					- d[interval][i] ;

			for (int j = 0; j < nTypes; j++){

			    if (i==j)
			        continue;


                gDot[i] += b_ij[interval][i][j]*g[i]; // birthAmongDemes_ij[intervals*i*(numberOfDemes-1)+(j<i?j:j-1)+indexTimeInterval]*g[i];

                gDot[i] -= b_ij[interval][i][j]*g[i]*g[j]; // birthAmongDemes_ij[intervals*i*(numberOfDemes-1)+(j<i?j:j-1)+indexTimeInterval]*g[i]*g[j];

                gDot[i] += M[interval][i][j] * g[i];
                gDot[i] -= M[interval][i][j] * g[j];
			}

			/*  ge equations: (dim .. 2*dim-1) */


			gDot[nTypes + i] = + (b[interval][i]+d[interval][i]+s[interval][i]
					- 2*b[interval][i]*g[i])*g[nTypes + i];


			for (int j = 0; j< nTypes; j++){

                if (i==j)
			        continue;

                gDot[nTypes +i] += b_ij[interval][i][j]*g[nTypes+i];
                gDot[nTypes +i] -= b_ij[interval][i][j]* ( g[i]*g[nTypes+j] + g[j]*g[nTypes+i]);

                gDot[nTypes + i] += M[interval][i][j] * g[nTypes + i];
                gDot[nTypes + i] -= M[interval][i][j] * g[nTypes + j];
			}

		}

	}
}

