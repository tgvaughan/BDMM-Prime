package bdmm.distributions;


import bdmm.util.Utils;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
import org.apache.commons.math3.ode.FirstOrderIntegrator;

import java.util.Arrays;

/**
 * User: Denise
 * Date: Jul 11, 2013
 * Time: 5:56:21 PM
 */


public class p0ge_ODE implements FirstOrderDifferentialEquations {

	public FirstOrderIntegrator p_integrator;
    public p0_ODE P;

	public double[][] b, d, s;
	public double[][] b_ij, M;

	public double origin;

	public int nTypes; /* ODE numberOfDemes = stateNumber */
	public int nIntervals;
	public double[] times;

	int maxEvals;
	public int maxEvalsUsed;
	public static double globalPrecisionThreshold;


	public p0ge_ODE(Parameterization parameterization, p0_ODE P, int maxEvals){


		this.b = parameterization.getBirthRates();
		this.d = parameterization.getDeathRates();
		this.s = parameterization.getSamplingRates();

        this.b_ij = parameterization.getCrossBirthRates();
		this.M = parameterization.getMigRates();

		this.nTypes = parameterization.getNTypes();
		this.nIntervals = parameterization.getTotalIntervalCount();

		this.maxEvals = maxEvals;
		maxEvalsUsed = 0;

		this.origin = parameterization.getOrigin();
		this.times= parameterization.getIntervalStartTimes();

		this.P = P;
	}

	public int getDimension() {
		return 2*this.nTypes;
	}

	public void computeDerivatives(double t, double[] g, double[] gDot) {

		int interval = Utils.index(t, times, nIntervals);

		for (int i = 0; i<nTypes; i++){

			/*  p0 equations (0 .. dim-1) */

			gDot[i] = + (b[interval][i]+d[interval][i]+s[interval][i]
					- b[interval][i] * g[i]) * g[i]
					- d[interval][i] ;

			for (int j = 0; j < nTypes; j++){

			    if (i==j)
			        continue;

			    int ijOffset = Utils.getArrayOffset(i, j, nTypes);

                gDot[i] += b_ij[interval][ijOffset]*g[i]; // birthAmongDemes_ij[intervals*i*(numberOfDemes-1)+(j<i?j:j-1)+indexTimeInterval]*g[i];

                gDot[i] -= b_ij[interval][ijOffset]*g[i]*g[j]; // birthAmongDemes_ij[intervals*i*(numberOfDemes-1)+(j<i?j:j-1)+indexTimeInterval]*g[i]*g[j];

                gDot[i] += M[interval][ijOffset] * g[i];
                gDot[i] -= M[interval][ijOffset] * g[j];
			}

			/*  ge equations: (dim .. 2*dim-1) */


			gDot[nTypes + i] = + (b[interval][i]+d[interval][i]+s[interval][i]
					- 2*b[interval][i]*g[i])*g[nTypes + i];


			for (int j = 0; j< nTypes; j++){

                if (i==j)
			        continue;

			    int ijOffset = Utils.getArrayOffset(i, j, nTypes);

                gDot[nTypes +i] += b_ij[interval][ijOffset]*g[nTypes+i];
                gDot[nTypes +i] -= b_ij[interval][ijOffset]* ( g[i]*g[nTypes+j] + g[j]*g[nTypes+i]);

                gDot[nTypes + i] += M[interval][ijOffset] * g[nTypes + i];
                gDot[nTypes + i] -= M[interval][ijOffset] * g[nTypes + j];
			}

		}

		//        gDot[2] = -(-(birth[0]+birthAmongDemes_ij[0]+death[0]+sampling[0])*g[2] + 2*birth[0]*g[0]*g[2] + birthAmongDemes_ij[0]*g[0]*g[3] + birthAmongDemes_ij[0]*g[1]*g[2]);
		//        gDot[3] = -(-(birth[1]+birthAmongDemes_ij[1]+death[1]+sampling[1])*g[3] + 2*birth[1]*g[1]*g[3] + birthAmongDemes_ij[1]*g[1]*g[2] + birthAmongDemes_ij[1]*g[0]*g[3]);

	}

	/**
	 * Perform integration on differential equations p
	 * @param t
	 * @param rhoSampling
	 * @param rho
	 * @return
	 */
	public double[] getP(double t, double[]P0, double t0, Boolean rhoSampling, Double[] rho){


		if (Math.abs(origin -t)<globalPrecisionThreshold || Math.abs(t0-t)<globalPrecisionThreshold ||   origin < t)
			return P0;

		double[] result = new double[P0.length];

		try {

			System.arraycopy(P0, 0, result, 0, P0.length);
			double from = t;
			double to = t0;
			double oneMinusRho;

			int indexFrom = Utils.index(from, times, times.length);
			int index = Utils.index(to, times, times.length);

			int steps = index - indexFrom;

			index--;
			if (Math.abs(from-times[indexFrom])<globalPrecisionThreshold) steps--;
			if (index>0 && Math.abs(to-times[index-1])<globalPrecisionThreshold) {
				steps--;
				index--;
			}

			while (steps > 0){

				from = times[index];

				// TODO: putting the if(rhosampling) in there also means the 1-rho may never be actually used so a workaround is potentially needed
				if (Math.abs(from-to)>globalPrecisionThreshold){
					p_integrator.integrate(P, to, result, from, result); // solve diffEquationOnP , store solution in y

					if (rhoSampling){
						for (int i = 0; i< nTypes; i++){
							oneMinusRho = (1-rho[i* nIntervals + index]);
							result[i] *= oneMinusRho;

							/*
							System.out.println("In getP, multiplying with oneMinusRho: " + oneMinusRho + ", from = " + from);
							*/
						}
					}
				}

				to = times[index];

				steps--;
				index--;
			}

			p_integrator.integrate(P, to, result, t, result); // solve diffEquationOnP, store solution in y

			// TO DO
			// check that both rateChangeTimes are really overlapping
			// but really not sure that this is enough, i have to build appropriate tests
			if(Math.abs(t-times[indexFrom])<globalPrecisionThreshold) {
				if (rhoSampling){
					for (int i = 0; i< nTypes; i++){
						oneMinusRho = (1-rho[i* nIntervals + indexFrom]);
						result[i] *= oneMinusRho;
						//	System.out.println("In getP, multiplying as the final step with oneMinusRho: " + oneMinusRho + ",  = " + t);

					}
				}
			}


		}catch(Exception e){

			throw new RuntimeException("couldn't calculate p");
		}

		return result;
	}


	public double[] getP(double t, Boolean rhoSampling, Double[] rho){

		double[] y = new double[nTypes];

		if (!rhoSampling)
			Arrays.fill(y,1.);   // initial condition: y_i[T]=1 for all i

		else{
			for (int i = 0; i< nTypes; i++) {
				y[i] = (1 - rho[i * nIntervals + Utils.index(origin, times, nIntervals)]);    // initial condition: y_i[T]=1-rho_i

				/*
				System.out.println("In getP, multiplying with oneMinusRho: " + (1 - rho[i * intervals + Utils.indexTimeInterval(T, rateChangeTimes, intervals)]) + ", t = " + t + ", to = " + T);
				*/
			}
		}

		if (Math.abs(origin -t)<globalPrecisionThreshold ||  origin < t) {
			return y;
		}

		return getP(t, y, origin, rhoSampling, rho);

	}

}

