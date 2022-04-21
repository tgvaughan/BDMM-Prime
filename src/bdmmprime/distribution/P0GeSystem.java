package bdmmprime.distribution;


import bdmmprime.parameterization.Parameterization;
import bdmmprime.util.Utils;
import org.apache.commons.math3.ode.FirstOrderIntegrator;
import org.apache.commons.math3.ode.nonstiff.DormandPrince54Integrator;

/**
 * User: Denise
 * Date: Jul 11, 2013
 * Time: 5:56:21 PM
 */

public class P0GeSystem extends P0System {

	public P0GeSystem(Parameterization parameterization,
                      double absoluteTolerance,
                      double relativeTolerance) {

	    super(parameterization, absoluteTolerance, relativeTolerance);
	}

	@Override
	public int getDimension() {
		return 2*this.nTypes;
	}

	@Override
	public void computeDerivatives(double t, double[] y, double[] yDot) {

		for (int i = 0; i<nTypes; i++){

			/*  p0 equations (0 .. dim-1) */

			yDot[i] = (d[interval][i]+s[interval][i])*y[i] - d[interval][i];

			for (int j = 0; j < nTypes; j++){

			    if (i != j) {
                    yDot[i] += M[interval][i][j] * y[i];
                    yDot[i] -= M[interval][i][j] * y[j];
                }

                for (int k=0; k<=j; k++) {
                    yDot[i] += b[interval][i][j][k] * y[i];
                    yDot[i] -= b[interval][i][j][k] * y[j] * y[k];
                }
			}

			/*  ge equations: (dim .. 2*dim-1) */

			yDot[nTypes + i] = (d[interval][i]+s[interval][i])*y[nTypes+i];

			for (int j = 0; j<nTypes; j++) {
                if (j!=i) {
                    yDot[nTypes+i] += M[interval][i][j]*y[nTypes+i]
                            - M[interval][i][j]*y[nTypes+j];
                }

                for (int k=0; k<=j; k++) {
                    yDot[nTypes+i] += b[interval][i][j][k]*y[nTypes+i]
                            - b[interval][i][j][k]*(y[j]*y[nTypes+k] + y[k]*y[nTypes+j]);
                }

			}
		}
	}

    /**
     * Perform the integration of PG with initial conds in pgScaled between to and from
     * Use an adaptive-step-size integrator
     * "Safe" because it divides the integration interval in two
     * if the interval is (arbitrarily) judged to be too big to give reliable results
     *
     * @param pgScaled
     * @param tStart
     * @param tEnd
     * @return result of integration
     */
    public ScaledNumbers safeIntegrate(ScaledNumbers pgScaled, double tStart, double tEnd) {

        // if the integration interval is too small, nothing is done (to prevent infinite looping)
        if (Utils.equalWithPrecision(tEnd, tStart))
            return pgScaled;

        //TODO make threshold a class field
        if (totalProcessLength > 0 && Math.abs(tEnd - tStart) > totalProcessLength / 6) {
            pgScaled = safeIntegrate(pgScaled, tStart, tEnd + (tStart - tEnd) / 2);
            pgScaled = safeIntegrate(pgScaled, tEnd + (tStart - tEnd) / 2, tEnd);
        } else {

            //setup of the relativeTolerance and absoluteTolerance input of the adaptive integrator
            //TODO set these two as class fields
            double relativeToleranceConstant = 1e-7;
            double absoluteToleranceConstant = 1e-100;
            double[] absoluteToleranceVector = new double[2 * nTypes];
            double[] relativeToleranceVector = new double[2 * nTypes];

            for (int i = 0; i < nTypes; i++) {
                absoluteToleranceVector[i] = absoluteToleranceConstant;
                if (pgScaled.getEquation()[i + nTypes] > 0) { // adapt absoluteTolerance to the values stored in pgScaled
                    absoluteToleranceVector[i + nTypes] = Math.max(1e-310, pgScaled.getEquation()[i + nTypes] * absoluteToleranceConstant);
                } else {
                    absoluteToleranceVector[i + nTypes] = absoluteToleranceConstant;
                }
                relativeToleranceVector[i] = relativeToleranceConstant;
                relativeToleranceVector[i + nTypes] = relativeToleranceConstant;
            }

            double[] integrationResults = new double[pgScaled.getEquation().length];
            int a = pgScaled.getScalingFactor(); // store scaling factor
            int n = pgScaled.getEquation().length / 2; // dimension of the ODE system


            FirstOrderIntegrator integrator = new DormandPrince54Integrator(
                    integrationMinStep, integrationMaxStep,
                    absoluteToleranceVector, relativeToleranceVector);
            integrator.integrate(this, tStart, pgScaled.getEquation(), tEnd, integrationResults); // perform the integration step

            double[] pConditions = new double[n];
            SmallNumber[] geConditions = new SmallNumber[n];
            for (int i = 0; i < n; i++) {
                pConditions[i] = integrationResults[i];
                geConditions[i] = new SmallNumber(integrationResults[i + n]);
            }
            pgScaled = (new P0GeState(pConditions, geConditions)).getScaledState();
            pgScaled.augmentFactor(a);
        }

        return pgScaled;
    }
}

