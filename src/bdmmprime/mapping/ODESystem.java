package bdmmprime.mapping;

import bdmmprime.parameterization.Parameterization;
import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.exception.MaxCountExceededException;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
import org.apache.commons.math3.ode.events.EventHandler;

/**
 * System of p0 and ge ODEs used in stochastic mapping.
 *
 * This class is kept distinct from the similar classes in the distributions
 * package just because the stochastic mapper needs are more basic.  We can
 * potentially merge these later on though.
 */
public class ODESystem implements FirstOrderDifferentialEquations, EventHandler {

    private final Parameterization param;
    private int interval;

    public ODESystem(Parameterization parameterization) {
        this.param = parameterization;
    }

    public void setInterval(int interval) {
        this.interval = interval;
    }


    /* FirstOrderDifferentialEquations implementation */

    @Override
    public int getDimension() {
        return param.getNTypes()*2;
    }

    @Override
    public void computeDerivatives(double t, double[] y, double[] yDot)
            throws MaxCountExceededException, DimensionMismatchException {

        int nTypes = param.getNTypes();

        for (int i = 0; i<nTypes; i++){

			/*  p0 equations (0 .. dim-1) */

			yDot[i] = (param.getDeathRates()[interval][i]
                    + param.getSamplingRates()[interval][i]) * y[i]
					- param.getDeathRates()[interval][i];

			for (int j = 0; j < nTypes; j++){

			    if (i!=j)
                    yDot[i] += param.getMigRates()[interval][i][j] * (y[i] - y[j]);

                for (int k=0; k<=j; k++) {
                    yDot[i] += param.getBirthRates()[interval][i][j][k]
                            * (y[i] - y[j] * y[k]);
                }
			}

			/*  ge equations: (dim .. 2*dim-1) */

			yDot[nTypes+i] = (param.getDeathRates()[interval][i]
                    + param.getSamplingRates()[interval][i])*y[nTypes+i];

			for (int j = 0; j<nTypes; j++) {
                if (i!=j) {
                    yDot[nTypes + i] += param.getMigRates()[interval][i][j]
                            * (y[nTypes + i] - y[nTypes + j]);
                }

                for (int k=0; k<=j; k++) {
                    yDot[nTypes + i] += param.getBirthRates()[interval][i][j][k]
                            * (y[nTypes+i] - (y[nTypes+j]*y[k] + y[nTypes+k]*y[j]));
                }
            }
		}

    }

    /* EventHandler Implementation */

    @Override
    public void init(double v, double[] doubles, double v1) { }

    @Override
    public double g(double v, double[] doubles) {
        double res = 1.0;
        for (double boundary : param.getIntervalEndTimes())
            res *= v - boundary;

        return res;
    }

    @Override
    public Action eventOccurred(double t, double[] y, boolean increasing) {
        return Action.RESET_STATE;
    }

    @Override
    public void resetState(double t, double[] y) {

        interval -= 1;

        // Include effect of any rho sampling times we pass:

        for (int type=0; type<param.getNTypes(); type++) {
            y[type] *= 1.0 - param.getRhoValues()[interval][type];
            y[type+param.getNTypes()] *= 1.0 - param.getRhoValues()[interval][type];
        }
    }
}
