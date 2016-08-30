package beast.math;

import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
import org.apache.commons.math3.ode.FirstOrderIntegrator;
import org.apache.commons.math3.ode.nonstiff.ClassicalRungeKuttaIntegrator;
import org.apache.commons.math3.ode.nonstiff.DormandPrince853Integrator;

import beast.core.util.Utils;


/**
 * @author dkuh004
 *         Date: May 24, 2012
 *         Time: 6:42:00 PM
 */

public class p0_ODE implements FirstOrderDifferentialEquations {

	Double[] b;
	Double[] b_ij;
	Double[] d;
	Double[] s;

	Double[] M;

	int dimension;
	int intervals;
	Double[] times;
	int index;



	// smallNumber switches on (if true) or off (if false) the implementation of likelihood calculations that prevents underflowing
	boolean smallNumber;

	// factor represents the order of magnitude by which initial conditions for the ode were increased (in order to prevent underflowing)
	// only used if smallNumber = true
	int factor;

	public p0_ODE(Double[] b, Double[] b_ij, Double[] d, Double[] s, Double[] M, int dimension , int intervals, Double[] times, int factor, boolean smallNumber){

		this.b = b;
		this.b_ij = b_ij;
		this.d = d;
		this.s = s;
		this.M = M;
		this.dimension = dimension;
		this.intervals = intervals;

		this.times = times;

		this.factor = factor;
		this.smallNumber = smallNumber;
	}

	public p0_ODE(Double[] b, Double[] b_ij, Double[] d, Double[] s, Double[] M, int dimension , int intervals, Double[] times) {

		this.b = b;
		this.b_ij = b_ij;
		this.d = d;
		this.s = s;
		this.M = M;
		this.dimension = dimension;
		this.intervals = intervals;

		this.times = times;

		this.factor = 0;
		this.smallNumber = false;

	}


	public void updateRates(Double[] b, Double[] b_ij, Double[] d, Double[] s, Double[] M, Double[] times){

		this.b = b;
		this.b_ij = b_ij;
		this.d = d;
		this.s = s;
		this.M = M;
		this.times = times;

	}

	public int getDimension() {
		return this.dimension;
	}

	public void setFactor(int factorP){
		this.factor = factorP;
	}

	public void computeDerivatives(double t, double[] y, double[] yDot) {

		index = Utils.index(t, times, intervals); //finds the index of the time interval t lies in 
		int k, l;

		double scale = 0;

		// compute the actual scale factor from 'factor'
		if (smallNumber) scale = Math.pow(10, factor);

		// if boolean smallNumber is false, use the classic method without scaling/unscaling

		for (int i = 0; i<dimension; i++){

			k = i*intervals + index;

			if (!smallNumber) yDot[i] = + (b[k]+d[k]+s[k])*y[i] - d[k] - b[k]*y[i]*y[i] ;
			else yDot[i] = + (b[k]+d[k]+s[k])*y[i] - d[k]*scale - b[k]*y[i]*y[i]/scale ;

			for (int j=0; j<dimension; j++){

				l = (i*(dimension-1)+(j<i?j:j-1))*intervals + index;

				if (i!=j){

					if (b_ij!=null){     // infection among demes


						yDot[i] += b_ij[l]*y[i]; 
						if (!smallNumber) yDot[i] -= b_ij[l]*y[i]*y[j];
						else yDot[i] -= b_ij[l]*y[i]*y[j]/scale;
					}

					// migration:
					yDot[i] += M[l]*y[i];;
					yDot[i] -= M[l]*y[j];
				}
			}
		}

	}


	public static void main(String[] args) throws Exception{

		// 2d test
		Double[] b = {1.03,1.06};
		Double[] d = {1.,1.};
		Double[] s = {0.02,0.04};
		Double[] M = new Double[]{3.,4.};

		FirstOrderIntegrator integrator = new DormandPrince853Integrator(1.0e-8, 100.0, 1.0e-10, 1.0e-10);//new ClassicalRungeKuttaIntegrator(.01); //
		FirstOrderDifferentialEquations ode = new p0_ODE(b,null,d,s,M, 2, 1, new Double[]{0.}, 0 , false);
		double[] y0 = new double[]{1.,1.};
		double[] y = new double[2];

		//        StepHandler stepHandler = new StepHandler() {
		//            public void init(double t0, double[] y0, double t) {
		//            }
		//
		//            public void handleStep(StepInterpolator interpolator, boolean isLast) {
		//                double   t = interpolator.getCurrentTime();
		//                double[] y = interpolator.getInterpolatedState();
		//                System.out.println(t + " " + y[0] + " " + y[1]);
		//            }
		//        };
		//        integrator.addStepHandler(stepHandler);
		//
		//
		integrator.integrate(ode, 10, y0, 1, y);

		System.out.println("Solution: " +y[0]+" "+y[1]);
		//
		//
		//        // 3d test
		//        Double[] b = {1.03,1.06, 1.5};
		//        Double[] d = {1.,1., 1.2};
		//        Double[] s = {0.02,0.04, 0.1};
		//        Double[] M = {3., 1.,4.,1.,2., 2.};
		//
		//        FirstOrderIntegrator integrator = new DormandPrince853Integrator(1.0e-8, 100.0, 1.0e-10, 1.0e-10);//new ClassicalRungeKuttaIntegrator(.01); //
		//        FirstOrderDifferentialEquations ode = new p0_ODE(b,d,s,M, 3);
		//        double[] y0 = new double[]{1.,1.,1.};
		//        double[] y = new double[3];
		//
		//        integrator.integrate(ode, 10, y0, 1, y);
		//
		//        System.out.println("Solution: " +y[0]+" "+y[1]+" "+y[2]);
	}

}


