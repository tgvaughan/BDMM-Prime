package bdmm.distributions;


import bdmm.util.Utils;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
import org.apache.commons.math3.ode.FirstOrderIntegrator;
import org.apache.commons.math3.ode.nonstiff.ClassicalRungeKuttaIntegrator;
import org.apache.commons.math3.ode.nonstiff.DormandPrince54Integrator;
import org.apache.commons.math3.ode.nonstiff.DormandPrince853Integrator;
import org.apache.commons.math3.ode.nonstiff.HighamHall54Integrator;

import java.util.Arrays;

/**
 * User: Denise
 * Date: Jul 11, 2013
 * Time: 5:56:21 PM
 */


public class p0ge_ODE implements FirstOrderDifferentialEquations {

	p0_ODE P;
	public FirstOrderIntegrator p_integrator;

	double[] b;
	Double[] b_ij;
	double[] d;
	Double[] s;

	Boolean augmented;
	Boolean birthAmongDemes;

	Double[] M;
	double T;

	int dimension; /* ODE numberOfDemes = stateNumber */
	int intervals;
	Double[] times;
	int index;

	int maxEvals;
	public int maxEvalsUsed;
	public static double globalPrecisionThreshold;


	public p0ge_ODE(double[] b, Double[] b_ij, double[] d, Double[] s, Double[] M, int dimension, int intervals, double T, Double[] times, p0_ODE P, int maxEvals, Boolean augmented){


		this.b = b;
		this.b_ij = b_ij;
		this.d = d;
		this.s = s;
		this.M = M;
		this.dimension = dimension;
		this.intervals = intervals;
		this.maxEvals = maxEvals;
		maxEvalsUsed = 0;

		this.T = T;
		this.times= times;
		this.P = P;

		this.augmented = augmented;
		this.birthAmongDemes = b_ij!=null;

	}

	public int getDimension() {
		return 2*this.dimension;
	}

	public void computeDerivatives(double t, double[] g, double[] gDot) {


		index = Utils.index(t, times, intervals);

		int k, l;

		for (int i=0; i<dimension; i++){

			/*  p0 equations (0 .. dim-1) */

			k = i*intervals + index;

			gDot[i] = + (b[k]+d[k]+s[k]
					- b[k] * g[i]) * g[i]
					- d[k] ;

			for (int j=0; j<dimension; j++){

				l = (i*(dimension-1)+(j<i?j:j-1))*intervals + index;

				if (i!=j){
					if (birthAmongDemes){     // infection among demes

						gDot[i] += b_ij[l]*g[i]; // birthAmongDemes_ij[intervals*i*(numberOfDemes-1)+(j<i?j:j-1)+indexTimeInterval]*g[i];

						gDot[i] -= b_ij[l]*g[i]*g[j]; // birthAmongDemes_ij[intervals*i*(numberOfDemes-1)+(j<i?j:j-1)+indexTimeInterval]*g[i]*g[j];

					}

					if (M[0]!=null) {// migration:
						gDot[i] += M[l] * g[i];
						gDot[i] -= M[l] * g[j];
					}
				}
			}

			/*  ge equations: (dim .. 2*dim-1) */


			gDot[dimension+i] = + (b[k]+d[k]+s[k]
					- 2*b[k]*g[i])*g[dimension+i];


			for (int j=0; j<dimension; j++){

				l = (i*(dimension-1)+(j<i?j:j-1))*intervals + index;

				if (i!=j){

					if (b_ij!=null){     // infection among demes

						gDot[dimension+i] += b_ij[l]*g[dimension+i];
						if (!augmented) {
							gDot[dimension+i] -= b_ij[l]* ( g[i]*g[dimension+j] + g[j]*g[dimension+i]);
						}
					}

					if (M[0]!=null) {// migration:
						gDot[dimension + i] += M[l] * g[dimension + i];
						if (!augmented) gDot[dimension + i] -= M[l] * g[dimension + j];
					}
				}
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


		if (Math.abs(T-t)<globalPrecisionThreshold || Math.abs(t0-t)<globalPrecisionThreshold ||   T < t)
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
						for (int i=0; i<dimension; i++){
							oneMinusRho = (1-rho[i*intervals + index]);
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
					for (int i=0; i<dimension; i++){
						oneMinusRho = (1-rho[i*intervals + indexFrom]);
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

		double[] y = new double[dimension];

		if (!rhoSampling)
			Arrays.fill(y,1.);   // initial condition: y_i[T]=1 for all i

		else{
			for (int i = 0; i<dimension; i++) {
				y[i] = (1 - rho[i * intervals + Utils.index(T, times, intervals)]);    // initial condition: y_i[T]=1-rho_i

				/*
				System.out.println("In getP, multiplying with oneMinusRho: " + (1 - rho[i * intervals + Utils.indexTimeInterval(T, rateChangeTimes, intervals)]) + ", t = " + t + ", to = " + T);
				*/
			}
		}

		if (Math.abs(T-t)<globalPrecisionThreshold ||  T < t) {
			return y;
		}

		return getP(t, y, T, rhoSampling, rho);

	}


	/**
	 * This method serves as a comparison for various integrators
	 */
	public static void testCorrelations(){

		double[] b;
		double[] d = {1.,1.};
		Double[] s;
		Double[] M;// = {3.,3.};

		Double psi;

		Double c1 = 0.01;
		Double c2 = 0.1;


		for (double i =1.1; i<2; i+=0.125){

			b = new double[]{i, i};

			//            psi = 0.5 * ((i - death[0]) - Math.sqrt((death[0] - i) * (death[0] - i) - .04));  // assume birth*sampling*m=constant


			M = new Double[]{b[0]-d[0]-c2/b[0], b[0]-d[0]-c2/b[0]};     // assume birth-death-sampling=migration

			psi = c2/c1 * M[0]; // assume birth*m = c1 and birth*sampling = c2
			s = new Double[] {psi,psi};


			FirstOrderIntegrator integrator1 = new ClassicalRungeKuttaIntegrator(.01);
			FirstOrderIntegrator integrator2 = new DormandPrince853Integrator(1.0e-10, 1., 1.0e-100, 1.0e-10);
			FirstOrderIntegrator integrator3 = new DormandPrince853Integrator(1.0e-10, 1., 1.0e-100, 1.0e-7);
			FirstOrderIntegrator integrator4 = new DormandPrince54Integrator(1.0e-10, 1., 1.0e-100, 1.0e-7);
			FirstOrderIntegrator integrator5 = new HighamHall54Integrator(1.0e-10, 1., 1.0e-14, 1.0e-7);
			FirstOrderIntegrator integrator6 = new HighamHall54Integrator(1.0e-20, 1., 1.0e-320, 1.0e-10);
			//FirstOrderIntegrator integrator7 = new DormandPrince54Integrator(1.0e-20, 1., 1.0e-10, 1.0e-9);
			//FirstOrderIntegrator integrator5 = new DormandPrince853Integrator(1.0e-10, 1., 1.0e-320, 1.0e-20);//new ClassicalRungeKuttaIntegrator(.01); //


			double T = 1;
			Boolean augmented = true;

			p0_ODE p_ode = new p0_ODE(b,null, d,s,M, 2, 1, new Double[]{0.});
			p0ge_ODE pg_ode = new p0ge_ODE(b,null, d,s,M, 2, 1, T, new Double[]{0.}, p_ode, Integer.MAX_VALUE,augmented);

			System.out.println("birth[0] = "+b[0]+ ", death[0] = " + Math.round(d[0]*100.)/100.+ "\t\t");

			double[] y0 = new double[]{1.1,1.1,4.,9.};
			double[] y = new double[4];
			//
			long start = System.nanoTime();
			integrator1.integrate(pg_ode, T, y0, 0., y);
			long end = System.nanoTime();
			long microseconds = (end - start) / 1000;
			System.out.println(y[0]+"\t"+y[1]+"\t" +y[2]+"\t"+y[3]+"\t"+microseconds+" \u03BCs");

			y0 = new double[]{1.1,1.1,4.,9.};
			y = new double[4];

			start = System.nanoTime();
			integrator2.integrate(pg_ode, T, y0, 0., y);
			end = System.nanoTime();
			microseconds = (end - start) / 1000;
			System.out.println(y[0]+"\t"+y[1]+"\t" +y[2]+"\t"+y[3]+"\t"+microseconds+" \u03BCs");

			y0 = new double[]{1.1,1.1,4.,9.};
			y = new double[4];

			start = System.nanoTime();
			integrator3.integrate(pg_ode, T, y0, 0., y);
			end = System.nanoTime();
			microseconds = (end - start) / 1000;
			System.out.println(y[0]+"\t"+y[1]+"\t" +y[2]+"\t"+y[3]+"\t"+microseconds+" \u03BCs");



			y0 = new double[]{1.1,1.1,4.,9.};
			y = new double[4];

			start = System.nanoTime();
			integrator4.integrate(pg_ode, T, y0, 0., y);
			end = System.nanoTime();
			microseconds = (end - start) / 1000;
			System.out.println(y[0]+"\t"+y[1]+"\t" +y[2]+"\t"+y[3]+"\t"+microseconds+" \u03BCs");

			y0 = new double[]{1.1,1.1,4.,9.};
			y = new double[4];

			start = System.nanoTime();
			integrator5.integrate(pg_ode, T, y0, 0., y);
			end = System.nanoTime();
			microseconds = (end - start) / 1000;
			System.out.println(y[0]+"\t"+y[1]+"\t" +y[2]+"\t"+y[3]+"\t"+microseconds+" \u03BCs");

			y0 = new double[]{1.1,1.1,4.,9.};
			y = new double[4];

			start = System.nanoTime();
			integrator6.integrate(pg_ode, T, y0, 0., y);
			end = System.nanoTime();
			microseconds = (end - start) / 1000;
			System.out.println(y[0]+"\t"+y[1]+"\t" +y[2]+"\t"+y[3]+"\t"+microseconds+" \u03BCs");

		}
	}

	public static void main(String args[]){

		testCorrelations();

		double[] birth = {2.,2.};
		double[] b;
		double[] d = {.5,.5};
		Double[] s = {.5,.5};
		Double[] M = {0.,0.};

		System.out.println("birth\tp\tg");

		//         for (double i = 0.; i<10.01; i+=0.01){

		int i = 1;

		b = new double[]{i*birth[0], i*birth[1]};

		FirstOrderIntegrator integrator = new DormandPrince853Integrator(1.0e-4, 1., 1.0e-6, 1.0e-6);//new ClassicalRungeKuttaIntegrator(.01); //

		double T = 10.;
		Boolean augmented = false;

		p0_ODE p_ode = new p0_ODE(b,new Double[]{1.,1.}, d,s,M, 2, 1, new Double[]{0.});
		p0ge_ODE pg_ode = new p0ge_ODE(b,new Double[]{1.,1.}, d,s,M, 2, 1, T, new Double[]{0.}, p_ode, Integer.MAX_VALUE,augmented);

		pg_ode.p_integrator = integrator;
		double[] p0 = new double[]{1.,1.};
		double[] p = new double[2];
		double[] y0 = new double[]{1.,1.,1.,1.};
		double[] y = new double[4];

		integrator.integrate(pg_ode, T, y0, 0., y);
		integrator.integrate(p_ode, T, p0, 0., p);

		double[] res2 = new double[] {4};
		double[] res  = pg_ode.getP(8, false, new Double[]{0.});
		System.out.println(b[0] + "\t" + res[0]+"\t"+res[1]);
		res2 = pg_ode.getP(5, res, 8, false, new Double[]{0.});
		System.out.println(b[0] + "\t" + res[0]+"\t"+res[1]);
		res2 = pg_ode.getP(0, res, 5, false, new Double[]{0.});
		System.out.println(b[0] + "\t" + res[0]+"\t"+res[1]);

		//             System.out.print("birth[0] = "+birth[0]+ ", death[0] = " + Math.round(death[0]*100.)/100.+ "\t\t");
		System.out.println(b[0] + "\t" + p[0]+"\t"+p[1]);
		System.out.println(b[0] + "\t" + y[0]+"\t"+y[1]+"\t"+y[2]+"\t"+y[3]);
		//         }

	}
}

