package beast.math;


import java.util.Arrays;

import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
import org.apache.commons.math3.ode.FirstOrderIntegrator;
import org.apache.commons.math3.ode.nonstiff.ClassicalRungeKuttaIntegrator;
import org.apache.commons.math3.ode.nonstiff.DormandPrince54Integrator;
import org.apache.commons.math3.ode.nonstiff.DormandPrince853Integrator;
import org.apache.commons.math3.ode.nonstiff.*;

import beast.core.util.Utils;

/**
 * User: Denise
 * Date: Jul 11, 2013
 * Time: 5:56:21 PM
 */


public class p0ge_ODE implements FirstOrderDifferentialEquations {

	p0_ODE P;
	public FirstOrderIntegrator p_integrator;

	Double[] b;
	Double[] b_ij;
	Double[] d;
	Double[] s;

	Boolean augmented;

	Double[] M;
	double T;

	int dimension; /* ODE dimension = stateNumber */
	int intervals;
	Double[] times;
	int index;

	int maxEvals;
	public int maxEvalsUsed;

	// implementation with a scale factor switched on or off depending on smallNumber
	int factorP;
	boolean smallNumber;

	public p0ge_ODE(Double[] b, Double[] b_ij, Double[] d, Double[] s, Double[] M, int dimension, int intervals, double T, Double[] times, p0_ODE P, int maxEvals, Boolean augmented, boolean SN){


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

		// factorP represents the order of magnitude by which initial conditions for the ode (on p equations) were increased (in order to prevent underflowing)
		// only used if smallNumber = true
		this.factorP = P.factor;
		this.smallNumber = SN;
	}

	public int getDimension() {
		return 2*this.dimension;
	}

	public void computeDerivatives(double t, double[] g, double[] gDot) {

		double scaleP = 0;

		// compute the actual scale factor from 'factorP'
		if (smallNumber) scaleP = Math.pow(10, factorP);

		index = Utils.index(t, times, intervals);

		int k, l;

		for (int i=0; i<dimension; i++){

			/*  p0 equations (0 .. dim-1) */

			k = i*intervals + index;

			if (smallNumber) {
				// scaling factor applied to master equations p0
				gDot[i] = + (b[k]+d[k]+s[k])*g[i]
						- (b[k] * g[i] * g[i])/scaleP
						- d[k]*scaleP ;	
			} else {
				gDot[i] = + (b[k]+d[k]+s[k]//)*g[i]
						- b[k] * g[i]) * g[i]
								- d[k] ;
			}


			for (int j=0; j<dimension; j++){

				l = (i*(dimension-1)+(j<i?j:j-1))*intervals + index;

				if (i!=j){
					if (b_ij!=null){     // infection among demes


						gDot[i] += b_ij[l]*g[i]; // b_ij[intervals*i*(dimension-1)+(j<i?j:j-1)+index]*g[i];
						if (!smallNumber)
							gDot[i] -= b_ij[l]*g[i]*g[j]; // b_ij[intervals*i*(dimension-1)+(j<i?j:j-1)+index]*g[i]*g[j];
						else
							gDot[i] -= (b_ij[l]*g[i]*g[j])/scaleP;
					}

					// migration:
					gDot[i] += M[l]*g[i];
					gDot[i] -= M[l]*g[j];

				}
			}

			/*  ge equations: (dim .. 2*dim-1) */

			if (smallNumber) {
				gDot[dimension+i] = + (b[k]+d[k]+s[k]
						- 2*b[k]*g[i]/scaleP)*g[dimension+i];
			} else {
				gDot[dimension+i] = + (b[k]+d[k]+s[k] //)*g[dimension+i]
						- 2*b[k]*g[i])*g[dimension+i];
			}


			for (int j=0; j<dimension; j++){

				l = (i*(dimension-1)+(j<i?j:j-1))*intervals + index;

				if (i!=j){

					if (b_ij!=null){     // infection among demes


						gDot[dimension+i] += b_ij[l]*g[dimension+i];
						if (!augmented) {
							if (smallNumber)
								gDot[dimension+i] -= b_ij[l]* ( (g[i]*g[dimension+j])/scaleP + (g[j]*g[dimension+i])/scaleP);
							else
								gDot[dimension+i] -= b_ij[l]* ( g[i]*g[dimension+j] + g[j]*g[dimension+i]);
						}
					}

					// migration:
					gDot[dimension+i] +=  M[l]*g[dimension+i];
					if (!augmented) gDot[dimension+i] -=  M[l]*g[dimension+j];
				}
			}

		}

		//        gDot[2] = -(-(b[0]+b_ij[0]+d[0]+s[0])*g[2] + 2*b[0]*g[0]*g[2] + b_ij[0]*g[0]*g[3] + b_ij[0]*g[1]*g[2]);
		//        gDot[3] = -(-(b[1]+b_ij[1]+d[1]+s[1])*g[3] + 2*b[1]*g[1]*g[3] + b_ij[1]*g[1]*g[2] + b_ij[1]*g[0]*g[3]);

	}

	/**
	 * Perform integration on differential equations p using a classical implementation (using double[] for the initial conditions and the output).
	 * WARNING: getP and getPSmallNumber are very similar. A correction made in one of the two would likely be needed in the other also.
	 * @param t
	 * @param rhoSampling
	 * @param rho
	 * @return
	 */
	public double[] getP(double t, Boolean rhoSampling, Double[] rho){

		double[] y = new double[dimension];

		Arrays.fill(y,1.);   // initial condition: y_i[T]=1 for all i

		if (rhoSampling)
			for (int i = 0; i<dimension; i++)
				y[i]-= rho[i*intervals + Utils.index(t, times, intervals)];    // initial condition: y_i[T]=1-rho_i

		if (Math.abs(T-t)<1e-10 ||  T < t) {
			return y;
		}

		try {

			double from = t;
			double to = T;
			double oneMinusRho;

			int indexFrom = Utils.index(from, times, times.length);
			int index = Utils.index(to, times, times.length);

			int steps = index - indexFrom;
			index--;
			if (Math.abs(from-times[indexFrom])<1e-10) steps--;
			if (index>0 && Math.abs(to-times[index-1])<1e-10) {
				steps--;
				index--;
			}

			while (steps > 0){

				from = times[index];//  + 1e-14;


				p_integrator.integrate(P, to, y, from, y); // solve P , store solution in y

				if (rhoSampling){
					for (int i=0; i<dimension; i++){
						oneMinusRho = (1-rho[i*intervals + Utils.index(times[index], times, intervals)]);
						y[i] *= oneMinusRho;
					}
				}

				to = times[index];


				steps--;
				index--;
			}


			p_integrator.integrate(P, to, y, t, y); // solve P, store solution in y

		}catch(Exception e){

			throw new RuntimeException("couldn't calculate p");
		}

		return y;
	}

	/**
	 * Implementation of getP with Small Numbers, to avoid potential underflowing.
	 * WARNING: getP and getPSmallNumber are very similar. A correction made in one of the two would likely be needed in the other also.
	 * @param t
	 * @param rhoSampling
	 * @param rho
	 * @return
	 */
	public SmallNumber[] getPSmallNumber(double t, Boolean rhoSampling, Double[] rho){


		double[] y = new double[dimension];

		Arrays.fill(y,1.);   // initial condition: y_i[T]=1 for all i


		if (rhoSampling)
			for (int i = 0; i<dimension; i++)
				y[i]-= rho[i*intervals + Utils.index(t, times, intervals)];    // initial condition: y_i[T]=1-rho_i

		// create an array of Small Numbers which will encapsulate the result of getP, to avoid any risk of underflowing
		SmallNumber[] ySmall = new SmallNumber[dimension];

		for (int i=0; i<dimension; i++){
			ySmall[i] = new SmallNumber(y[i]);
		}

		if (Math.abs(T-t)<1e-10 ||  T < t) {
			return ySmall;
		}

		try {
			ScaledNumbers yScaled = new ScaledNumbers();

			// yScaled contains the set of initial conditions scaled to fit the requirements on the values 'double' can represent, with the factor by which the numbers were multiplied 
			yScaled = SmallNumberScaler.scale(ySmall, EquationType.EquationOnP);

			double from = t;
			double to = T;
			double oneMinusRho;

			int indexFrom = Utils.index(from, times, times.length);
			int index = Utils.index(to, times, times.length);

			int steps = index - indexFrom;
			index--;
			if (Math.abs(from-times[indexFrom])<1e-10) steps--;
			if (index>0 && Math.abs(to-times[index-1])<1e-10) {
				steps--;
				index--;
			}

			// integrationResults is used temporarily to store the results of each integration step
			double[] integrationResults = new double[dimension];

			while (steps > 0){

				from = times[index];//  + 1e-14;

				P.setFactor(yScaled.getScalingFactor()[0]); // set the right scaling factor for P equations
				p_integrator.integrate(P, to, yScaled.getEquation(), from, integrationResults); // solve P , store solution in integrationResults then in y, when unscaled
				// 'unscale' values in integrationResults so as to retrieve accurate values after the integration.
				ySmall = SmallNumberScaler.unscale(integrationResults, yScaled.getScalingFactor(), EquationType.EquationOnP);



				if (rhoSampling){
					for (int i=0; i<dimension; i++){
						oneMinusRho = (1-rho[i*intervals + Utils.index(times[index], times, intervals)]);
						ySmall[i] = ySmall[i].scalarMultiply(oneMinusRho);
					}
				}

				to = times[index];

				steps--;
				index--;

				// 'rescale' the results of the last integration to prepare for the next integration step 
				yScaled = SmallNumberScaler.scale(ySmall, EquationType.EquationOnP);
			}

			P.setFactor(yScaled.getScalingFactor()[0]); // set the right scaling factor for P equations
			p_integrator.integrate(P, to, yScaled.getEquation(), t, integrationResults); // solve P, store solution in y
			ySmall = SmallNumberScaler.unscale(integrationResults, yScaled.getScalingFactor(), EquationType.EquationOnP);


		}catch(Exception e){
			// e.printStackTrace();
			throw new RuntimeException("couldn't calculate p");
		}

		return ySmall;
	}

	public void setFactor(int factorP){
		this.factorP = factorP;
		this.P.setFactor(factorP);
	}

	public int getFactor(){
		return this.factorP;
	}

	public void setSmallNumber(boolean SN){
		this.smallNumber = SN;
		this.P.smallNumber = SN;
	}

	/**
	 * This method serves as a comparison for various integrators
	 */
	public static void testCorrelations(){

		Double[] b;
		Double[] d = {1.,1.};
		Double[] s;
		Double[] M;// = {3.,3.};

		Double psi;

		Double c1 = 0.01;
		Double c2 = 0.1;
		
		
		for (double i =1.1; i<2; i+=0.125){

			b = new Double[]{i, i};

			//            psi = 0.5 * ((i - d[0]) - Math.sqrt((d[0] - i) * (d[0] - i) - .04));  // assume b*s*m=constant


			M = new Double[]{b[0]-d[0]-c2/b[0], b[0]-d[0]-c2/b[0]};     // assume b-d-s=M

			psi = c2/c1 * M[0]; // assume b*m = c1 and b*s = c2
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

			p0_ODE p_ode = new p0_ODE(b,null, d,s,M, 2, 1, new Double[]{0.}, 0, false);
			p0ge_ODE pg_ode = new p0ge_ODE(b,null, d,s,M, 2, 1, T, new Double[]{0.}, p_ode, Integer.MAX_VALUE,augmented, false);
			
			System.out.println("b[0] = "+b[0]+ ", d[0] = " + Math.round(d[0]*100.)/100.+ "\t\t");

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
		
		//testCorrelations();

		Double[] birth = {2.,2.};
		Double[] b;
		Double[] d = {.5,.5};
		Double[] s = {.5,.5};
		Double[] M = {0.,0.};

		System.out.println("b\tp\tg");

		//         for (double i = 0.; i<10.01; i+=0.01){

		int i = 1;

		b = new Double[]{i*birth[0], i*birth[1]};

		FirstOrderIntegrator integrator = new DormandPrince853Integrator(1.0e-4, 1., 1.0e-6, 1.0e-6);//new ClassicalRungeKuttaIntegrator(.01); //

		double T = 10.;
		Boolean augmented = false;

		p0_ODE p_ode = new p0_ODE(b,new Double[]{1.,1.}, d,s,M, 2, 1, new Double[]{0.}, 0, false);
		p0ge_ODE pg_ode = new p0ge_ODE(b,new Double[]{1.,1.}, d,s,M, 2, 1, T, new Double[]{0.}, p_ode, Integer.MAX_VALUE,augmented, false);


		double[] p0 = new double[]{1.,1.};
		double[] p = new double[2];
		double[] y0 = new double[]{1.,1.,1.,1.};
		double[] y = new double[4];

		integrator.integrate(pg_ode, T, y0, 0., y);
		integrator.integrate(p_ode, T, p0, 0., p);

		//             System.out.print("b[0] = "+b[0]+ ", d[0] = " + Math.round(d[0]*100.)/100.+ "\t\t");
		System.out.println(b[0] + "\t" + p[0]+"\t"+p[1]);
		System.out.println(b[0] + "\t" + y[0]+"\t"+y[1]+"\t"+y[2]+"\t"+y[3]);
		//         }
	}
}

