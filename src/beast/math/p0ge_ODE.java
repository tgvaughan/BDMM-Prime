package beast.math;


import org.apache.commons.math3.ode.FirstOrderIntegrator;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
import org.apache.commons.math3.ode.nonstiff.DormandPrince853Integrator;

import java.util.Arrays;

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
//    Boolean birthChanges; // Does birth rate change over time, i.e. among intervals?
//    Boolean deathChanges; // Does death rate change over time, i.e. among intervals?
//    Boolean samplingChanges; // Does sampling rate change over time, i.e. between intervals?

    Boolean augmented;

    Double[] M;
    double T;

    int dimension; /* ODE dimension = stateNumber */
    int intervals;
    Double[] times;
    int index;

    int maxEvals;
    public int maxEvalsUsed;

    public p0ge_ODE(Double[] b, Double[] b_ij, Double[] d, Double[] s, Double[] M, int dimension, int intervals, double T, Double[] times, p0_ODE P, int maxEvals, Boolean augmented){


        this.b = b;
        this.b_ij = b_ij;
        this.d = d;
        this.s = s;
        this.M = M;
        this.dimension = dimension;
        this.intervals = intervals;
        this.maxEvals = maxEvals;
        maxEvalsUsed = 0;

//        birthChanges = (b.length == dimension*intervals);
//        deathChanges = (d.length == dimension*intervals);
//        samplingChanges = (s.length == dimension*intervals);

        this.T = T;
        this.times= times;
        this.P = P;

        this.augmented = augmented;

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

            gDot[i] = + (b[k]+d[k]+s[k]//)*g[i]
                    - b[k] * g[i]) * g[i]
                    - d[k] ;

            for (int j=0; j<dimension; j++){

                l = (i*(dimension-1)+(j<i?j:j-1))*intervals + index;

                if (i!=j){
                    if (b_ij!=null){     // infection among demes


                            gDot[i] += b_ij[l]*g[i]; // b_ij[intervals*i*(dimension-1)+(j<i?j:j-1)+index]*g[i];
                            gDot[i] -= b_ij[l]*g[i]*g[j]; // b_ij[intervals*i*(dimension-1)+(j<i?j:j-1)+index]*g[i]*g[j];
                    }

                    // migration:
                    gDot[i] += M[l]*g[i]; //M[i*(dimension-1)+(j<i?j:j-1)]*g[i];
                    gDot[i] -= M[l]*g[j]; //M[i*(dimension-1)+(j<i?j:j-1)]*g[j];

                }
            }

            /*  ge equations: (dim .. 2*dim-1) */


            gDot[dimension+i] = + (b[k]+d[k]+s[k] //)*g[dimension+i]
                    - 2*b[k]*g[i])*g[dimension+i];

            for (int j=0; j<dimension; j++){

                l = (i*(dimension-1)+(j<i?j:j-1))*intervals + index;

                if (i!=j){

                    if (b_ij!=null){     // infection among demes

                            gDot[dimension+i] += b_ij[l]*g[dimension+i];
                            if (!augmented)
                                gDot[dimension+i] -= b_ij[l]* ( g[i]*g[dimension+j] + g[j]*g[dimension+i]);
                    }

                    // migration:
                    gDot[dimension+i] +=  M[l]*g[dimension+i]; //M[i*(dimension-1)+(j<i?j:j-1)]*g[dimension+i];
                    if (!augmented) gDot[dimension+i] -=  M[l]*g[dimension+j]; // M[i*(dimension-1)+(j<i?j:j-1)]*g[dimension+j];
                }
            }

        }

//        gDot[2] = -(-(b[0]+b_ij[0]+d[0]+s[0])*g[2] + 2*b[0]*g[0]*g[2] + b_ij[0]*g[0]*g[3] + b_ij[0]*g[1]*g[2]);
//        gDot[3] = -(-(b[1]+b_ij[1]+d[1]+s[1])*g[3] + 2*b[1]*g[1]*g[3] + b_ij[1]*g[1]*g[2] + b_ij[1]*g[0]*g[3]);

    }

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



    public void testCorrelations(){

        Double[] b;
        Double[] d = {1.,1.};
        Double[] s;
        Double[] M;// = {3.,3.};

        Double psi;

        Double c1 = 0.01;
        Double c2 = 0.1;

        for (double i = 1.5; i<5; i+=0.25){

            b = new Double[]{i, i};

//            psi = 0.5 * ((i - d[0]) - Math.sqrt((d[0] - i) * (d[0] - i) - .04));  // assume b*s*m=constant


            M = new Double[]{b[0]-d[0]-c2/b[0], b[0]-d[0]-c2/b[0]};     // assume b-d-s=M

            psi = c2/c1 * M[0]; // assume b*m = c1 and b*s = c2
            s = new Double[] {psi,psi};

            FirstOrderIntegrator integrator = new DormandPrince853Integrator(1.0e-4, 1., 1.0e-10, 1.0e-10);//new ClassicalRungeKuttaIntegrator(.01); //

            double T = 10.;
            Boolean augmented = true;

            p0_ODE p_ode = new p0_ODE(b,null, d,s,M, 2, 1, new Double[]{0.});
            p0ge_ODE pg_ode = new p0ge_ODE(b,null, d,s,M, 2, 1, T, new Double[]{0.}, p_ode, Integer.MAX_VALUE,augmented);


            double[] y0 = new double[]{1.,1.,1.,1.};
            double[] y = new double[4];
//
            integrator.integrate(pg_ode, T, y0, 0., y);

            System.out.print("b[0] = "+b[0]+ ", d[0] = " + Math.round(d[0]*100.)/100.+ "\t\t");
            System.out.println(y[0]+"\t"+y[1]+"\t" +y[2]+"\t"+y[3]);
        }
    }

    public static void main(String args[]){

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

        p0_ODE p_ode = new p0_ODE(b,new Double[]{1.,1.}, d,s,M, 2, 1, new Double[]{0.});
        p0ge_ODE pg_ode = new p0ge_ODE(b,new Double[]{1.,1.}, d,s,M, 2, 1, T, new Double[]{0.}, p_ode, Integer.MAX_VALUE,augmented);


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

