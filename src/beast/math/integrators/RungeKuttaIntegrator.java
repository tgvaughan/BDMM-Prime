package beast.math.integrators;

import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;

import java.util.Arrays;

/**
 * @author dkuh004
 *         Date: Jun 26, 2012
 *         Time: 12:07:29 PM
 */
public class RungeKuttaIntegrator {


    double stepsize;
    Boolean forward;


    public RungeKuttaIntegrator(double stepsize){

        this.stepsize = stepsize;
    }

    // Computation by 4th order Runge-Kutta
    public void integrate(FirstOrderDifferentialEquations ODE, double t0, double[] y0, double t, double[] y) throws Exception{

        if (Math.abs(t-t0) < 1e-20 ) y = y0;

        else {

            Boolean forward = (t0 < t);  // make sure the steps go into the right direction
            if ( (forward && stepsize<0) || (!forward && stepsize>0) ) stepsize = -stepsize;

            double k1;
            double k2;
            double k3;
            double k4;
            double[][] yDot = new double [4][ODE.getDimension()];
            double[] yk = new double[ODE.getDimension()];
            double step = stepsize;


            for (double x=t0; (forward && x<t) || (!forward && x>t) ; x+=stepsize){

                if ((forward && x+stepsize>t) || (!forward && x+stepsize<t))
                    step = t-x; // adjust final step size

                ODE.computeDerivatives(x,y,yDot[0]);

                for (int i=0; i<ODE.getDimension(); i++){


                    // Computing all of the trial values
                    k1 = step * yDot[0][i];

                    fill(yk, y, k1/2);
                    ODE.computeDerivatives(x + step/2, yk, yDot[1]);
                    k2 = step * yDot[1][i];

                    fill(yk, y, k2/2);
                    ODE.computeDerivatives(x + step/2, yk, yDot[2]);
                    k3 = step * yDot[2][i];

                    fill(yk, y, k3);
                    ODE.computeDerivatives(x + step, yk, yDot[3]);
                    k4 = step * yDot[3][i];

                    // Incrementing y
                    y[i] += k1/6 + k2/3+ k3/3 + k4/6;
                }
            }
        }
    }

    void fill(double[] yk, double[] y, double k){

        for (int i =0; i < yk.length; i++){

            yk[i] = y[i] + k;

        }

    }

}
