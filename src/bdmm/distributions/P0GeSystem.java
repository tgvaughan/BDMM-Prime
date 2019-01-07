package bdmm.distributions;


import bdmm.parameterization.Parameterization;

/**
 * User: Denise
 * Date: Jul 11, 2013
 * Time: 5:56:21 PM
 */

public class P0GeSystem extends P0System {

	public P0GeSystem(Parameterization parameterization) {
	    super(parameterization);
	}

	@Override
	public int getDimension() {
		return 2*this.nTypes;
	}

	@Override
	public void computeDerivatives(double t, double[] y, double[] yDot) {

		for (int i = 0; i<nTypes; i++){

			/*  p0 equations (0 .. dim-1) */

			yDot[i] = + (b[interval][i]+d[interval][i]+s[interval][i]
					- b[interval][i] * y[i]) * y[i]
					- d[interval][i] ;

			for (int j = 0; j < nTypes; j++){

			    if (i==j)
			        continue;


                yDot[i] += b_ij[interval][i][j]*y[i];
                yDot[i] -= b_ij[interval][i][j]*y[i]*y[j];

                yDot[i] += M[interval][i][j] * y[i];
                yDot[i] -= M[interval][i][j] * y[j];
			}

			/*  ge equations: (dim .. 2*dim-1) */

			yDot[nTypes + i] = + (b[interval][i]+d[interval][i]+s[interval][i]
					- 2*b[interval][i]*y[i])*y[nTypes + i];


			for (int j = 0; j< nTypes; j++){

                if (i==j)
			        continue;

                yDot[nTypes + i] += b_ij[interval][i][j]*y[nTypes + i];
                yDot[nTypes + i] -= b_ij[interval][i][j] *
                        (y[i]*y[nTypes + j] + y[j]*y[nTypes + i]);

                yDot[nTypes + i] += M[interval][i][j] * y[nTypes + i];
                yDot[nTypes + i] -= M[interval][i][j] * y[nTypes + j];
			}
		}
	}
}

