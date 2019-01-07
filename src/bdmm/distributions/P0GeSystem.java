package bdmm.distributions;


import bdmm.parameterization.Parameterization;
import bdmm.util.Utils;

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


                gDot[i] += b_ij[interval][i][j]*g[i];
                gDot[i] -= b_ij[interval][i][j]*g[i]*g[j];

                gDot[i] += M[interval][i][j] * g[i];
                gDot[i] -= M[interval][i][j] * g[j];
			}

			/*  ge equations: (dim .. 2*dim-1) */

			gDot[nTypes + i] = + (b[interval][i]+d[interval][i]+s[interval][i]
					- 2*b[interval][i]*g[i])*g[nTypes + i];


			for (int j = 0; j< nTypes; j++){

                if (i==j)
			        continue;

                gDot[nTypes + i] += b_ij[interval][i][j]*g[nTypes + i];
                gDot[nTypes + i] -= b_ij[interval][i][j] *
                        (g[i]*g[nTypes + j] + g[j]*g[nTypes + i]);

                gDot[nTypes + i] += M[interval][i][j] * g[nTypes + i];
                gDot[nTypes + i] -= M[interval][i][j] * g[nTypes + j];
			}
		}
	}
}

