package bdmm.distributions;

/**
 * Created by Jeremie Scire (jscire)
 */

/**
 * Class containing the values of P0 and Ge.
 */
public class P0GeState {
	
	int dimension;
	public SmallNumber[] ge;
	public double[] p0;

	public P0GeState(int nTypes) {
	    dimension = nTypes;
        p0 = new double[nTypes];
		ge = new SmallNumber[nTypes];
		for (int i = 0; i<nTypes; i++)
		    ge[i] = new SmallNumber();

    }
	
	public P0GeState(double[] p0, SmallNumber[] ge) {
		if(p0.length != ge.length) {
			throw new RuntimeException("Incorrect initialization: difference of size between ge and p0");
		}
		dimension = p0.length;
		this.p0 = p0;
		this.ge = ge;
	}
	
	public P0GeState() {
		dimension = 1;
		p0 = new double[] {0};
		ge = new SmallNumber[] {new SmallNumber()};
	}
}
