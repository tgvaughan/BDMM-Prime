package bdmm.distributions;

/**
 * Class containing the values of P0.
 */
public class P0State {

    int dimension;
    public double[] p0;

    public P0State(int nTypes) {
        dimension = nTypes;
        p0 = new double[nTypes];
    }

    public P0State(double[] p0) {
        dimension = p0.length;
        this.p0 = p0;
    }

    public P0State() {
        dimension = 1;
        p0 = new double[] {0};
    }
}
