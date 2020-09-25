package bdmmprime.distribution;

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

    @Override
    public String toString() {
	    StringBuilder sb = new StringBuilder();

        for (int type=0; type<dimension; type++) {
            if (type>0)
                sb.append(" ");

            sb.append("p0[").append(type).append("]=").append(p0[type]);
        }

        return sb.toString();
    }
}
