package bdmm.parameterization;

import beast.core.parameter.RealParameter;

public class SkylineMatrixParameter extends SkylineParameter {

    double[][][] values, storedValues;
    double[][] valuesAtTime;

    boolean inputIsScalar;

    public SkylineMatrixParameter() { }

    public SkylineMatrixParameter(RealParameter changeTimesParam,
                                  RealParameter rateValuesParam) {
        super(changeTimesParam, rateValuesParam);
    }

    public SkylineMatrixParameter(RealParameter changeTimesParam,
                                  RealParameter rateValuesParam,
                                  int nTypes) {
        super(changeTimesParam, rateValuesParam, nTypes);
    }

    @Override
    public void initAndValidate() {
        super.initAndValidate();

        int totalElementCount = rateValuesInput.get() != null
                ? rateValuesInput.get().getDimension()
                : 0;

        if (totalElementCount % nIntervals != 0)
            throw new IllegalArgumentException("Value parameter dimension must " +
                    "be a multiple of the number of intervals.");

        int elementsPerMatrix = totalElementCount/nIntervals;
        inputIsScalar = elementsPerMatrix == 1;

        if (typeSetInput.get() != null) {
            nTypes = typeSetInput.get().getNTypes();

            if (!inputIsScalar && elementsPerMatrix != nTypes*(nTypes-1)) {
                throw new IllegalArgumentException("SkylineMatrix parameter has " +
                        "an incorrect number of elements.");
            }
        } else {
            nTypes = (int) Math.round((1 + Math.sqrt(1 + 4 * elementsPerMatrix)) / 2);
        }

        values = new double[nIntervals][nTypes][nTypes];
        storedValues = new double[nIntervals][nTypes][nTypes];

        valuesAtTime = new double[nTypes][nTypes];
    }

    @Override
    protected void updateValues() {
        int idx=0;
        for (int interval=0; interval<nIntervals; interval++) {
            for (int i=0; i<nTypes; i++) {
                for (int j=0; j<nTypes; j++) {
                    if (i==j) {
                        values[interval][i][j] = 0.0;
                        continue;
                    }

                    if (inputIsScalar)
                        values[interval][i][j] = rateValuesInput.get().getValue(interval);
                    else
                        values[interval][i][j] = rateValuesInput.get().getValue(idx);

                    idx += 1;
                }
            }
        }
    }

    /**
     * Retrieve value of matrix parameter at particular time (not age).
     *
     * @param time when to evaluate the parameter.
     * @return the matrix value at the chosen time.
     */
    protected double[][] getValuesAtTime(double time) {
        update();

        int intervalIdx = getIntervalIdx(time);

        for (int i=0; i<nTypes; i++) {
            System.arraycopy(values[intervalIdx][i], 0,
                    valuesAtTime[i], 0, nTypes);
        }

        return valuesAtTime;
    }

    public int getNTypes() {
        return nTypes;
    }

    @Override
    protected void store() {
        super.store();

        for (int interval=0; interval<nIntervals; interval++) {
            for (int i=0; i<nTypes; i++) {
                System.arraycopy(values[interval][i], 0,
                        storedValues[interval][i], 0, nTypes);
            }
        }
    }

    @Override
    protected void restore() {
        super.restore();

        double[][][] tmp;
        tmp = values;
        values = storedValues;
        storedValues = tmp;
    }
}
