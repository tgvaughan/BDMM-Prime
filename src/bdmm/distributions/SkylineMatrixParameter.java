package bdmm.distributions;

import beast.core.parameter.RealParameter;

import java.util.Arrays;

public class SkylineMatrixParameter extends SkylineParameter {

    int elementsPerMatrix;


    double[][][] values, storedValues;
    double[][] valuesAtTime;

    public SkylineMatrixParameter() { }

    public SkylineMatrixParameter(RealParameter changeTimesParam,
                                  RealParameter rateValuesParam) {
        super(changeTimesParam, rateValuesParam);
    }

    @Override
    public void initAndValidate() {
        super.initAndValidate();

        if (rateValuesInput.get().getDimension() % nIntervals != 0)
            throw new IllegalArgumentException("Value parameter dimension must be a multiple of the number of intervals.");

        elementsPerMatrix = rateValuesInput.get().getDimension()/nIntervals;
        nTypes = (int)Math.round((1 + Math.sqrt(1 + 4*elementsPerMatrix))/2);

        if (elementsPerMatrix != nTypes*(nTypes-1))
            throw new IllegalArgumentException("Wrong number of elements in matrix parameter: shoud be nTypes*(nTypes-1).");

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

                    values[interval][i][j] = rateValuesInput.get().getValue(idx);

                    idx += 1;
                }
            }
        }
    }

    public double[][] getValuesAtTime(double time) {
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
