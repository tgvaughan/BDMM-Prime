package bdmm.distributions;

import bdmm.util.Utils;
import beast.core.parameter.RealParameter;

public class SkylineVectorParameter extends SkylineParameter {

    int nTypes;

    double[][] values, storedValues;
    double[] valuesAtTime;


    public SkylineVectorParameter() { }

    public SkylineVectorParameter(RealParameter changeTimesParam,
                                  RealParameter rateValuesParam) {
        super(changeTimesParam, rateValuesParam);
    }


    @Override
    public void initAndValidate() {
        super.initAndValidate();

        if (rateValuesInput.get().getDimension() % nIntervals != 0)
            throw new IllegalArgumentException("Value parameter dimension must be a multiple of the number of intervals.");

        nTypes = rateValuesInput.get().getDimension()/nIntervals;

        values = new double[nIntervals][nTypes];
        storedValues = new double[nIntervals][nTypes];

        valuesAtTime = new double[nTypes];
    }

    @Override
    protected void updateValues() {

        for (int interval=0; interval<nIntervals; interval++) {
            for (int i=0; i<nTypes; i++)
                values[interval][i] = rateValuesInput.get().getValue(interval*nTypes + i);
        }

        if (timesAreAges)
            Utils.reverseArray(values);
    }

    protected double[] getValuesAtTime(double time) {
        update();

        int intervalIdx = getIntervalIdx(time);

        System.arraycopy(values[intervalIdx], 0, valuesAtTime, 0, nTypes);

        return valuesAtTime;
    }

    public int getNTypes() {
        return nTypes;
    }

    @Override
    protected void store() {
        super.store();

        for (int interval=0; interval<nIntervals; interval++)
            System.arraycopy(values[interval], 0, storedValues[interval], 0, nTypes);
    }

    @Override
    protected void restore() {
        super.restore();

        double [][] tmp;
        tmp = values;
        values = storedValues;
        storedValues = tmp;
    }
}
