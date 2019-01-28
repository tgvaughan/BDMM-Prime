package bdmm.parameterization;

import bdmm.util.Utils;
import beast.core.parameter.RealParameter;

public class SkylineVectorParameter extends SkylineParameter {

    double[][] values, storedValues;
    double[] valuesAtTime;

    boolean inputIsScalar;

    public SkylineVectorParameter() { }

    public SkylineVectorParameter(RealParameter changeTimesParam,
                                  RealParameter rateValuesParam) {
        super(changeTimesParam, rateValuesParam);
    }

    public SkylineVectorParameter(RealParameter changeTimesParam,
                                  RealParameter rateValuesParam,
                                  int nTypes) {
        super(changeTimesParam, rateValuesParam, nTypes);
    }


    @Override
    public void initAndValidate() {
        super.initAndValidate();

        if (rateValuesInput.get().getDimension() % nIntervals != 0)
            throw new IllegalArgumentException("Value parameter dimension must " +
                    "be a multiple of the number of intervals.");

        int valsPerInterval = rateValuesInput.get().getDimension()/nIntervals;
        inputIsScalar = valsPerInterval==1;

        if (typeSetInput.get() != null) {
            nTypes = typeSetInput.get().getNTypes();

            if (!inputIsScalar && nTypes != valsPerInterval)
                throw new IllegalArgumentException("SkylineVector has an incorrect " +
                        "number of elements.");
        } else {
            nTypes = valsPerInterval;
        }

        values = new double[nIntervals][nTypes];
        storedValues = new double[nIntervals][nTypes];

        valuesAtTime = new double[nTypes];
    }

    @Override
    protected void updateValues() {

        for (int interval=0; interval<nIntervals; interval++) {
            for (int i=0; i<nTypes; i++) {
                if (inputIsScalar)
                    values[interval][i] = rateValuesInput.get().getValue(interval);
                else
                    values[interval][i] = rateValuesInput.get().getValue(interval * nTypes + i);
            }
        }

        if (timesAreAges)
            Utils.reverseArray(values);
    }

    /**
     * Retrieve value of vector at a chosen time (not age).
     *
     * @param time when to evaluate the skyline parameter.
     * @return value of the vector at the chosen time.
     */
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
