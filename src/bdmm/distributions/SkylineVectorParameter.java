package bdmm.distributions;

public class SkylineVectorParameter extends SkylineParameter {

    int nTypes;

    double[] valuesAtTime;

    @Override
    public void initAndValidate() {
        super.initAndValidate();

        if (rateValuesInput.get().getDimension() % nIntervals != 0)
            throw new IllegalArgumentException("Value parameter dimension must be a multiple of the number of intervals.");

        nTypes = rateValuesInput.get().getDimension()/nIntervals;

        values = new double[nTypes*nIntervals];
        storedValues = new double[nTypes*nIntervals];

        valuesAtTime = new double[nTypes];
    }

    @Override
    protected void updateValues() {

        for (int interval=0; interval<nIntervals; interval++) {
            int destOffset = interval*nTypes;

            int srcOffset = timesAreAges ? (nIntervals-1-interval)*nTypes : destOffset;

            for (int i=0; i<nTypes; i++)
                values[destOffset + i] = rateValuesInput.get().getValue(srcOffset+i);
        }
    }

    public double[] getValuesAtTime(double time) {
        update();

        int intervalIdx = getIntervalIdx(time);

        System.arraycopy(values, intervalIdx*nTypes, valuesAtTime, 0, nTypes);

        return valuesAtTime;
    }

    public int getNTypes() {
        return nTypes;
    }
}
