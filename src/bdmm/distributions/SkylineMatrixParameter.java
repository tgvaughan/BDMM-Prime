package bdmm.distributions;

public class SkylineMatrixParameter extends SkylineParameter {

    int elementsPerMatrix, nTypes;

    double[] valuesAtTime;

    @Override
    public void initAndValidate() {
        super.initAndValidate();

        if (rateValuesInput.get().getDimension() % nIntervals != 0)
            throw new IllegalArgumentException("Value parameter dimension must be a multiple of the number of intervals.");

        elementsPerMatrix = rateValuesInput.get().getDimension()/nIntervals;
        nTypes = (int)Math.round((1 + Math.sqrt(1 + 4*elementsPerMatrix))/2);

        if (elementsPerMatrix != nTypes*(nTypes-1))
            throw new IllegalArgumentException("Wrong number of elements in matrix parameter: shoud be nTypes*(nTypes-1).");

        values = new double[elementsPerMatrix*nIntervals];
        storedValues = new double[elementsPerMatrix*nIntervals];

        valuesAtTime = new double[elementsPerMatrix];
    }

    @Override
    protected void updateValues() {

        for (int interval=0; interval<nIntervals; interval++) {
            int destOffset = interval*elementsPerMatrix;

            int srcOffset = timesAreAges ? (nIntervals-1-interval)*elementsPerMatrix : destOffset;

            for (int i=0; i<elementsPerMatrix; i++)
                values[destOffset + i] = rateValuesInput.get().getValue(srcOffset+i);
        }
    }

    public double[] getValuesAtTime(double time) {
        update();

        int intervalIdx = getIntervalIdx(time);

        System.arraycopy(values, intervalIdx*elementsPerMatrix, valuesAtTime, 0, elementsPerMatrix);

        return valuesAtTime;
    }

    public int getNTypes() {
        return nTypes;
    }
}
