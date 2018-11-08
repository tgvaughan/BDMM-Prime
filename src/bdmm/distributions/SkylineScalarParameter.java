package bdmm.distributions;

import bdmm.util.Utils;

import java.util.Arrays;

public class SkylineScalarParameter extends  SkylineParameter {

    @Override
    public void initAndValidate() {
        super.initAndValidate();

        if (rateValuesInput.get().getDimension() != nIntervals)
            throw new IllegalArgumentException("Incorrect number of values for given number of change times.");

        values = new double[nIntervals];
        storedValues = new double[nIntervals];
    }

    @Override
    protected void updateValues() {

        for (int i=0; i<nIntervals; i++)
            values[i] = rateValuesInput.get().getValue(i);

        if (nIntervals==1)
            return;

        if (timesAreAges) {
            Utils.reverseDoubleArray(values);
        }
    }

    public double getValueAtTime(double time) {
        update();

        return values[getIntervalIdx(time)];
    }
}
