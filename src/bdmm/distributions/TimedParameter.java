package bdmm.distributions;

import bdmm.util.Utils;
import beast.core.CalculationNode;
import beast.core.Input;
import beast.core.parameter.RealParameter;

import java.util.Arrays;

/**
 * Parameter in which the value of each element corresponds to a particular time.
 */
public class TimedParameter extends CalculationNode {

    public Input<RealParameter> timesInput = new Input<>(
            "times",
            "Times associated with probabilities.",
            Input.Validate.REQUIRED);

    public Input<Boolean> timesAreAgesInput = new Input<>(
            "timesAreAges",
            "True if times are ages (before most recent sample) instead of times after origin.",
            false);

    public Input<Boolean> timesAreRelativeInput = new Input<>(
            "timesAreRelative",
            "True if times are relative to the origin. (Default false.)",
            false);

    public Input<RealParameter> originInput = new Input<>("origin",
            "Parameter specifying origin of process.");

    public Input<RealParameter> valuesInput = new Input<>(
            "values",
            "Probability values associated with each time.",
            Input.Validate.REQUIRED);

    boolean timesAreAges, timesAreRelative;

    double[] times, storedTimes;
    double[][] values, storedValues;
    double[] valuesAtTime, zeroValuesAtTime;
    int nTimes, nTypes;

    boolean isDirty;

    public TimedParameter() { }

    public TimedParameter(RealParameter timesParam, RealParameter valuesParam) {
        timesInput.setValue(timesParam, this);
        valuesInput.setValue(valuesParam, this);
        initAndValidate();
    }

    @Override
    public void initAndValidate() {
        timesAreAges = timesAreAgesInput.get();
        timesAreRelative = timesAreRelativeInput.get();

        if ((timesAreAges || timesAreRelative) && originInput.get() == null)
            throw new IllegalArgumentException("Origin parameter must be supplied " +
                    "when times are given as ages and/or when times are relative.");

        nTimes = timesInput.get().getDimension();
        times = new double[nTimes];
        storedTimes = new double[nTimes];

        if (valuesInput.get().getDimension() % nTimes != 0)
            throw new IllegalArgumentException("Values parameter must be a " +
                    "multiple of the number of times.");

        nTypes = valuesInput.get().getDimension()/nTimes;

        values = new double[nTimes][nTypes];
        storedValues = new double[nTimes][nTypes];

        valuesAtTime = new double[nTypes];

        zeroValuesAtTime = new double[nTypes];
        Arrays.fill(zeroValuesAtTime, 0.0);

        isDirty = true;
    }

    public int getNTypes() {
        update();

        return nTypes;
    }

    public double[] getTimes() {
        update();

        return times;
    }

    public double[] getValuesAtTime(double time) {
        update();

        int intervalIdx = Arrays.binarySearch(times, time);

        if (intervalIdx<0)
            return zeroValuesAtTime;

        System.arraycopy(values[intervalIdx], 0, valuesAtTime, 0, nTypes);

        return valuesAtTime;
    }

    private void update() {
        if (!isDirty)
            return;

        updateTimes();
        updateValues();

        isDirty = false;
    }

    private void updateTimes() {
        for (int i=0; i<nTimes; i++)
            times[i] = timesInput.get().getValue(i);

        if (timesAreRelative) {
            for (int i=0; i<nTimes; i++)
                times[i] *= originInput.get().getValue();
        }

        if (timesAreAges) {
            Utils.reverseDoubleArray(times);

            for (int i=0; i<times.length; i++) {
                times[i] = originInput.get().getValue()-times[i];
            }
        }
    }

    private void updateValues() {
        for (int timeIdx=0; timeIdx<nTimes; timeIdx++) {
            for (int typeIdx=0; typeIdx<nTypes; typeIdx++) {
                values[timeIdx][typeIdx] = valuesInput.get().getValue(timeIdx*nTypes + typeIdx);
            }
        }

        if (timesAreAges)
            Utils.reverseArray(values);
    }

    @Override
    protected void store() {
        super.store();

        System.arraycopy(times, 0, storedTimes, 0, nTimes);

        for (int timeIdx=0; timeIdx<nTimes; timeIdx++)
            System.arraycopy(values[timeIdx], 0, storedValues[timeIdx], 0, nTypes);
    }

    @Override
    protected void restore() {
        super.restore();
    }

    @Override
    protected boolean requiresRecalculation() {
        isDirty = true;

        return true;
    }
}
