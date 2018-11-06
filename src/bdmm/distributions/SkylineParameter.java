package bdmm.distributions;

import bdmm.util.Utils;
import beast.core.CalculationNode;
import beast.core.Input;
import beast.core.parameter.RealParameter;

import java.util.*;

public class SkylineParameter extends CalculationNode {

    public Input<RealParameter> changeTimesInput = new Input<>("changeTimes",
            "Parameter containing change times for skyline function.");

    public Input<Boolean> timesAreRelativeInput = new Input<>("timesAreRelative",
            "True if times are relative to tree height. (Default false.)",
            false);

    public Input<Boolean> timesAreAgesInput = new Input<>("timesAreAges",
            "True if times are ages (before most recent sample) instead of times after origin.",
            false);

    public Input<RealParameter> rateValuesInput = new Input<>("rateValues",
            "Parameter specifying rate values through time.",
            Input.Validate.REQUIRED);

    public Input<RealParameter> originInput = new Input<>("origin",
            "Parameter specifying origin of process.");

    boolean timesAreAges, timesAreRelative, isScalar;

    double[] times, storedTimes;
    double[] values, storedValues;


    int nIntervals;

    boolean isDirty;

    @Override
    public void initAndValidate() {
        timesAreAges = timesAreAgesInput.get();
        timesAreRelative = timesAreRelativeInput.get();

        if (timesAreAges && originInput.get() == null)
            throw new IllegalArgumentException("Origin parameter must be supplied when times are given as ages.");

        nIntervals = rateValuesInput.get().getDimension();
        int nChangeTimes = changeTimesInput.get() == null ? 0 : changeTimesInput.get().getDimension();

        if (nIntervals != nChangeTimes + 1)
            throw new IllegalArgumentException("Incorrect number of values for given number of change times.");

        times = new double[nIntervals-1];
        values = new double[nIntervals];

        storedTimes = new double[nIntervals-1];
        storedValues = new double[nIntervals];

        isDirty = true;
    }

    public double[] getChangeTimes() {
        update();

        return times;
	}

	public double[] getValues() {
        update();

        return values;
    }

	public int getChangeCount() {
        update();

        return changeTimesInput.get().getDimension();
    }

    public double getValueAtTime(double time) {
        update();

        if (isScalar)
            return values[0];

		int idx = Arrays.binarySearch(times, time);

		if (idx < 0) {
			idx = -idx - 1;
		}

		return values[idx];
    }

	private void update() {
	    if (!isDirty)
	        return;

	    if (isScalar) {
	        values[0] = rateValuesInput.get().getValue();
	        isDirty = false;
	        return;
        }

        for (int i=0; i<nIntervals; i++)
            values[i] = rateValuesInput.get().getValue(i);

        for (int i=0; i<nIntervals-1; i++)
            times[i] = changeTimesInput.get().getValue(i);

        if (timesAreRelative) {
            for (int i=0; i<times.length; i++)
                times[i] *= originInput.get().getValue();
        }

        if (timesAreAges) {
            Utils.reverseDoubleArray(times);
            Utils.reverseDoubleArray(values);

            for (int i=0; i<times.length; i++)
                times[i] = originInput.get().getValue()-times[i];
        }

        isDirty = false;
    }

    @Override
    protected void restore() {
        super.restore();

        double[] tmp;

        if (!isScalar) {
            tmp = times;
            times = storedTimes;
            storedTimes = tmp;
        }

        tmp = values;
        values = storedValues;
        storedValues = tmp;

        isDirty = false;
    }

    @Override
    protected void store() {
        super.store();

        if (!isScalar)
            System.arraycopy(times, 0, storedTimes, 0, times.length);

        System.arraycopy(values, 0, storedValues, 0, values.length);
    }

    @Override
    protected boolean requiresRecalculation() {
	    isDirty = true;
	    return true;
    }
}
