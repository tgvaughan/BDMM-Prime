package bdmm.parameterization;

import bdmm.util.Utils;
import beast.core.CalculationNode;
import beast.core.Input;
import beast.core.parameter.RealParameter;

import java.util.*;

public abstract class SkylineParameter extends CalculationNode {

    public Input<RealParameter> changeTimesInput = new Input<>("changeTimes",
            "Parameter containing change times for skyline function.");

    public Input<Boolean> timesAreRelativeInput = new Input<>("timesAreRelative",
            "True if times are relative to origin. (Default false.)",
            false);

    public Input<Boolean> timesAreAgesInput = new Input<>("timesAreAges",
            "True if times are ages (before most recent sample) instead of times after origin.",
            false);

    public Input<RealParameter> rateValuesInput = new Input<>("rateValues",
            "Parameter specifying rate values through time.",
            Input.Validate.REQUIRED);

    public Input<RealParameter> originInput = new Input<>("origin",
            "Parameter specifying origin of process.");

    boolean timesAreAges, timesAreRelative;

    double[] times, storedTimes;


    int nIntervals, nTypes;

    boolean isDirty;

    public SkylineParameter() { }

    public SkylineParameter(RealParameter changeTimesParam,
                            RealParameter rateValuesParam) {
        changeTimesInput.setValue(changeTimesParam, this);
        rateValuesInput.setValue(rateValuesParam, this);
        initAndValidate();
    }

    @Override
    public void initAndValidate() {
        timesAreAges = timesAreAgesInput.get();
        timesAreRelative = timesAreRelativeInput.get();

        if ((timesAreAges || timesAreRelative) && originInput.get() == null)
            throw new IllegalArgumentException("Origin parameter must be supplied " +
                    "when times are given as ages and/or when times are relative.");

        int nChangeTimes = changeTimesInput.get() == null ? 0 : changeTimesInput.get().getDimension();
        nIntervals = nChangeTimes + 1;

        times = new double[nIntervals-1];

        storedTimes = new double[nIntervals-1];

        isDirty = true;
    }

    public int getNTypes() {
        update();

        return nTypes;
    }

    /**
     * @return times (not ages) when parameter changes.
     */
    public double[] getChangeTimes() {
        update();

        return times;
	}

	public int getChangeCount() {
        update();

        return changeTimesInput.get().getDimension();
    }


	protected void update() {
	    if (!isDirty)
	        return;

	    updateTimes();
	    updateValues();

        isDirty = false;
    }

    protected void updateTimes() {

 	    if (nIntervals==1)
	        return;

        for (int i=0; i<nIntervals-1; i++)
            times[i] = changeTimesInput.get().getValue(i);

        if (timesAreRelative) {
            for (int i=0; i<times.length; i++)
                times[i] *= originInput.get().getValue();
        }

        if (timesAreAges) {
            Utils.reverseDoubleArray(times);

            for (int i=0; i<times.length; i++)
                times[i] = originInput.get().getValue()-times[i];
        }
    }

    protected abstract void updateValues();

    protected int getIntervalIdx(double time) {
        if (nIntervals==1)
            return 0;

		int idx = Arrays.binarySearch(times, time);

		if (idx < 0)
			idx = -idx - 1;
		else
		    idx += 1;

		return idx;
    }

    @Override
    protected void restore() {
        super.restore();

        double[] tmp;

        if (nIntervals>1) {
            tmp = times;
            times = storedTimes;
            storedTimes = tmp;
        }

        isDirty = false;
    }

    @Override
    protected void store() {
        super.store();

        if (nIntervals>1)
            System.arraycopy(times, 0, storedTimes, 0, times.length);
    }

    @Override
    protected boolean requiresRecalculation() {
	    isDirty = true;
	    return true;
    }
}
