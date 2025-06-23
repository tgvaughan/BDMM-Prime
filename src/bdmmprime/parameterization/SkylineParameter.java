/*
 * Copyright (C) 2019-2024 ETH Zurich
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package bdmmprime.parameterization;

import bdmmprime.util.Utils;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.inference.CalculationNode;

import java.io.PrintStream;
import java.util.*;

public abstract class SkylineParameter extends CalculationNode implements Loggable {

    public Input<Function> changeTimesInput = new Input<>("changeTimes",
            "Parameter containing change times for skyline function.");

    public Input<Boolean> timesAreRelativeInput = new Input<>("timesAreRelative",
            "True if times are relative to process start. (Default false.)",
            false);

    public Input<Boolean> timesAreAgesInput = new Input<>("timesAreAges",
            "True if times are ages (before the end of the sampling process) instead " +
                    "of times after birth-death process start.",
            false);

    public Input<Function> skylineValuesInput = new Input<>("skylineValues",
            "Parameter specifying parameter values through time.",
            Input.Validate.REQUIRED);

    public Input<Function> processLengthInput = new Input<>("processLength",
            "Time between start of process and the end.");

    public Input<TypeSet> typeSetInput = new Input<>("typeSet",
            "Type set defining distinct types in model. Used when a" +
                    "single value is to be shared amongst several types.");

    public Input<Boolean> isScalarInput = new Input<>("isScalar",
            "BEAUti hint to regard parameter as scalar. WARNING: " +
                    "the value of this input has no impact outside of BEAUti!",
            true);

    public Input<Boolean> linkIdenticalValuesInput = new Input<>("linkIdenticalValues",
            "BEAUti hint to create XMLs in which identical values are considered " +
                    "linked. WARNING: The value of this input usually has no impact " +
                    "outside of BEAUti!",
            false);

    boolean timesAreAges, timesAreRelative;

    double[] times, storedTimes;


    int nIntervals, nTypes;

    boolean isDirty;

    public SkylineParameter() { }

    public SkylineParameter(Function changeTimesParam,
                            Function skylineValuesParam) {
        changeTimesInput.setValue(changeTimesParam, this);
        skylineValuesInput.setValue(skylineValuesParam, this);
        initAndValidate();
    }

    public SkylineParameter(Function changeTimesParam,
                            Function skylineValuesParam,
                            int nTypes,
                            Function processLength) {

        changeTimesInput.setValue(changeTimesParam, this);
        skylineValuesInput.setValue(skylineValuesParam, this);
        typeSetInput.setValue(new TypeSet(nTypes), this);

        if (processLength != null) {
            this.timesAreAgesInput.setValue(true, this);
            this.processLengthInput.setValue(processLength, this);
        }

        initAndValidate();
    }

    @Override
    public void initAndValidate() {
        timesAreAges = timesAreAgesInput.get();
        timesAreRelative = timesAreRelativeInput.get();

        if ((timesAreAges || timesAreRelative) && processLengthInput.get() == null)
            throw new IllegalArgumentException("Process length parameter or tree must be supplied " +
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

        return times.length;
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
            times[i] = changeTimesInput.get().getArrayValue(i);

        if (timesAreRelative) {

            double startAge = processLengthInput.get().getArrayValue();

            for (int i=0; i<times.length; i++)
                times[i] *= startAge;
        }

        if (timesAreAges) {
            Utils.reverseDoubleArray(times);

            double startAge = processLengthInput.get().getArrayValue();

            for (int i=0; i<times.length; i++)
                times[i] = startAge-times[i];
        }
    }

    protected abstract void updateValues();

    /**
     * Retrieve index of interval containing time.  Note that
     * index returned by the the time at a boundary between two intervals
     * is the index of the earlier interval.
     *
     * @param time time at which to determine index
     * @return index
     */
    protected int getIntervalIdx(double time) {
        if (nIntervals==1)
            return 0;

		int idx = Arrays.binarySearch(times, time);

		if (idx < 0)
			idx = -idx - 1;

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

    /**
     * Beauti hack
     */

    public boolean epochVisualizerDisplayed = false;
}
