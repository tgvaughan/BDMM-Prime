package bdmmprime.parameterization;

import bdmmprime.util.Utils;
import beast.core.CalculationNode;
import beast.core.Function;
import beast.core.Input;
import beast.core.Loggable;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Tree;

import java.io.PrintStream;
import java.util.Arrays;

/**
 * Parameter in which the value of each element corresponds to a particular time.
 */
public class TimedParameter extends CalculationNode implements Loggable {

    public Input<Function> timesInput = new Input<>(
            "times",
            "Times associated with probabilities.");

    public Input<Boolean> timesAreAgesInput = new Input<>(
            "timesAreAges",
            "True if times are ages (before most recent sample) instead of times after origin.",
            false);

    public Input<Boolean> timesAreRelativeInput = new Input<>(
            "timesAreRelative",
            "True if times are relative to the origin. (Default false.)",
            false);

    public Input<Function> originInput = new Input<>("origin",
            "Parameter specifying origin of process.");

    public Input<Tree> treeInput = new Input<>("tree",
            "Tree when root time is used to identify the start of the process.");

    public Input<TypeSet> typeSetInput = new Input<>("typeSet",
            "Type set defining distinct types in model. Used when a" +
                    "single value is to be shared amongst several types.");


    public Input<Function> valuesInput = new Input<>(
            "values",
            "Probability values associated with each time.");

    boolean timesAreAges, timesAreRelative, inputIsScalar;

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

    public TimedParameter(RealParameter timesParam, RealParameter valuesParam, int nTypes) {
        timesInput.setValue(timesParam, this);
        valuesInput.setValue(valuesParam, this);
        typeSetInput.setValue(new TypeSet(nTypes), this);
        initAndValidate();
    }

    public TimedParameter(RealParameter timesParam, RealParameter valuesParam, RealParameter originParam) {
        timesInput.setValue(timesParam, this);
        valuesInput.setValue(valuesParam, this);
        originInput.setValue(originParam, this);
        timesAreAgesInput.setValue(true, this);
        initAndValidate();
    }

    public TimedParameter(RealParameter timesParam, RealParameter valuesParam, Tree tree) {
        timesInput.setValue(timesParam, this);
        valuesInput.setValue(valuesParam, this);
        treeInput.setValue(tree, this);
        timesAreAgesInput.setValue(true, this);
        initAndValidate();
    }

    @Override
    public void initAndValidate() {
        timesAreAges = timesAreAgesInput.get();
        timesAreRelative = timesAreRelativeInput.get();

        if ((timesAreAges || timesAreRelative) && (originInput.get() == null && treeInput.get() == null))
            throw new IllegalArgumentException("Origin parameter or tree must be supplied " +
                    "when times are given as ages and/or when times are relative.");

        if (originInput.get() != null && treeInput.get() != null)
            throw new IllegalArgumentException("Only one of origin or tree " +
                    "should be specified.");

        if (timesInput.get() != null)
            nTimes = timesInput.get().getDimension();
        else
            nTimes = 0;

        if (nTimes == 0 && valuesInput.get() != null)
            throw new IllegalArgumentException("When timed parameter is empty, " +
                    "both times and values inputs must be omitted.");

        times = new double[nTimes];
        storedTimes = new double[nTimes];

        int valsPerInterval = nTimes>0
                ? valuesInput.get().getDimension() / nTimes
                : 1;

        inputIsScalar = valsPerInterval == 1;

        if (typeSetInput.get() != null) {
            nTypes = typeSetInput.get().getNTypes();

            if (!inputIsScalar && nTypes != valsPerInterval)
                throw new IllegalArgumentException("TimedParameter has an incorrect " +
                        "number of elements.");
        } else {
            nTypes = valsPerInterval;
        }

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

    public int getTimeCount() {
        return times.length;
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
            times[i] = timesInput.get().getArrayValue(i);

        if (timesAreRelative) {
            double startAge = originInput.get() != null
                    ? originInput.get().getArrayValue()
                    : treeInput.get().getRoot().getHeight();

            for (int i=0; i<nTimes; i++)
                times[i] *= startAge;
        }

        if (timesAreAges) {
            Utils.reverseDoubleArray(times);

            double startAge = originInput.get() != null
                    ? originInput.get().getArrayValue()
                    : treeInput.get().getRoot().getHeight();

            for (int i=0; i<times.length; i++) {
                times[i] = startAge-times[i];
            }
        }
    }

    private void updateValues() {
        for (int timeIdx=0; timeIdx<nTimes; timeIdx++) {
            for (int typeIdx=0; typeIdx<nTypes; typeIdx++) {
                if (inputIsScalar)
                    values[timeIdx][typeIdx] = valuesInput.get().getArrayValue(timeIdx);
                else
                    values[timeIdx][typeIdx] = valuesInput.get().getArrayValue(timeIdx*nTypes + typeIdx);
            }
        }

        if (timesAreAges)
            Utils.reverseArray(values);
    }

    @Override
    protected void store() {
        super.store();

        if (nTimes>0)
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

    /*
     * Loggable implementation
     */

    @Override
    public void init(PrintStream out) {

        for (int timeIdx=0; timeIdx<nTimes; timeIdx++) {

            out.print(getID() + "t" + timeIdx + "_time\t");

            if (inputIsScalar) {
                out.print(getID() + "t" + timeIdx);

            } else {
                for (int type = 0; type < nTypes; type++) {
                    out.print(getID() + "t" + timeIdx);

                    if (typeSetInput.get() != null)
                        out.print(typeSetInput.get().getTypeName(type));
                    else
                        out.print("type" + type);

                    out.print("\t");
                }
            }
        }
    }

    @Override
    public void log(long sample, PrintStream out) {

        for (int timeIdx=0; timeIdx<nTimes; timeIdx++) {
            out.print(times[timeIdx] + "\t");

            if (inputIsScalar) {
                out.print(values[timeIdx][0] + "\t");
            } else {
                for (int type = 0; type < nTypes; type++)
                    out.print(values[timeIdx][type] + "\t");
            }
        }
    }

    @Override
    public void close(PrintStream out) {

    }
}
