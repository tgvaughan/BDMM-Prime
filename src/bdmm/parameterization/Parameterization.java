package bdmm.parameterization;

import beast.core.CalculationNode;
import beast.core.Input;
import beast.core.parameter.RealParameter;

import java.util.*;

public abstract class Parameterization extends CalculationNode {

    public Input<Integer> nTypesInput = new Input<>("nTypes",
            "Number of types in model.",
            Input.Validate.REQUIRED);

    public Input<RealParameter> originInput = new Input<>("origin",
            "Time between start of process and most recent sample.",
            Input.Validate.REQUIRED);

    private boolean dirty;

    SortedSet<Double> intervalStartTimesSet = new TreeSet<>();

    private double[] intervalStartTimes;

    private double[][] birthRates, deathRates, samplingRates, removalProbs, rhoValues;
    private double[][][] migRates, crossBirthRates;

    private double[][] storedBirthRates, storedDeathRates, storedSamplingRates,
            storedRemovalProbs, storedRhoValues;
    private double[][][] storedMigRates, storedCrossBirthRates;

    final static double[] EMPTY_TIME_ARRAY = new double[0];
    double[] ZERO_VALUE_ARRAY;

    int nTypes;

    @Override
    public void initAndValidate() {
        nTypes = nTypesInput.get();
        ZERO_VALUE_ARRAY = new double[nTypes];

        dirty = true;
    }

    public abstract double[] getBirthRateChangeTimes();
    public abstract double[] getMigRateChangeTimes();
    public abstract double[] getCrossBirthRateChangeTimes();
    public abstract double[] getDeathRateChangeTimes();
    public abstract double[] getSamplingRateChangeTimes();
    public abstract double[] getRemovalProbChangeTimes();
    public abstract double[] getRhoSamplingTimes();

    protected abstract double[] getBirthRateValues(double time);
    protected abstract double[][] getMigRateValues(double time);
    protected abstract double[][] getCrossBirthRateValues(double time);
    protected abstract double[] getDeathRateValues(double time);
    protected abstract double[] getSamplingRateValues(double time);
    protected abstract double[] getRemovalProbValues(double time);
    protected abstract double[] getRhoValues(double time);

    protected abstract void validateParameterTypeCounts();

    public int getNTypes() {
        return nTypes;
    }

    public double getOrigin() {
        return originInput.get().getValue();
    }

    private void update() {
        if (!dirty)
            return;

        updateModelEventTimes();

        if (birthRates == null) {
            validateParameterTypeCounts();

            birthRates = new double[intervalStartTimes.length][nTypes];
            migRates = new double[intervalStartTimes.length][nTypes][nTypes];
            crossBirthRates = new double[intervalStartTimes.length][nTypes][nTypes];
            deathRates = new double[intervalStartTimes.length][nTypes];
            samplingRates = new double[intervalStartTimes.length][nTypes];
            removalProbs = new double[intervalStartTimes.length][nTypes];
            rhoValues = new double[intervalStartTimes.length][nTypes];

            storedBirthRates = new double[intervalStartTimes.length][nTypes];
            storedMigRates = new double[intervalStartTimes.length][nTypes][nTypes];
            storedCrossBirthRates = new double[intervalStartTimes.length][nTypes][nTypes];
            storedDeathRates = new double[intervalStartTimes.length][nTypes];
            storedSamplingRates = new double[intervalStartTimes.length][nTypes];
            storedRemovalProbs = new double[intervalStartTimes.length][nTypes];
            storedRhoValues = new double[intervalStartTimes.length][nTypes];
        }

        updateValues();

        dirty = false;
    }

    private void addTimes(double[] times) {
        if (times == null)
            return;

        for (double time : times)
            intervalStartTimesSet.add(time);
    }

    private void updateModelEventTimes() {

        intervalStartTimesSet.clear();

        intervalStartTimesSet.add(0.0); // Start time of first interval

        addTimes(getMigRateChangeTimes());
        addTimes(getBirthRateChangeTimes());
        addTimes(getCrossBirthRateChangeTimes());
        addTimes(getDeathRateChangeTimes());
        addTimes(getSamplingRateChangeTimes());
        addTimes(getRemovalProbChangeTimes());
        addTimes(getRhoSamplingTimes());

        if (intervalStartTimes == null)
            intervalStartTimes = new double[intervalStartTimesSet.size()];

        List<Double> timeList = new ArrayList<>(intervalStartTimesSet);
        for (int i=0; i<intervalStartTimesSet.size(); i++)
            intervalStartTimes[i] = timeList.get(i);
    }

    public double[] getIntervalStartTimes() {
        update();

        return intervalStartTimes;
    }

    public int getTotalIntervalCount() {
        update();

        return intervalStartTimes.length;
    }

    void updateValues() {

        for (int interval = 0; interval< intervalStartTimes.length; interval++) {

            double t = intervalStartTimes[interval];
            System.arraycopy(getBirthRateValues(t), 0, birthRates[interval], 0, nTypes);
            System.arraycopy(getDeathRateValues(t), 0, deathRates[interval], 0, nTypes);
            System.arraycopy(getSamplingRateValues(t), 0, samplingRates[interval], 0, nTypes);
            System.arraycopy(getRemovalProbValues(t), 0, removalProbs[interval], 0, nTypes);
            System.arraycopy(getRhoValues(t), 0, rhoValues[interval], 0, nTypes);

            double[][] migRateMatrix = getMigRateValues(t);
            double[][] crossBirthRateMatrix = getCrossBirthRateValues(t);
            for (int i=0; i<nTypes; i++) {
                System.arraycopy(migRateMatrix[i], 0, migRates[interval][i], 0, nTypes);
                System.arraycopy(crossBirthRateMatrix[i], 0, crossBirthRates[interval][i], 0, nTypes);
            }
        }
    }

    public double[][] getBirthRates() {
        update();

        return birthRates;
    }

    public double[][] getDeathRates() {
        update();

        return deathRates;
    }

    public double[][] getSamplingRates() {
        update();

        return samplingRates;
    }

    public double[][] getRemovalProbs() {
        update();

        return removalProbs;
    }

    public double[][] getRhoValues() {
        updateValues();

        return rhoValues;
    }

    public double[][][] getMigRates() {
        update();

        return migRates;
    }

    public double[][][] getCrossBirthRates() {
        update();

        return crossBirthRates;
    }


    protected SortedSet<Double> changeTimeSet = new TreeSet<>();

    /**
     * Combine times from individual time arrays, removing duplicates.
     *
     * @param changeTimeArrays One or more arrays to combine.
     * @return combined time array
     */
    protected double[] combineAndSortTimes(double[] destArray, double[] ... changeTimeArrays) {
        changeTimeSet.clear();

        for (double[] changeTimeArray : changeTimeArrays) {
            for (double t : changeTimeArray)
                changeTimeSet.add(t);
        }

        if (destArray == null || destArray.length != changeTimeSet.size())
            destArray = new double[changeTimeSet.size()];

        int i=0;
        for (double changeTime : changeTimeSet)
            destArray[i++] = changeTime;

        return destArray;
    }

    @Override
    protected boolean requiresRecalculation() {
        dirty = true;
        return true;
    }

    @Override
    protected void restore() {
        dirty = true;
        super.restore();
    }
}
