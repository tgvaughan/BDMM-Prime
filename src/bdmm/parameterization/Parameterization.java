package bdmm.parameterization;

import beast.core.CalculationNode;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Tree;

import java.util.*;

public abstract class Parameterization extends CalculationNode {

    public Input<Integer> nTypesInput = new Input<>("nTypes",
            "Number of types in model.",
            Input.Validate.REQUIRED);

    public Input<RealParameter> originInput = new Input<>("origin",
            "Time between start of process and most recent sample.");

    public Input<Tree> treeInput = new Input<>("tree",
            "If specified, condition on root time rather than origin time.",
            Input.Validate.XOR, originInput);

    private boolean dirty;

    private SortedSet<Double> intervalEndTimesSet = new TreeSet<>();

    private double[] intervalEndTimes;

    private double[][] birthRates, deathRates, samplingRates, removalProbs, rhoValues;
    private double[][][] migRates, crossBirthRates;

    private double[][] storedBirthRates, storedDeathRates, storedSamplingRates,
            storedRemovalProbs, storedRhoValues;
    private double[][][] storedMigRates, storedCrossBirthRates;

    final static double[] EMPTY_TIME_ARRAY = new double[0];
    double[] ZERO_VALUE_ARRAY;
    double[][] ZERO_VALUE_MATRIX;

    int nTypes;

    @Override
    public void initAndValidate() {
        nTypes = nTypesInput.get();
        ZERO_VALUE_ARRAY = new double[nTypes];
        ZERO_VALUE_MATRIX = new double[nTypes][nTypes];

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

    public double getTotalProcessLength() {
        if (originInput.get() != null)
            return originInput.get().getValue();
        else
            return treeInput.get().getRoot().getHeight();
    }

    public boolean conditionedOnRoot() {
        return originInput.get() == null;
    }

    private void update() {
        if (!dirty)
            return;

        updateModelEventTimes();

        if (birthRates == null) {
            validateParameterTypeCounts();

            birthRates = new double[intervalEndTimes.length][nTypes];
            migRates = new double[intervalEndTimes.length][nTypes][nTypes];
            crossBirthRates = new double[intervalEndTimes.length][nTypes][nTypes];
            deathRates = new double[intervalEndTimes.length][nTypes];
            samplingRates = new double[intervalEndTimes.length][nTypes];
            removalProbs = new double[intervalEndTimes.length][nTypes];
            rhoValues = new double[intervalEndTimes.length][nTypes];

            storedBirthRates = new double[intervalEndTimes.length][nTypes];
            storedMigRates = new double[intervalEndTimes.length][nTypes][nTypes];
            storedCrossBirthRates = new double[intervalEndTimes.length][nTypes][nTypes];
            storedDeathRates = new double[intervalEndTimes.length][nTypes];
            storedSamplingRates = new double[intervalEndTimes.length][nTypes];
            storedRemovalProbs = new double[intervalEndTimes.length][nTypes];
            storedRhoValues = new double[intervalEndTimes.length][nTypes];
        }

        updateValues();

        dirty = false;
    }

    private void addTimes(double[] times) {
        if (times == null)
            return;

        for (double time : times)
            intervalEndTimesSet.add(time);
    }

    private void updateModelEventTimes() {

        intervalEndTimesSet.clear();

        addTimes(getMigRateChangeTimes());
        addTimes(getBirthRateChangeTimes());
        addTimes(getCrossBirthRateChangeTimes());
        addTimes(getDeathRateChangeTimes());
        addTimes(getSamplingRateChangeTimes());
        addTimes(getRemovalProbChangeTimes());
        addTimes(getRhoSamplingTimes());

        intervalEndTimesSet.add(getTotalProcessLength()); // End time of final interval

        if (intervalEndTimes == null)
            intervalEndTimes = new double[intervalEndTimesSet.size()];

        List<Double> timeList = new ArrayList<>(intervalEndTimesSet);
        for (int i = 0; i< intervalEndTimesSet.size(); i++)
            intervalEndTimes[i] = timeList.get(i);
    }

    public double[] getIntervalEndTimes() {
        update();

        return intervalEndTimes;
    }

    public int getTotalIntervalCount() {
        update();

        return intervalEndTimes.length;
    }

    void updateValues() {

        for (int interval = 0; interval < intervalEndTimes.length; interval++) {

            double t = intervalEndTimes[interval];
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
    protected void store() {
        for (int interval=0; interval<intervalEndTimes.length; interval++) {

            System.arraycopy(birthRates[interval], 0, storedBirthRates[interval], 0, nTypes);
            System.arraycopy(deathRates[interval], 0, storedDeathRates[interval], 0, nTypes);
            System.arraycopy(samplingRates[interval], 0, storedSamplingRates[interval], 0, nTypes);
            System.arraycopy(removalProbs[interval], 0, storedRemovalProbs[interval], 0, nTypes);
            System.arraycopy(rhoValues[interval], 0, storedRhoValues[interval], 0, nTypes);

            for (int fromType=0; fromType<nTypes; fromType++) {
                System.arraycopy(migRates[interval][fromType], 0, storedMigRates[interval][fromType], 0, nTypes);
                System.arraycopy(crossBirthRates[interval][fromType], 0, storedCrossBirthRates[interval][fromType], 0, nTypes);
            }
        }

        super.store();
    }

    @Override
    protected void restore() {
        double[][] vectorTmp;
        double[][][] matrixTmp;

        vectorTmp = birthRates;
        birthRates = storedBirthRates;
        storedBirthRates = vectorTmp;

        vectorTmp = deathRates;
        deathRates = storedDeathRates;
        storedDeathRates = vectorTmp;

        vectorTmp = samplingRates;
        samplingRates = storedSamplingRates;
        storedSamplingRates = vectorTmp;

        vectorTmp = removalProbs;
        removalProbs = storedRemovalProbs;
        storedRemovalProbs = vectorTmp;

        vectorTmp = rhoValues;
        rhoValues = storedRhoValues;
        storedRhoValues = vectorTmp;

        matrixTmp = migRates;
        migRates = storedMigRates;
        storedMigRates = matrixTmp;

        matrixTmp = crossBirthRates;
        crossBirthRates = storedCrossBirthRates;
        storedCrossBirthRates = matrixTmp;

        super.restore();
    }
}
