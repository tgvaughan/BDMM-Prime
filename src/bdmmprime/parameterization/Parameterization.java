package bdmmprime.parameterization;

import bdmmprime.util.Utils;
import beast.core.CalculationNode;
import beast.core.Function;
import beast.core.Input;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;

import java.util.*;

/**
 * Full parameterization for a multi-type birth-death skyline model with sampled ancestors.
 *
 * Objects of this class expose a variety of methods for interrogating
 * the canonical (lambda, mu, M, psi, rho, r, t_origin) parameters at different
 * times.  By "time" we mean a number that increases into the future from
 * some point in the past defined by the start of the birth-death process.
 * (We us "ages" to refer to numbers that increase into the past.)
 *
 * In accordance with birth-death skyline convention, the time period spanned by
 * the birth-death process is broken up into intervals at boundaries t_i.
 * Interval i includes the times (t_i-1,t_i], i.e. it does NOT include the time
 * at the earlier boundary.
 */
public abstract class Parameterization extends CalculationNode {

    public Input<TypeSet> typeSetInput = new Input<>("typeSet",
            "Type set containing types in model.",
            new TypeSet(1));

    public Input<Function> originInput = new Input<>("origin",
            "Time between start of process and the end.");

    public Input<Tree> treeInput = new Input<>("tree",
            "If specified, condition on root time rather than origin time.",
            Input.Validate.XOR, originInput);

    private boolean dirty;

    private SortedSet<Double> intervalEndTimesSet = new TreeSet<>(Utils::precisionLimitedComparator);

    private double[] intervalEndTimes, storedIntervalEndTimes;

    private double[][] deathRates, samplingRates, removalProbs, rhoValues;
    private double[][][] migRates;
    private double[][][][] birthRates;

    private double[][] storedDeathRates, storedSamplingRates,
            storedRemovalProbs, storedRhoValues;
    private double[][][] storedMigRates;
    private double[][][][] storedBirthRates;

    final static double[] EMPTY_TIME_ARRAY = new double[0];
    double[] ZERO_VALUE_ARRAY;
    double[][] ZERO_VALUE_MATRIX;
    double[][][] ZERO_VALUE_ARRAY3D;

    TypeSet typeSet;

    int nTypes;

    @Override
    public void initAndValidate() {
        typeSet = typeSetInput.get();
        nTypes = typeSet.getNTypes();
        ZERO_VALUE_ARRAY = new double[nTypes];
        ZERO_VALUE_MATRIX = new double[nTypes][nTypes];
        ZERO_VALUE_ARRAY3D = new double[nTypes][nTypes][nTypes];

        dirty = true;
        update();
    }

    public abstract double[] getBirthRateChangeTimes();
    public abstract double[] getMigRateChangeTimes();
    public abstract double[] getDeathRateChangeTimes();
    public abstract double[] getSamplingRateChangeTimes();
    public abstract double[] getRemovalProbChangeTimes();
    public abstract double[] getRhoSamplingTimes();

    protected abstract double[][][] getBirthRateValues(double time);
    protected abstract double[][] getMigRateValues(double time);
    protected abstract double[] getDeathRateValues(double time);
    protected abstract double[] getSamplingRateValues(double time);
    protected abstract double[] getRemovalProbValues(double time);
    protected abstract double[] getRhoValues(double time);

    protected abstract void validateParameterTypeCounts();

    public int getNTypes() {
        return nTypes;
    }

    public TypeSet getTypeSet() {
        return typeSet;
    }

    public double getTotalProcessLength() {
        if (originInput.get() != null)
            return originInput.get().getArrayValue();
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

            birthRates = new double[intervalEndTimes.length][nTypes][nTypes][nTypes];
            migRates = new double[intervalEndTimes.length][nTypes][nTypes];
            deathRates = new double[intervalEndTimes.length][nTypes];
            samplingRates = new double[intervalEndTimes.length][nTypes];
            removalProbs = new double[intervalEndTimes.length][nTypes];
            rhoValues = new double[intervalEndTimes.length][nTypes];

            storedBirthRates = new double[intervalEndTimes.length][nTypes][nTypes][nTypes];
            storedMigRates = new double[intervalEndTimes.length][nTypes][nTypes];
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
        addTimes(getDeathRateChangeTimes());
        addTimes(getSamplingRateChangeTimes());
        addTimes(getRemovalProbChangeTimes());
        addTimes(getRhoSamplingTimes());

        intervalEndTimesSet.add(getTotalProcessLength()); // End time of final interval

        if (intervalEndTimes == null) {
            intervalEndTimes = new double[intervalEndTimesSet.size()];
            storedIntervalEndTimes = new double[intervalEndTimesSet.size()];
        }

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

    /**
     * Finds the index of the time interval t lies in.  Note that if t
     * lies on a boundary between intervals, the interval returned will be
     * the _earlier_ of these two intervals.
     *
     * @param t time for which to identify interval
     * @return index identifying interval.
     */
    public int getIntervalIndex(double t) {
        update();

        int index = Arrays.binarySearch(intervalEndTimes, t);

        if (index < 0)
            index = -index - 1;

        // return at most the index of the last interval (m-1)
        return Math.max(0, Math.min(index, intervalEndTimes.length-1));
    }

    void updateValues() {

        for (int interval = 0; interval < intervalEndTimes.length; interval++) {

            double t = intervalEndTimes[interval];
            System.arraycopy(getDeathRateValues(t), 0, deathRates[interval], 0, nTypes);
            System.arraycopy(getSamplingRateValues(t), 0, samplingRates[interval], 0, nTypes);
            System.arraycopy(getRemovalProbValues(t), 0, removalProbs[interval], 0, nTypes);
            System.arraycopy(getRhoValues(t), 0, rhoValues[interval], 0, nTypes);

            double[][] migRateMatrix = getMigRateValues(t);
            double[][][] birthRateArray3D = getBirthRateValues(t);
            for (int i=0; i<nTypes; i++) {
                System.arraycopy(migRateMatrix[i], 0, migRates[interval][i], 0, nTypes);
                for (int j=0; j<nTypes; j++) {
                    System.arraycopy(birthRateArray3D[i][j], 0, birthRates[interval][i][j], 0, nTypes);
                }
            }
        }
    }

    public double[][][][] getBirthRates() {
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
        update();

        return rhoValues;
    }

    public double[][][] getMigRates() {
        update();

        return migRates;
    }

    /**
     * Return time of node, i.e. T - node_age.
     *
     * @param node node whose time to query.
     * @return time of node.
     */
    public double getNodeTime(Node node, double finalSampleOffset) {
        return getTotalProcessLength() - node.getHeight() - finalSampleOffset;
    }

    /**
     * Return the age corresponding to the given time relative to the most recent sample.
     *
     * @param time time to query age for
     * @return age corresponding to time
     */
    public double getAge(double time, double finalSampleOffset) {
        return getTotalProcessLength() - time - finalSampleOffset;
    }

    /**
     * Retrieve index of interval containing node.
     *
     * @param node node whose interval to query.
     * @return index of interval.
     */
    public int getNodeIntervalIndex(Node node, double finalSampleOffset) {
        return getIntervalIndex(getNodeTime(node, finalSampleOffset));
    }


    protected SortedSet<Double> changeTimeSet = new TreeSet<>(Utils::precisionLimitedComparator);

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
       System.arraycopy(intervalEndTimes, 0, storedIntervalEndTimes, 0, intervalEndTimes.length);

        for (int interval=0; interval<intervalEndTimes.length; interval++) {

            System.arraycopy(deathRates[interval], 0, storedDeathRates[interval], 0, nTypes);
            System.arraycopy(samplingRates[interval], 0, storedSamplingRates[interval], 0, nTypes);
            System.arraycopy(removalProbs[interval], 0, storedRemovalProbs[interval], 0, nTypes);
            System.arraycopy(rhoValues[interval], 0, storedRhoValues[interval], 0, nTypes);

            for (int fromType=0; fromType<nTypes; fromType++) {
                System.arraycopy(migRates[interval][fromType], 0,
                        storedMigRates[interval][fromType], 0, nTypes);
                for (int toType1=0; toType1<nTypes; toType1++) {
                    System.arraycopy(birthRates[interval][fromType][toType1], 0,
                            storedBirthRates[interval][fromType][toType1], 0, nTypes);
                }
            }
        }

        super.store();
    }

    @Override
    protected void restore() {

        double[] scalarTmp;
        double[][] vectorTmp;
        double[][][] matrixTmp;
        double[][][][] array3dTmp;

        scalarTmp = intervalEndTimes;
        intervalEndTimes = storedIntervalEndTimes;
        storedIntervalEndTimes = scalarTmp;

        array3dTmp = birthRates;
        birthRates = storedBirthRates;
        storedBirthRates = array3dTmp;

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

        super.restore();
    }
}
