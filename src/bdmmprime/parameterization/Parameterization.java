/*
 * Copyright (C) 2019-2024 Tim Vaughan
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
import beast.base.evolution.tree.Node;
import beast.base.inference.CalculationNode;

import java.util.*;

/**
 * Full parameterization for a multi-type birth-death skyline model with sampled ancestors.
 *
 * Objects of this class expose a variety of methods for interrogating
 * the canonical (lambda, mu, M, psi, rho, r, t_origin) parameters at different
 * times.  By "time" we mean a number that increases into the future from
 * some point in the past defined by the start of the birth-death process.
 * (We use "ages" to refer to numbers that increase into the past.)
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

    public Input<Function> processLengthInput = new Input<>("processLength",
            "Time between start of process and the end.",
            Input.Validate.REQUIRED);

    private boolean dirty;

    private SortedSet<Double> intervalEndTimesSet = new TreeSet<>(Utils::precisionLimitedComparator);

    private double[] intervalEndTimes, storedIntervalEndTimes;

    private double[][] birthRates, deathRates, samplingRates, removalProbs, rhoValues;
    private double[][][] migRates, crossBirthRates2;
    private double[][][][] crossBirthRates3;

    private double[][] storedBirthRates, storedDeathRates, storedSamplingRates,
            storedRemovalProbs, storedRhoValues;
    private double[][][] storedMigRates, storedCrossBirthRates2;
    private double[][][][] storedCrossBirthRates3;

    final static double[] EMPTY_TIME_ARRAY = new double[0];
    double[] ZERO_VALUE_ARRAY;
    double[][] ZERO_VALUE_ARRAY2;
    double[][][] ZERO_VALUE_ARRAY3;

    TypeSet typeSet;

    int nTypes;

    @Override
    public void initAndValidate() {
        typeSet = typeSetInput.get();
        nTypes = typeSet.getNTypes();
        ZERO_VALUE_ARRAY = new double[nTypes];
        ZERO_VALUE_ARRAY2 = new double[nTypes][nTypes];
        ZERO_VALUE_ARRAY3 = new double[nTypes][nTypes][nTypes];

        dirty = true;
        update();
    }

    public abstract double[] getBirthRateChangeTimes();
    public abstract double[] getMigRateChangeTimes();
    public abstract double[] getCrossBirthRate2ChangeTimes();
    public abstract double[] getCrossBirthRate3ChangeTimes();
    public abstract double[] getDeathRateChangeTimes();
    public abstract double[] getSamplingRateChangeTimes();
    public abstract double[] getRemovalProbChangeTimes();
    public abstract double[] getRhoSamplingTimes();

    protected abstract double[] getBirthRateValues(double time);
    protected abstract double[][] getMigRateValues(double time);
    protected abstract double[][] getCrossBirthRate2Values(double time);
    protected abstract double[][][] getCrossBirthRate3Values(double time);
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
        return processLengthInput.get().getArrayValue();
    }

    private void update() {
        if (!dirty)
            return;

        updateModelEventTimes();

        if (birthRates == null) {
            validateParameterTypeCounts();

            birthRates = new double[intervalEndTimes.length][nTypes];
            migRates = new double[intervalEndTimes.length][nTypes][nTypes];
            crossBirthRates2 = new double[intervalEndTimes.length][nTypes][nTypes];
            crossBirthRates3 = new double[intervalEndTimes.length][nTypes][nTypes][nTypes];
            deathRates = new double[intervalEndTimes.length][nTypes];
            samplingRates = new double[intervalEndTimes.length][nTypes];
            removalProbs = new double[intervalEndTimes.length][nTypes];
            rhoValues = new double[intervalEndTimes.length][nTypes];

            storedBirthRates = new double[intervalEndTimes.length][nTypes];
            storedMigRates = new double[intervalEndTimes.length][nTypes][nTypes];
            storedCrossBirthRates2 = new double[intervalEndTimes.length][nTypes][nTypes];
            storedCrossBirthRates3 = new double[intervalEndTimes.length][nTypes][nTypes][nTypes];
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


        addTimes(getBirthRateChangeTimes());
        addTimes(getDeathRateChangeTimes());
        addTimes(getSamplingRateChangeTimes());
        addTimes(getRemovalProbChangeTimes());
        addTimes(getRhoSamplingTimes());
        if (nTypes > 1) {
            addTimes(getMigRateChangeTimes());
            addTimes(getCrossBirthRate2ChangeTimes());
            addTimes(getCrossBirthRate3ChangeTimes());
        }

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
            System.arraycopy(getBirthRateValues(t), 0, birthRates[interval], 0, nTypes);
            System.arraycopy(getDeathRateValues(t), 0, deathRates[interval], 0, nTypes);
            System.arraycopy(getSamplingRateValues(t), 0, samplingRates[interval], 0, nTypes);
            System.arraycopy(getRemovalProbValues(t), 0, removalProbs[interval], 0, nTypes);
            System.arraycopy(getRhoValues(t), 0, rhoValues[interval], 0, nTypes);

            double[][] migRateMatrix = getMigRateValues(t);
            double[][] crossBirthRate2Matrix = getCrossBirthRate2Values(t);
            double[][][] crossBirthRate3Matrix = getCrossBirthRate3Values(t);

            if (nTypes>1) {
                for (int i = 0; i < nTypes; i++) {
                    System.arraycopy(migRateMatrix[i], 0, migRates[interval][i], 0, nTypes);
                    System.arraycopy(crossBirthRate2Matrix[i], 0, crossBirthRates2[interval][i], 0, nTypes);

                    for (int j=0; j < nTypes; j++) {
                        System.arraycopy(crossBirthRate3Matrix[i][j], 0, crossBirthRates3[interval][i][j], 0, nTypes);
                    }
                }
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
        update();

        return rhoValues;
    }

    public double[][][] getMigRates() {
        update();

        return migRates;
    }

    public double[][][] getCrossBirthRates2() {
        update();

        return crossBirthRates2;
    }

    public double[][][][] getCrossBirthRates3() {
        update();

        return crossBirthRates3;
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
    public double getNodeAge(double time, double finalSampleOffset) {
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

            System.arraycopy(birthRates[interval], 0, storedBirthRates[interval], 0, nTypes);
            System.arraycopy(deathRates[interval], 0, storedDeathRates[interval], 0, nTypes);
            System.arraycopy(samplingRates[interval], 0, storedSamplingRates[interval], 0, nTypes);
            System.arraycopy(removalProbs[interval], 0, storedRemovalProbs[interval], 0, nTypes);
            System.arraycopy(rhoValues[interval], 0, storedRhoValues[interval], 0, nTypes);

            if (nTypes>1) {
                for (int fromType = 0; fromType < nTypes; fromType++) {
                    System.arraycopy(migRates[interval][fromType], 0,
                            storedMigRates[interval][fromType], 0, nTypes);
                    System.arraycopy(crossBirthRates2[interval][fromType], 0,
                            storedCrossBirthRates2[interval][fromType], 0, nTypes);

                    for (int toType1 = 0; toType1 < nTypes; toType1++) {
                        System.arraycopy(crossBirthRates3[interval][fromType][toType1], 0,
                                storedCrossBirthRates3[interval][fromType][toType1], 0, nTypes);
                    }
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

        scalarTmp = intervalEndTimes;
        intervalEndTimes = storedIntervalEndTimes;
        storedIntervalEndTimes = scalarTmp;

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

        if (nTypes>1) {
            matrixTmp = migRates;
            migRates = storedMigRates;
            storedMigRates = matrixTmp;

            matrixTmp = crossBirthRates2;
            crossBirthRates2 = storedCrossBirthRates2;
            storedCrossBirthRates2 = matrixTmp;
        }

        super.restore();
    }
}
