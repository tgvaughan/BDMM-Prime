package bdmm.distributions;

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

    double[] intervalStartTimes;

    double[][] birthRates, deathRates, samplingRates, removalProbs, rhoValues;
    double[][][] migRates, crossBirthRates;

    double[][] storedBirthRates, storedDeathRates, storedSamplingRates,
            storedRemovalProbs, storedRhoValues;
    double[][][] storedMigRates, storedCrossBirthRates;

    int nTypes;

    @Override
    public void initAndValidate() {
        nTypes = nTypesInput.get();

        dirty = true;
    }

    protected abstract double[] getBirthRateChangeTimes();
    protected abstract double[] getMigRateChangeTimes();
    protected abstract double[] getCrossBirthRateChangeTimes();
    protected abstract double[] getDeathRateChangeTimes();
    protected abstract double[] getSamplingRateChangeTimes();
    protected abstract double[] getRemovalProbChangeTimes();
    protected abstract double[] getRhoSamplingTimes();

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
