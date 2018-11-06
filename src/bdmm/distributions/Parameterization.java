package bdmm.distributions;

import beast.core.CalculationNode;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Tree;
import beast.evolution.tree.TreeInterface;

import java.util.*;

public abstract class Parameterization extends CalculationNode {

    private boolean dirty;

    SortedSet<Double> modelEventTimesSet = new TreeSet<>();
    List<Double> modelEventTimes = new ArrayList<>();

    double[] birthRates, deathRates, crossBirthRates, samplingRates, removalProbs, rhoValues;
    double[] storedBirthRates, storedDeathRates, storedCrossBirthRates,
            storedSamplingRates, storedRemovalProbs, storedRhoValues;

    @Override
    public void initAndValidate() {
        dirty = true;
    }

    public abstract double[] getMigRateChangeTimes();
    public abstract double[] getBirthRateChangeTimes();
    public abstract double[] getCrossBirthRateChangeTimes();
    public abstract double[] getDeathRateChangeTimes();
    public abstract double[] getSamplingRateChangeTimes();
    public abstract double[] getRemovalProbChangeTimes();
    public abstract double[] getRhoSamplingTimes();

    public abstract double getMigRateValue(double time);
    public abstract double getBirthRateValue(double time);
    public abstract double getCrossBirthRateValue(double time);
    public abstract double getDeathRateValue(double time);
    public abstract double getSamplingRateValue(double time);
    public abstract double getRemovalProbValue(double time);
    public abstract double getRhoValue(double time);

    public abstract int getMigRateChangeCount();
    public abstract int getBirthRateChangeCount();
    public abstract int getCrossBirthRateChangeCount();
    public abstract int getDeathRateChangeCount();
    public abstract int getSamplingRateChangeCount();
    public abstract int getRemovalProbChangeCount();
    public abstract int getRhoChangeCount();

    private void update() {
        if (!dirty)
            return;

        updateModelEventTimes();

        if (birthRates == null) {
            birthRates = new double[modelEventTimes.size()+1];
            crossBirthRates = new double[modelEventTimes.size()+1];
            deathRates = new double[modelEventTimes.size()+1];
            samplingRates = new double[modelEventTimes.size()+1];
            removalProbs = new double[modelEventTimes.size()+1];
            rhoValues = new double[modelEventTimes.size()+1];
        }
        updateBirthDeathPsiParams();

        dirty = false;
    }

    private void addTimes(double[] times) {
        for (double time : times)
            modelEventTimes.add(time);
    }

    private void updateModelEventTimes() {

        modelEventTimesSet.clear();
        addTimes(getMigRateChangeTimes());
        addTimes(getBirthRateChangeTimes());
        addTimes(getCrossBirthRateChangeTimes());
        addTimes(getDeathRateChangeTimes());
        addTimes(getSamplingRateChangeTimes());
        addTimes(getRemovalProbChangeTimes());
        addTimes(getRhoSamplingTimes());

        modelEventTimes.clear();
        modelEventTimes.addAll(modelEventTimesSet);
    }

    public List<Double> getModelEventTimes() {
        update();

        return modelEventTimes;
    }

    public int getTotalModelEvents() {
        update();

        return modelEventTimes.size();
    }

    void updateBirthDeathPsiParams() {

        if (SAModel) {
            removalProbabilities = removalProbability.get().getValues();
            r =  new Double[n*totalIntervals];
        }

        int state;

        for (int i = 0; i < n*totalIntervals; i++) {

            state =  i/totalIntervals;

            birth[i] = (identicalRatesForAllTypes[0]) ? birthRates[index(times[i%totalIntervals], birthRateChangeTimes)] :
                    birthRates[birthRates.length > n ? (birthChanges+1)*state+index(times[i%totalIntervals], birthRateChangeTimes) : state];
            death[i] = (identicalRatesForAllTypes[1]) ? deathRates[index(times[i%totalIntervals], deathRateChangeTimes)] :
                    deathRates[deathRates.length > n ? (deathChanges+1)*state+index(times[i%totalIntervals], deathRateChangeTimes) : state];
            psi[i] = (identicalRatesForAllTypes[2]) ? samplingRates[index(times[i%totalIntervals], samplingRateChangeTimes)] :
                    samplingRates[samplingRates.length > n ? (samplingChanges+1)*state+index(times[i%totalIntervals], samplingRateChangeTimes) : state];
            if (SAModel) r[i] = (identicalRatesForAllTypes[4]) ? removalProbabilities[index(times[i%totalIntervals], rChangeTimes)] :
                    removalProbabilities[removalProbabilities.length > n ? (rChanges+1)*state+index(times[i%totalIntervals], rChangeTimes) : state];

        }

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
