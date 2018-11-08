package bdmm.distributions;

import beast.core.CalculationNode;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Tree;
import beast.evolution.tree.TreeInterface;

import java.util.*;

public abstract class Parameterization extends CalculationNode {

    public Input<Integer> nTypesInput = new Input<>("nTypes",
            "Number of types in model.",
            Input.Validate.REQUIRED);

    private boolean dirty;

    SortedSet<Double> modelEventTimesSet = new TreeSet<>();
    List<Double> modelEventTimes = new ArrayList<>();

    double[] migRates, birthRates, deathRates, crossBirthRates, samplingRates, removalProbs, rhoValues;
    double[] storedMigRates, storedBirthRates, storedDeathRates, storedCrossBirthRates,
            storedSamplingRates, storedRemovalProbs, storedRhoValues;

    int nTypes;

    @Override
    public void initAndValidate() {
        nTypes = nTypesInput.get();

        dirty = true;
    }

    public abstract double[] getMigRateChangeTimes();
    public abstract double[] getBirthRateChangeTimes();
    public abstract double[] getCrossBirthRateChangeTimes();
    public abstract double[] getDeathRateChangeTimes();
    public abstract double[] getSamplingRateChangeTimes();
    public abstract double[] getRemovalProbChangeTimes();
    public abstract double[] getRhoSamplingTimes();

    public abstract double[] getMigRateValues(double time);
    public abstract double[] getBirthRateValues(double time);
    public abstract double[] getCrossBirthRateValues(double time);
    public abstract double[] getDeathRateValues(double time);
    public abstract double[] getSamplingRateValues(double time);
    public abstract double[] getRemovalProbValues(double time);
    public abstract double getRhoValue(double time);

    public abstract int getMigRateChangeCount();
    public abstract int getBirthRateChangeCount();
    public abstract int getCrossBirthRateChangeCount();
    public abstract int getDeathRateChangeCount();
    public abstract int getSamplingRateChangeCount();
    public abstract int getRemovalProbChangeCount();
    public abstract int getRhoChangeCount();

    public int getNTypes() {
        return nTypes;
    }

    private void update() {
        if (!dirty)
            return;

        updateModelEventTimes();

        if (birthRates == null) {

            migRates = new double[nTypes*(nTypes-1)*(modelEventTimes.size()+1)];
            birthRates = new double[nTypes*(modelEventTimes.size()+1)];
            crossBirthRates = new double[nTypes*(nTypes-1)*(modelEventTimes.size()+1)];
            deathRates = new double[nTypes*(modelEventTimes.size()+1)];
            samplingRates = new double[nTypes*(modelEventTimes.size()+1)];
            removalProbs = new double[modelEventTimes.size()+1];
            rhoValues = new double[modelEventTimes.size()+1];

            storedMigRates = new double[nTypes*(nTypes-1)*(modelEventTimes.size()+1)];
            storedBirthRates = new double[modelEventTimes.size()+1];
            storedCrossBirthRates = new double[nTypes*(nTypes-1)*(modelEventTimes.size()+1)];
            storedDeathRates = new double[modelEventTimes.size()+1];
            storedSamplingRates = new double[modelEventTimes.size()+1];
            storedRemovalProbs = new double[modelEventTimes.size()+1];
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

        int state;

        for (int i=0; i<modelEventTimes.size()+1; i++) {

        }
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
