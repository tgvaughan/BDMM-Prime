package bdmm.distributions;

import beast.core.CalculationNode;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Tree;
import beast.evolution.tree.TreeInterface;

import java.util.ArrayList;
import java.util.List;
import java.util.SortedSet;
import java.util.TreeSet;

public abstract class Parameterization extends CalculationNode {

    public Input<RealParameter> originInput = new Input<>("origin",
            "Length of the BDMM process.");

    public Input<Tree> treeInput = new Input<>("treeForMRCA",
            "If provided, condition BDMM on time of MRCA of this tree.",
            Input.Validate.XOR, originInput);

    private boolean dirty;

    SortedSet<Double> modelEventTimesSet = new TreeSet<>();
    List<Double> modelEventTimes = new ArrayList<>();

    @Override
    public void initAndValidate() {
        dirty = true;
    }

    public double getOrigin() {
        return originInput.get().getValue();
    }

    public boolean hasOrigin() {
        return originInput.get() != null;
    }

    public double getMaxTime() {
        if (hasOrigin())
            return getOrigin();
        else
            return treeInput.get().getRoot().getHeight();
    }

    public double getRootEdgeLength(TreeInterface tree) {
        return getMaxTime() - tree.getRoot().getHeight();
    }

    public abstract List<Double> getMigRateChangeTimes(double maxTime);
    public abstract List<Double> getBirthRateChangeTimes(double maxTime);
    public abstract List<Double> getCrossBirthRateChangeTimes(double maxTime);
    public abstract List<Double> getDeathRateChangeTimes(double maxTime);
    public abstract List<Double> getSamplingRateChangeTimes(double maxTime);
    public abstract List<Double> getRemovalProbChangeTimes(double maxTime);
    public abstract List<Double> getRhoSamplingTimes(double maxTime);

    private void update() {
        if (!dirty)
            return;

        updateModelEventTimes();

        dirty = false;
    }

    private void updateModelEventTimes() {

        modelEventTimesSet.clear();
        modelEventTimesSet.addAll(getMigRateChangeTimes(getMaxTime()));
        modelEventTimesSet.addAll(getBirthRateChangeTimes(getMaxTime()));
        modelEventTimesSet.addAll(getCrossBirthRateChangeTimes(getMaxTime()));
        modelEventTimesSet.addAll(getDeathRateChangeTimes(getMaxTime()));
        modelEventTimesSet.addAll(getSamplingRateChangeTimes(getMaxTime()));
        modelEventTimesSet.addAll(getRemovalProbChangeTimes(getMaxTime()));
        modelEventTimesSet.addAll(getRhoSamplingTimes(getMaxTime()));

        modelEventTimes.clear();
        modelEventTimes.addAll(modelEventTimesSet);
    }

    public List<Double> getModelEventTimes() {
        update();

        return modelEventTimes;
    }

    public int getTotalIntervals() {
        update();

        return modelEventTimes.size();
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
