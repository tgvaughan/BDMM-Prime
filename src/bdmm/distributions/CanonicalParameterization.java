package bdmm.distributions;

import beast.core.Input;

public class CanonicalParameterization extends Parameterization {

    public Input<SkylineParameter> migRateInput = new Input<>("migrationRate",
            "Migration rate skyline.", Input.Validate.REQUIRED);

    public Input<SkylineParameter> birthRateInput = new Input<>("birthRate",
            "Birth rate skyline.", Input.Validate.REQUIRED);

    public Input<SkylineParameter> crossBirthRateInput = new Input<>("birthRateAmongDemes",
            "Birth rate among demes skyline.", Input.Validate.REQUIRED);

    public Input<SkylineParameter> deathRateInput = new Input<>("deathRate",
            "Birth rate among demes skyline.", Input.Validate.REQUIRED);

    public Input<SkylineParameter> samplingRateInput = new Input<>("samplingRate",
            "Birth rate among demes skyline.", Input.Validate.REQUIRED);

    public Input<SkylineParameter> removalProbInput = new Input<>("removalRate",
            "Birth rate among demes skyline.", Input.Validate.REQUIRED);

    public Input<SkylineParameter> rhoSamplingInput = new Input<>("rhoSampling",
            "Birth rate among demes skyline.", Input.Validate.REQUIRED);


    @Override
    public double[] getMigRateChangeTimes() {
        return migRateInput.get().getChangeTimes();
    }

    @Override
    public double[] getBirthRateChangeTimes() {
        return birthRateInput.get().getChangeTimes();
    }

    @Override
    public double[] getCrossBirthRateChangeTimes() {
        return crossBirthRateInput.get().getChangeTimes();
    }

    @Override
    public double[] getDeathRateChangeTimes() {
        return deathRateInput.get().getChangeTimes();
    }

    @Override
    public double[] getSamplingRateChangeTimes() {
        return samplingRateInput.get().getChangeTimes();
    }

    @Override
    public double[] getRemovalProbChangeTimes() {
        return removalProbInput.get().getChangeTimes();
    }

    @Override
    public double[] getRhoSamplingTimes() {
        return rhoSamplingInput.get().getChangeTimes();
    }


    @Override
    public int getMigRateChangeCount() {
        return migRateInput.get().getChangeCount();
    }

    @Override
    public int getBirthRateChangeCount() {
        return birthRateInput.get().getChangeCount();
    }

    @Override
    public int getCrossBirthRateChangeCount() {
        return crossBirthRateInput.get().getChangeCount();
    }

    @Override
    public int getDeathRateChangeCount() {
        return deathRateInput.get().getChangeCount();
    }

    @Override
    public int getSamplingRateChangeCount() {
        return samplingRateInput.get().getChangeCount();
    }

    @Override
    public int getRemovalProbChangeCount() {
        return removalProbInput.get().getChangeCount();
    }

    @Override
    public int getRhoChangeCount() {
        return rhoSamplingInput.get().getChangeCount();
    }

    @Override
    public double getMigRateValue(double time) {
        return migRateInput.get().getValueAtTime(time);
    }

    @Override
    public double getBirthRateValue(double time) {
        return birthRateInput.get().getValueAtTime(time);
    }

    @Override
    public double getCrossBirthRateValue(double time) {
        return crossBirthRateInput.get().getValueAtTime(time);
    }

    @Override
    public double getDeathRateValue(double time) {
        return deathRateInput.get().getValueAtTime(time);
    }

    @Override
    public double getSamplingRateValue(double time) {
        return samplingRateInput.get().getValueAtTime(time);
    }

    @Override
    public double getRemovalProbValue(double time) {
        return removalProbInput.get().getValueAtTime(time);
    }

    @Override
    public double getRhoValue(double time) {
        return rhoSamplingInput.get().getValueAtTime(time);
    }
}
