package bdmm.distributions;

import beast.core.Input;

public class CanonicalParameterization extends Parameterization {

    public Input<SkylineMatrixParameter> migRateInput = new Input<>("migrationRate",
            "Migration rate skyline.", Input.Validate.REQUIRED);

    public Input<SkylineVectorParameter> birthRateInput = new Input<>("birthRate",
            "Birth rate skyline.", Input.Validate.REQUIRED);

    public Input<SkylineMatrixParameter> crossBirthRateInput = new Input<>("birthRateAmongDemes",
            "Birth rate among demes skyline.", Input.Validate.REQUIRED);

    public Input<SkylineVectorParameter> deathRateInput = new Input<>("deathRate",
            "Birth rate among demes skyline.", Input.Validate.REQUIRED);

    public Input<SkylineVectorParameter> samplingRateInput = new Input<>("samplingRate",
            "Birth rate among demes skyline.", Input.Validate.REQUIRED);

    public Input<SkylineVectorParameter> removalProbInput = new Input<>("removalRate",
            "Birth rate among demes skyline.", Input.Validate.REQUIRED);

    public Input<SkylineScalarParameter> rhoSamplingInput = new Input<>("rhoSampling",
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
    public double[] getMigRateValues(double time) {
        return migRateInput.get().getValuesAtTime(time);
    }

    @Override
    public double[] getBirthRateValues(double time) {
        return birthRateInput.get().getValuesAtTime(time);
    }

    @Override
    public double[] getCrossBirthRateValues(double time) {
        return crossBirthRateInput.get().getValuesAtTime(time);
    }

    @Override
    public double[] getDeathRateValues(double time) {
        return deathRateInput.get().getValuesAtTime(time);
    }

    @Override
    public double[] getSamplingRateValues(double time) {
        return samplingRateInput.get().getValuesAtTime(time);
    }

    @Override
    public double[] getRemovalProbValues(double time) {
        return removalProbInput.get().getValuesAtTime(time);
    }

    @Override
    public double getRhoValue(double time) {
        return rhoSamplingInput.get().getValueAtTime(time);
    }
}
