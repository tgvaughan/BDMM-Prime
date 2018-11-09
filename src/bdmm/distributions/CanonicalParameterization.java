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
            "Death rate skyline.", Input.Validate.REQUIRED);

    public Input<SkylineVectorParameter> samplingRateInput = new Input<>("samplingRate",
            "Sampling rate skyline.", Input.Validate.REQUIRED);

    public Input<SkylineVectorParameter> removalProbInput = new Input<>("removalRate",
            "Removal prob skyline.", Input.Validate.REQUIRED);

    public Input<TimedParameter> rhoSamplingInput = new Input<>("rhoSampling",
            "Contemporaneous sampling times and probabilities.", Input.Validate.REQUIRED);

    @Override
    protected double[] getMigRateChangeTimes() {
        return migRateInput.get().getChangeTimes();
    }

    @Override
    protected double[] getBirthRateChangeTimes() {
        return birthRateInput.get().getChangeTimes();
    }

    @Override
    protected double[] getCrossBirthRateChangeTimes() {
        return crossBirthRateInput.get().getChangeTimes();
    }

    @Override
    protected double[] getDeathRateChangeTimes() {
        return deathRateInput.get().getChangeTimes();
    }

    @Override
    protected double[] getSamplingRateChangeTimes() {
        return samplingRateInput.get().getChangeTimes();
    }

    @Override
    protected double[] getRemovalProbChangeTimes() {
        return removalProbInput.get().getChangeTimes();
    }

    @Override
    protected double[] getRhoSamplingTimes() {
        return rhoSamplingInput.get().getTimes();
    }

    @Override
    protected double[][] getMigRateValues(double time) {
        return migRateInput.get().getValuesAtTime(time);
    }

    @Override
    protected double[] getBirthRateValues(double time) {
        return birthRateInput.get().getValuesAtTime(time);
    }

    @Override
    protected double[][] getCrossBirthRateValues(double time) {
        return crossBirthRateInput.get().getValuesAtTime(time);
    }

    @Override
    protected double[] getDeathRateValues(double time) {
        return deathRateInput.get().getValuesAtTime(time);
    }

    @Override
    protected double[] getSamplingRateValues(double time) {
        return samplingRateInput.get().getValuesAtTime(time);
    }

    @Override
    protected double[] getRemovalProbValues(double time) {
        return removalProbInput.get().getValuesAtTime(time);
    }

    @Override
    protected double[] getRhoValues(double time) {
        return rhoSamplingInput.get().getValuesAtTime(time);
    }
}
