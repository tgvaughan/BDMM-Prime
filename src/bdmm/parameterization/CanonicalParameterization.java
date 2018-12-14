package bdmm.parameterization;

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

    public Input<SkylineVectorParameter> removalProbInput = new Input<>("removalProb",
            "Removal prob skyline.", Input.Validate.REQUIRED);

    public Input<TimedParameter> rhoSamplingInput = new Input<>("rhoSampling",
            "Contemporaneous sampling times and probabilities.", Input.Validate.REQUIRED);

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

    @Override
    protected void validateParameterTypeCounts() {
        if (birthRateInput.get().getNTypes() != getNTypes())
            throw new IllegalArgumentException("Birth rate skyline type count does not match type count of model.");

        if (deathRateInput.get().getNTypes() != getNTypes())
            throw new IllegalArgumentException("Death rate skyline type count does not match type count of model.");

        if (samplingRateInput.get().getNTypes() != getNTypes())
            throw new IllegalArgumentException("Sampling rate skyline type count does not match type count of model.");

        if (migRateInput.get().getNTypes() != getNTypes())
            throw new IllegalArgumentException("Migration rate skyline type count does not match type count of model.");

        if (crossBirthRateInput.get().getNTypes() != getNTypes())
            throw new IllegalArgumentException("Birth rate among demes skyline type count does not match type count of model.");

        if (removalProbInput.get().getNTypes() != getNTypes())
            throw new IllegalArgumentException("Removal prob skyline type count does not match type count of model.");

        if (rhoSamplingInput.get().getNTypes() != getNTypes())
            throw new IllegalArgumentException("Rho sampling type count does not match type count of model.");
    }
}
