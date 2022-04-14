package bdmmprime.parameterization;

import beast.core.Input;

public class CanonicalParameterization extends Parameterization {

    public Input<SkylineVectorParameter> birthRateInput = new Input<>("birthRate",
            "Birth rate skyline.", Input.Validate.REQUIRED);

    public Input<SkylineVectorParameter> deathRateInput = new Input<>("deathRate",
            "Death rate skyline.", Input.Validate.REQUIRED);

    public Input<SkylineVectorParameter> samplingRateInput = new Input<>("samplingRate",
            "Sampling rate skyline.", Input.Validate.REQUIRED);

    public Input<TimedParameter> rhoSamplingInput = new Input<>("rhoSampling",
            "Contemporaneous sampling times and probabilities.");

    public Input<SkylineVectorParameter> removalProbInput = new Input<>("removalProb",
            "Removal prob skyline.", Input.Validate.REQUIRED);

    public Input<SkylineMatrixParameter> migRateInput = new Input<>("migrationRate",
            "Migration rate skyline.");

    public Input<SkylineMatrixParameter> crossBirthRateInput = new Input<>("birthRateAmongDemes",
            "Birth rate among demes skyline.");

    @Override
    public double[] getMigRateChangeTimes() {
        if (migRateInput.get() == null)
            return EMPTY_TIME_ARRAY;

        return migRateInput.get().getChangeTimes();
    }

    private double[] birthRateChangeTimes;

    @Override
    public double[] getBirthRateChangeTimes() {

        if (crossBirthRateInput.get() != null)
            birthRateChangeTimes = combineAndSortTimes(birthRateChangeTimes,
                    birthRateInput.get().getChangeTimes(),
                    crossBirthRateInput.get().getChangeTimes());
        else
            birthRateChangeTimes = combineAndSortTimes(birthRateChangeTimes,
                    birthRateInput.get().getChangeTimes());

        return birthRateChangeTimes;
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
        if (rhoSamplingInput.get() != null)
            return rhoSamplingInput.get().getTimes();
        else
            return EMPTY_TIME_ARRAY;
    }

    @Override
    protected double[][] getMigRateValues(double time) {
        if (migRateInput.get() == null)
            return ZERO_VALUE_MATRIX;

        return migRateInput.get().getValuesAtTime(time);
    }

    private double[][][] birthRateValues;

    @Override
    protected double[][][] getBirthRateValues(double time) {
        if (birthRateValues == null)
            birthRateValues = new double[nTypes][nTypes][nTypes];

        double[] birthRates = birthRateInput.get().getValuesAtTime(time);

        for (int i = 0; i < nTypes; i++)
            birthRateValues[i][i][i] = birthRates[i];

        if (crossBirthRateInput.get() != null) {
            double[][] crossBirthRates = crossBirthRateInput.get().getValuesAtTime(time);

            for (int i = 0; i < nTypes; i++) {
                for (int j = 0; j < nTypes; j++) {
                    if (i == j)
                        continue;

                    birthRateValues[i][i][j] = crossBirthRates[i][j];
                }
            }
        }

        return birthRateValues;
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
        if (rhoSamplingInput.get() != null)
            return rhoSamplingInput.get().getValuesAtTime(time);
        else
            return ZERO_VALUE_ARRAY;
    }

    @Override
    protected void validateParameterTypeCounts() {
        if (birthRateInput.get().getNTypes() != getNTypes())
            throw new IllegalArgumentException("Birth rate skyline type count does not match type count of model.");

        if (deathRateInput.get().getNTypes() != getNTypes())
            throw new IllegalArgumentException("Death rate skyline type count does not match type count of model.");

        if (samplingRateInput.get().getNTypes() != getNTypes())
            throw new IllegalArgumentException("Sampling rate skyline type count does not match type count of model.");

        if (migRateInput.get() != null
                && migRateInput.get().getNTypes() != getNTypes())
            throw new IllegalArgumentException("Migration rate skyline type count does not match type count of model.");

        if (crossBirthRateInput.get() != null
                && crossBirthRateInput.get().getNTypes() != getNTypes())
            throw new IllegalArgumentException("Birth rate among demes skyline type count does not match type count of model.");

        if (removalProbInput.get().getNTypes() != getNTypes())
            throw new IllegalArgumentException("Removal prob skyline type count does not match type count of model.");

        if (rhoSamplingInput.get() != null
                && rhoSamplingInput.get().getNTypes() != getNTypes())
            throw new IllegalArgumentException("Rho sampling type count does not match type count of model.");
    }
}
