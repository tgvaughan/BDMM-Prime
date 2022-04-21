package bdmmprime.parameterization;

import beast.core.Input;

public class FBDParameterization extends Parameterization {

    public Input<SkylineVectorParameter> diversificationRateInput = new Input<>("diversificationRate",
            "Diversification rate skyline.", Input.Validate.REQUIRED);

    public Input<SkylineVectorParameter> turnoverInput = new Input<>("turnover",
            "Turnover skyline.", Input.Validate.REQUIRED);

    public Input<SkylineVectorParameter> samplingProportionInput = new Input<>("samplingProportion",
            "Sampling proportion skyline.", Input.Validate.REQUIRED);

    public Input<TimedParameter> rhoSamplingInput = new Input<>("rhoSampling",
            "Contemporaneous sampling times and probabilities.");

    public Input<SkylineMatrixParameter> migRateInput = new Input<>("migrationRate",
            "Migration rate skyline.");

    public Input<SkylineMatrixParameter> diversificationRateAmongDemesInput = new Input<>("diversificationRateAmongDemes",
            "Diversification rate among demes skyline.");

    @Override
    public double[] getMigRateChangeTimes() {
        if (migRateInput.get() == null)
            return EMPTY_TIME_ARRAY;

        return migRateInput.get().getChangeTimes();
    }

    private double[] birthRateChangeTimes;

    @Override
    public double[] getBirthRateChangeTimes() {
        if (diversificationRateAmongDemesInput.get() != null) {
            birthRateChangeTimes = combineAndSortTimes(birthRateChangeTimes,
                    diversificationRateInput.get().getChangeTimes(),
                    diversificationRateAmongDemesInput.get().getChangeTimes(),
                    turnoverInput.get().getChangeTimes());
        } else {
            birthRateChangeTimes = combineAndSortTimes(birthRateChangeTimes,
                    diversificationRateInput.get().getChangeTimes(),
                    turnoverInput.get().getChangeTimes());
        }

        return birthRateChangeTimes;
    }

    private double[] deathRateChangeTimes;

    @Override
    public double[] getDeathRateChangeTimes() {
        deathRateChangeTimes = combineAndSortTimes(deathRateChangeTimes,
                turnoverInput.get().getChangeTimes(),
                diversificationRateInput.get().getChangeTimes());

        return deathRateChangeTimes;
    }

    private double[] samplingRateChangeTimes;

    @Override
    public double[] getSamplingRateChangeTimes() {
        samplingRateChangeTimes = combineAndSortTimes(samplingRateChangeTimes,
                samplingProportionInput.get().getChangeTimes(),
                diversificationRateInput.get().getChangeTimes(),
                turnoverInput.get().getChangeTimes());

        return samplingRateChangeTimes;
    }

    @Override
    public double[] getRemovalProbChangeTimes() {
        return EMPTY_TIME_ARRAY;
    }

    @Override
    public double[] getRhoSamplingTimes() {
        if (rhoSamplingInput.get() == null)
            return EMPTY_TIME_ARRAY;

        return rhoSamplingInput.get().getTimes();
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

        double[] divRateVals = diversificationRateInput.get().getValuesAtTime(time);
        double[] toVals = turnoverInput.get().getValuesAtTime(time);

        for (int type = 0; type < nTypes; type++)
            birthRateValues[type][type][type] =
                    divRateVals[type] / (1.0 - toVals[type]);

        if (diversificationRateAmongDemesInput.get() != null) {
            double[][] divRateADVals = diversificationRateAmongDemesInput.get().getValuesAtTime(time);

            for (int sourceType = 0; sourceType < nTypes; sourceType++) {
                for (int destType = 0; destType < sourceType; destType++) {
                    birthRateValues[sourceType][sourceType][destType] =
                            divRateADVals[sourceType][destType]
                                    + divRateVals[sourceType]
                                    * toVals[sourceType] / (1.0 - toVals[sourceType]);
                }
                for (int destType = sourceType+1; destType < nTypes; destType++) {
                    birthRateValues[sourceType][destType][sourceType] =
                            divRateADVals[sourceType][destType]
                                    + divRateVals[sourceType]
                                    * toVals[sourceType] / (1.0 - toVals[sourceType]);
                }
            }
        }

        return birthRateValues;
    }

    @Override
    protected double[] getDeathRateValues(double time) {
        double[] res = diversificationRateInput.get().getValuesAtTime(time);
        double[] toVals = turnoverInput.get().getValuesAtTime(time);

        for (int type=0; type<nTypes; type++)
            res[type] *= toVals[type]/(1.0-toVals[type]);

        return res;
    }

    @Override
    protected double[] getSamplingRateValues(double time) {
        double[] res = diversificationRateInput.get().getValuesAtTime(time);
        double[] sVals = samplingProportionInput.get().getValuesAtTime(time);
        double[] toVals  = turnoverInput.get().getValuesAtTime(time);

        for (int type=0; type<nTypes; type++)
            res[type] *= sVals[type]/(1.0-sVals[type])
                    * toVals[type]/(1.0-toVals[type]);

        return res;
    }

    @Override
    protected double[] getRemovalProbValues(double time) {
        return ZERO_VALUE_ARRAY;
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
        if (diversificationRateInput.get().getNTypes() != getNTypes())
            throw new IllegalArgumentException("Diversification rate skyline type count does not match type count of model.");

        if (turnoverInput.get().getNTypes() != getNTypes())
            throw new IllegalArgumentException("Turnover skyline type count does not match type count of model.");

        if (samplingProportionInput.get().getNTypes() != getNTypes())
            throw new IllegalArgumentException("Sampling proportion skyline type count does not match type count of model.");

        if (migRateInput.get() != null
                && migRateInput.get().getNTypes() != getNTypes())
            throw new IllegalArgumentException("Migration rate skyline type count does not match type count of model.");

        if (diversificationRateAmongDemesInput.get() != null
                && diversificationRateAmongDemesInput.get().getNTypes() != getNTypes())
            throw new IllegalArgumentException("Diversification rate among demes skyline type count does not match type count of model.");

        if (rhoSamplingInput.get() != null
                && rhoSamplingInput.get().getNTypes() != getNTypes())
            throw new IllegalArgumentException("Rho sampling type count does not match type count of model.");
    }
}
