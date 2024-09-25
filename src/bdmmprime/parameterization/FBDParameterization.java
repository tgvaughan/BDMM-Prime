/*
 * Copyright (C) 2019-2024 Tim Vaughan
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package bdmmprime.parameterization;

import beast.base.core.Input;

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
        if (diversificationRateAmongDemesInput.get() == null){
            birthRateChangeTimes = combineAndSortTimes(birthRateChangeTimes,
                turnoverInput.get().getChangeTimes());
        }
        else {
            birthRateChangeTimes = combineAndSortTimes(birthRateChangeTimes,
                diversificationRateAmongDemesInput.get().getChangeTimes(),
                turnoverInput.get().getChangeTimes());
        }

        return birthRateChangeTimes;
    }


    private double[] crossBirthRateChangeTimes;

    @Override
    public double[] getCrossBirthRate2ChangeTimes() {
        if (diversificationRateAmongDemesInput.get() == null)
            return EMPTY_TIME_ARRAY;

        crossBirthRateChangeTimes = combineAndSortTimes(crossBirthRateChangeTimes,
                diversificationRateInput.get().getChangeTimes(),
                diversificationRateAmongDemesInput.get().getChangeTimes(),
                turnoverInput.get().getChangeTimes());

        return crossBirthRateChangeTimes;
    }

    @Override
    public double[] getCrossBirthRate3ChangeTimes() {
        return EMPTY_TIME_ARRAY;
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
            return ZERO_VALUE_ARRAY2;

        return migRateInput.get().getValuesAtTime(time);
    }

    @Override
    protected double[] getBirthRateValues(double time) {
        double[] res = diversificationRateInput.get().getValuesAtTime(time);
        double[] toVals = turnoverInput.get().getValuesAtTime(time);

        for (int type=0; type<nTypes; type++)
            res[type] /= 1.0 - toVals[type];

        return res;
    }

    @Override
    protected double[][] getCrossBirthRate2Values(double time) {
        if (diversificationRateAmongDemesInput.get() == null)
            return ZERO_VALUE_ARRAY2;

        double[][] res = diversificationRateAmongDemesInput.get().getValuesAtTime(time);
        double[] dVals = diversificationRateInput.get().getValuesAtTime(time);
        double[] toVals = turnoverInput.get().getValuesAtTime(time);

        for (int sourceType=0; sourceType<nTypes; sourceType++) {
            for (int destType=0; destType<nTypes; destType++) {
                if (sourceType==destType)
                    continue;

                res[sourceType][destType] += dVals[sourceType]
                        * toVals[sourceType]/(1.0-toVals[sourceType]);
            }
        }

        return res;
    }

    @Override
    protected double[][][] getCrossBirthRate3Values(double time) {
        return ZERO_VALUE_ARRAY3;
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
