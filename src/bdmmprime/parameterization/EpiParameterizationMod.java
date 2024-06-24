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

public class EpiParameterizationMod extends Parameterization {

    public Input<SkylineVectorParameter> R0Input = new Input<>("R0",
            "Basic reproduction number skyline.", Input.Validate.REQUIRED);

    public Input<SkylineVectorParameter> R0modInput = new Input<>("R0mod",
            "Modulation for R0.", Input.Validate.REQUIRED);

    public Input<SkylineVectorParameter> becomeUninfectiousRateInput = new Input<>("becomeUninfectiousRate",
            "Become uninfectious rate skyline.", Input.Validate.REQUIRED);

    public Input<SkylineVectorParameter> samplingProportionInput = new Input<>("samplingProportion",
            "Sampling proportion skyline.", Input.Validate.REQUIRED);

    public Input<TimedParameter> rhoSamplingInput = new Input<>("rhoSampling",
            "Contemporaneous sampling times and probabilities.");

    public Input<SkylineVectorParameter> removalProbInput = new Input<>("removalProb",
            "Removal prob skyline.", Input.Validate.REQUIRED);

    public Input<SkylineMatrixParameter> migRateInput = new Input<>("migrationRate",
            "Migration rate skyline.");

    public Input<SkylineMatrixParameter> R0AmongDemesInput = new Input<>("R0AmongDemes",
            "Basic reproduction number among demes skyline.");

    @Override
    public double[] getMigRateChangeTimes() {
        if (migRateInput.get() == null)
            return EMPTY_TIME_ARRAY;

        return migRateInput.get().getChangeTimes();
    }

    private double[] birthRateChangeTimes;

    @Override
    public double[] getBirthRateChangeTimes() {
        birthRateChangeTimes = combineAndSortTimes(birthRateChangeTimes,
                R0Input.get().getChangeTimes(),
                R0modInput.get().getChangeTimes(),
                becomeUninfectiousRateInput.get().getChangeTimes());

        return birthRateChangeTimes;
    }


    private double[] crossBirthRateChangeTimes;

    @Override
    public double[] getCrossBirthRateChangeTimes() {
        if (R0AmongDemesInput.get() == null)
            return EMPTY_TIME_ARRAY;

        crossBirthRateChangeTimes = combineAndSortTimes(crossBirthRateChangeTimes,
                R0AmongDemesInput.get().getChangeTimes(),
                R0modInput.get().getChangeTimes(),
                becomeUninfectiousRateInput.get().getChangeTimes());

        return crossBirthRateChangeTimes;
    }


    private double[] deathRateChangeTimes;

    @Override
    public double[] getDeathRateChangeTimes() {
        deathRateChangeTimes = combineAndSortTimes(deathRateChangeTimes,
                becomeUninfectiousRateInput.get().getChangeTimes(),
                samplingProportionInput.get().getChangeTimes(),
                removalProbInput.get().getChangeTimes());

        return deathRateChangeTimes;
    }

    private double[] samplingRateChangeTimes;

    @Override
    public double[] getSamplingRateChangeTimes() {
        samplingRateChangeTimes = combineAndSortTimes(samplingRateChangeTimes,
                samplingProportionInput.get().getChangeTimes(),
                becomeUninfectiousRateInput.get().getChangeTimes());

        return samplingRateChangeTimes;
    }

    @Override
    public double[] getRemovalProbChangeTimes() {
        return removalProbInput.get().getChangeTimes();
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

    @Override
    protected double[] getBirthRateValues(double time) {
        double[] res = R0Input.get().getValuesAtTime(time);
        double[] mod = R0modInput.get().getValuesAtTime(time);
        double[] buVals = becomeUninfectiousRateInput.get().getValuesAtTime(time);

        for (int type=0; type<nTypes; type++)
            res[type] *= mod[type]*buVals[type];

        return res;
    }

    @Override
    protected double[][] getCrossBirthRateValues(double time) {
        if (R0AmongDemesInput.get() == null)
            return ZERO_VALUE_MATRIX;

        double[][] res = R0AmongDemesInput.get().getValuesAtTime(time);
        double[] mod = R0modInput.get().getValuesAtTime(time);
        double[] buVals = becomeUninfectiousRateInput.get().getValuesAtTime(time);

        for (int sourceType=0; sourceType<nTypes; sourceType++) {
            for (int destType=0; destType<nTypes; destType++) {
                if (sourceType==destType)
                    continue;

                res[sourceType][destType] *= mod[sourceType]*buVals[sourceType];
            }
        }

        return res;
    }

    @Override
    protected double[] getDeathRateValues(double time) {
        double[] res = becomeUninfectiousRateInput.get().getValuesAtTime(time);
        double[] samplingProp = samplingProportionInput.get().getValuesAtTime(time);
        double[] removalProb = removalProbInput.get().getValuesAtTime(time);

        for (int type=0; type<nTypes; type++)
            res[type] *= (1 - samplingProp[type])
                    / (1.0 - (1.0-removalProb[type])*samplingProp[type]);

        return res;
    }

    @Override
    protected double[] getSamplingRateValues(double time) {
        double[] res = samplingProportionInput.get().getValuesAtTime(time);
        double[] buRate  = becomeUninfectiousRateInput.get().getValuesAtTime(time);
        double[] removalProb  = removalProbInput.get().getValuesAtTime(time);

        for (int type=0; type<nTypes; type++)
            res[type] = res[type]*buRate[type]/(1 - (1-removalProb[type])*res[type]);

        return res;
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
        if (R0Input.get().getNTypes() != getNTypes())
            throw new IllegalArgumentException("R0 skyline type count ("
                    + R0Input.get().getNTypes()
                    + ") does not match type count of model (" +
                    + getNTypes() + ").");

        if (R0modInput.get().getNTypes() != getNTypes())
            throw new IllegalArgumentException("R0mod skyline type count does not match type count of model.");

        if (becomeUninfectiousRateInput.get().getNTypes() != getNTypes())
            throw new IllegalArgumentException("Become uninfectious rate skyline type count does not match type count of model.");

        if (samplingProportionInput.get().getNTypes() != getNTypes())
            throw new IllegalArgumentException("Sampling proportion skyline type count does not match type count of model.");

        if (migRateInput.get() != null
                && migRateInput.get().getNTypes() != getNTypes())
            throw new IllegalArgumentException("Migration rate skyline type count does not match type count of model.");

        if (R0AmongDemesInput.get() != null
                && R0AmongDemesInput.get().getNTypes() != getNTypes())
            throw new IllegalArgumentException("R0 among demes skyline type count does not match type count of model.");

        if (removalProbInput.get().getNTypes() != getNTypes())
            throw new IllegalArgumentException("Removal prob skyline type count does not match type count of model.");

        if (rhoSamplingInput.get() != null
                && rhoSamplingInput.get().getNTypes() != getNTypes())
            throw new IllegalArgumentException("Rho sampling type count does not match type count of model.");
    }
}
