/*
 * Copyright (C) 2019-2024 ETH Zurich
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

    public Input<SkylineMatrixParameter> birthRateAmongDemesInput = new Input<>("birthRateAmongDemes",
            "Birth rate among demes skyline.");

    @Override
    public double[] getMigRateChangeTimes() {
        if (migRateInput.get() == null)
            return EMPTY_TIME_ARRAY;

        return migRateInput.get().getChangeTimes();
    }

    @Override
    public double[] getBirthRateChangeTimes() {
        return birthRateInput.get().getChangeTimes();
    }

    @Override
    public double[] getCrossBirthRate2ChangeTimes() {
        if (birthRateAmongDemesInput.get() == null)
            return EMPTY_TIME_ARRAY;

        return birthRateAmongDemesInput.get().getChangeTimes();
    }

    @Override
    public double[] getCrossBirthRate3ChangeTimes() {
        return EMPTY_TIME_ARRAY;
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
            return ZERO_VALUE_ARRAY2;

        return migRateInput.get().getValuesAtTime(time);
    }

    @Override
    protected double[] getBirthRateValues(double time) {
        return birthRateInput.get().getValuesAtTime(time);
    }

    @Override
    protected double[][] getCrossBirthRate2Values(double time) {
        if (birthRateAmongDemesInput.get() == null)
            return ZERO_VALUE_ARRAY2;

        return birthRateAmongDemesInput.get().getValuesAtTime(time);
    }

    @Override
    protected double[][][] getCrossBirthRate3Values(double time) {
        return null;
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

        if (birthRateAmongDemesInput.get() != null
                && birthRateAmongDemesInput.get().getNTypes() != getNTypes())
            throw new IllegalArgumentException("Birth rate among demes skyline type count does not match type count of model.");

        if (removalProbInput.get().getNTypes() != getNTypes())
            throw new IllegalArgumentException("Removal prob skyline type count does not match type count of model.");

        if (rhoSamplingInput.get() != null
                && rhoSamplingInput.get().getNTypes() != getNTypes())
            throw new IllegalArgumentException("Rho sampling type count does not match type count of model.");
    }
}
