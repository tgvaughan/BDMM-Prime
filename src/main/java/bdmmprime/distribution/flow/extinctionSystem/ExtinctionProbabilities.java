/*
 * Copyright (c) 2017-2026 ETH Zürich
 *
 * This file is part of bdmm-prime.
 *
 * bdmm-prime is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or (at
 * your option) any later version.
 *
 * bdmm-prime is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with bdmm-prime. If not, see <https://www.gnu.org/licenses/>.
 */

package bdmmprime.distribution.flow.extinctionSystem;

import org.apache.commons.math3.ode.ContinuousOutputModel;

/**
 * This class is a lightweight wrapper of the integration output of ExtinctionProbabilitiesODESystem. It allows
 * to conveniently query the extinction probability p0 or the survival probability 1 - p0 at a given time. The
 * survival probability is stored directly, so it keeps its relative accuracy even when p0 is close to 1.
 */
public class ExtinctionProbabilities {

    private final ContinuousOutputModel[] outputModels;
    private final int n;
    private boolean validateProbabilities = false;

    public ExtinctionProbabilities(ContinuousOutputModel[] outputModels, int n) {
        this.outputModels = outputModels;
        this.n = n;
    }

    /**
     * Returns the extinction probability for the given output model at the given time.
     * Note that this method is thread-safe.
     */
    public double[] getProbability(ContinuousOutputModel output, double time) {
        return this.convertToExtinctionProbability(this.getSurvivalProbability(output, time));
    }

    /**
     * Returns the survival probability (one minus the extinction probability) for the given output model at the
     * given time. Note that this method is thread-safe.
     */
    public double[] getSurvivalProbability(ContinuousOutputModel output, double time) {
        double[] state = new double[this.n];

        synchronized (output) {
            output.setInterpolatedTime(time);
            double[] interpolatedState = output.getInterpolatedState();
            System.arraycopy(interpolatedState, 0, state, 0, this.n);
        }

        if (this.validateProbabilities) {
            // check that all are valid probabilities
            // (interpolation can give rise to values outside [0, 1])

            for (int i = 0; i < state.length; i++) {
                if (state[i] < -0.01 || 1.01 < state[i]) {
                    throw new IllegalStateException("Invalid extinction probability found.");
                }
            }
        }

        return state;
    }

    /**
     * Returns the extinction probability for the given output model at the given time.
     * Note that this method is not thread-safe.
     */
    public double[] unsafeGetProbability(ContinuousOutputModel output, double time) {
        output.setInterpolatedTime(time);
        double[] state = output.getInterpolatedState();

        if (this.validateProbabilities) {
            // check that all are valid probabilities
            // (interpolation can give rise to values outside [0, 1])

            for (int i = 0; i < state.length; i++) {
                if (state[i] < -0.01 || 1.01 < state[i]) {
                    throw new IllegalStateException("Invalid extinction probability found.");
                }
            }
        }

        return this.convertToExtinctionProbability(state);
    }

    /**
     * Returns the extinction probability at the given time.
     * Note that this method is thread-safe.
     * @param time the time.
     * @return the extinction probability at the given time.
     */
    public double[] getProbability(double time) {
        return this.getProbability(this.getOutputModel(time), time);
    }

    /**
     * Returns the survival probability (one minus the extinction probability) at the given time.
     * Note that this method is thread-safe.
     * @param time the time.
     * @return the survival probability at the given time.
     */
    public double[] getSurvivalProbability(double time) {
        return this.getSurvivalProbability(this.getOutputModel(time), time);
    }

    /**
     * Returns the continuous output model for the given time.
     */
    public ContinuousOutputModel getOutputModel(double time) {
        if (time > this.outputModels[0].getInitialTime()) {
            return this.outputModels[0];
        }

        for (ContinuousOutputModel model : this.outputModels) {
            if (model.getInitialTime() >= time && time > model.getFinalTime()) {
                return model;
            }
        }

        return this.outputModels[this.outputModels.length - 1];
    }

    /**
     * Enables or disables that an error is thrown when getProbability does not return
     * a valid probability (outside [0, 1]).
     */
    public void validateProbabilities(boolean validateProbabilities) {
        this.validateProbabilities = validateProbabilities;
    }

    /**
     * Returns a new array with the extinction probabilities that correspond to the given survival probabilities.
     */
    private double[] convertToExtinctionProbability(double[] survivalProbability) {
        double[] extinctionProbability = new double[this.n];

        for (int i = 0; i < this.n; i++) {
            extinctionProbability[i] = 1.0 - survivalProbability[i];
        }

        return extinctionProbability;
    }

}
