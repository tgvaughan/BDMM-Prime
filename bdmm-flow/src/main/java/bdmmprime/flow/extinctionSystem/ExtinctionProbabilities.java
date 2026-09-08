package bdmmflow.extinctionSystem;

import org.apache.commons.math3.ode.ContinuousOutputModel;

/**
 * This class is a lightweight wrapper of the integration output of ExtinctionProbabilitiesODESystem. It allows
 * to conveniently query the extinction probability at a given time.
 */
public class ExtinctionProbabilities {
    ContinuousOutputModel[] outputModels;
    boolean validateProbabilities = false;
    int n;

    public ExtinctionProbabilities(ContinuousOutputModel[] outputModels, int n) {
        this.outputModels = outputModels;
        this.n = n;
    }

    /**
     * Returns the extinction probability for the given output model at the given time.
     * Note that this method is thread-safe.
     */
    public double[] getProbability(ContinuousOutputModel output, double time) {
        double[] state = new double[this.n];

        synchronized (output) {
            output.setInterpolatedTime(time);
            double[] interpolatedState = output.getInterpolatedState();
            System.arraycopy(interpolatedState, 0, state, 0, this.n);
        }

        if (this.validateProbabilities) {
            // check that all are valid probabilities
            // (interpolation can give rise to values outside [0, 1]

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
            // (interpolation can give rise to values outside [0, 1]

            for (int i = 0; i < state.length; i++) {
                if (state[i] < -0.01 || 1.01 < state[i]) {
                    throw new IllegalStateException("Invalid extinction probability found.");
                }
            }
        }

        return state;
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

}
