package bdmmflow.extinctionSystem;

import bdmmflow.intervals.Interval;
import bdmmflow.intervals.IntervalODESystem;
import bdmmprime.parameterization.Parameterization;

import java.util.List;

/**
 * This class represents the ODE system for extinction probabilities.
 */
public class ExtinctionProbabilitiesODESystem extends IntervalODESystem {

    private final double[][] birthRates;
    private final double[][] deathRates;
    private final double[][] samplingRates;
    private final double[][][] crossBirthRates;
    private final double[][][] migrationRates;

    private int currentInterval;

    public ExtinctionProbabilitiesODESystem(Parameterization parameterization, List<Interval> intervals, double absoluteTolerance, double relativeTolerance) {
        super(parameterization, intervals, absoluteTolerance, relativeTolerance);

        this.birthRates = this.parameterization.getBirthRates();
        this.deathRates = this.parameterization.getDeathRates();
        this.samplingRates = this.parameterization.getSamplingRates();
        this.crossBirthRates = this.parameterization.getCrossBirthRates();
        this.migrationRates = this.parameterization.getMigRates();

        this.currentInterval = intervals.size() - 1;
    }

    @Override
    public int getDimension() {
        return this.parameterization.getNTypes();
    }

    @Override
    public void computeDerivatives(double t, double[] y, double[] yDot) {
        if (Double.isNaN(t)) {
            throw new IllegalStateException("NaN detected during integration.");
        }

        for (int i = 0; i < this.parameterization.getNTypes(); i++) {
            yDot[i] = (
                    this.birthRates[this.currentInterval][i]
                            + this.deathRates[this.currentInterval][i]
                            + this.samplingRates[this.currentInterval][i]
                            - this.birthRates[this.currentInterval][i] * y[i]
            )* y[i] - this.deathRates[this.currentInterval][i];

            for (int j = 0; j < this.parameterization.getNTypes(); j++) {
                if (i == j) {
                    continue;
                }

                yDot[i] += (
                        this.crossBirthRates[this.currentInterval][i][j] * (y[i] - y[i] * y[j])
                                + this.migrationRates[this.currentInterval][i][j] * (y[i] - y[j])
                );
            }
        }
    }

    @Override
    protected void handleParameterizationIntervalBoundary(double boundaryTime, int oldInterval, int newInterval, double[] state) {
        super.handleParameterizationIntervalBoundary(boundaryTime, oldInterval, newInterval, state);

        // include rho sampling effects

        for (int type = 0; type < parameterization.getNTypes(); type++) {
            state[type] *= (1.0 - parameterization.getRhoValues()[newInterval][type]);
        }

        this.currentInterval--;
    }

}
