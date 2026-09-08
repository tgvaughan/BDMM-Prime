package bdmmflow.flowSystems;

import bdmmflow.extinctionSystem.ExtinctionProbabilities;
import bdmmflow.intervals.Interval;
import bdmmflow.intervals.IntervalODESystem;
import bdmmflow.utils.Utils;
import bdmmprime.parameterization.Parameterization;
import org.apache.commons.math3.linear.*;
import org.apache.commons.math3.ode.ContinuousOutputModel;


import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

/**
 * This class represents the classical backwards-in-time flow ODE.
 */
public class FlowODESystem extends IntervalODESystem implements IFlowODESystem {
    final ExtinctionProbabilities extinctionProbabilities;

    final RealMatrix[] timeInvariantSystemMatrices;

    final double[][] birthRates;
    final double[][] deathRates;
    final double[][] samplingRates;
    final double[][][] crossBirthRates;
    final double[][][] migrationRates;

    int seed;
    double maxConditionNumber;
    boolean useLoucaPennellIntervals;

    public FlowODESystem(
            Parameterization parameterization,
            ExtinctionProbabilities extinctionProbabilities,
            List<Interval> intervals,
            double absoluteTolerance,
            double relativeTolerance,
            int seed,
            double maxConditionNumber,
            boolean useLoucaPennellIntervals) {
        super(parameterization, intervals, absoluteTolerance, relativeTolerance);
        this.extinctionProbabilities = extinctionProbabilities;

        this.birthRates = this.parameterization.getBirthRates();
        this.deathRates = this.parameterization.getDeathRates();
        this.samplingRates = this.parameterization.getSamplingRates();
        this.crossBirthRates = this.parameterization.getCrossBirthRates();
        this.migrationRates = this.parameterization.getMigRates();
        this.useLoucaPennellIntervals = useLoucaPennellIntervals;

        this.seed = seed;
        this.maxConditionNumber = maxConditionNumber;

        this.timeInvariantSystemMatrices = new RealMatrix[this.parameterization.getTotalIntervalCount()];

        for (int i = 0; i < this.parameterization.getTotalIntervalCount(); i++) {
            this.timeInvariantSystemMatrices[i] = this.buildTimeInvariantSystemMatrix(i);
        }
    }

    @Override
    public int getDimension() {
        return parameterization.getNTypes() * parameterization.getNTypes();
    }

    /**
     * Builds the time-invariant part of the system matrix for a given interval. This can be reused.
     */
    RealMatrix buildTimeInvariantSystemMatrix(int interval) {
        RealMatrix system = new BlockRealMatrix(parameterization.getNTypes(), parameterization.getNTypes());

        for (int i = 0; i < parameterization.getNTypes(); i++) {
            system.addToEntry(
                    i,
                    i,
                    this.deathRates[interval][i] + this.samplingRates[interval][i]
            );

            for (int j = 0; j < parameterization.getNTypes(); j++) {
                system.addToEntry(
                        i,
                        i,
                        this.migrationRates[interval][i][j] + this.crossBirthRates[interval][i][j]
                );
                system.addToEntry(
                        i,
                        j,
                        -this.migrationRates[interval][i][j]
                );
            }
        }

        return system;
    }

    /**
     * Builds the time-varying part of the system matrix for a given interval. This has to be computed for every
     * time step.
     */
    void addTimeVaryingSystemMatrix(double t, RealMatrix system) {
        ContinuousOutputModel extinctionOutputModel = this.extinctionProbabilities.getOutputModel(t);

        double[] extinctProbabilities = this.extinctionProbabilities.getProbability(extinctionOutputModel, t);
        int interval = this.getCurrentParameterizationInterval(t);

        for (int i = 0; i < parameterization.getNTypes(); i++) {
            system.addToEntry(
                    i,
                    i,
                    -2 * this.birthRates[interval][i] * extinctProbabilities[i] + this.birthRates[interval][i]
            );

            for (int j = 0; j < parameterization.getNTypes(); j++) {
                system.addToEntry(
                        i,
                        i,
                        -this.crossBirthRates[interval][i][j] * extinctProbabilities[j]
                );

                system.addToEntry(
                        i,
                        j,
                        -this.crossBirthRates[interval][i][j] * extinctProbabilities[i]
                );
            }
        }

    }

    /**
     * Builds the system matrix for a given time point.
     */
    RealMatrix buildSystemMatrix(double t) {
        int interval = this.parameterization.getIntervalIndex(t);
        RealMatrix systemMatrix = this.timeInvariantSystemMatrices[interval].copy();
        this.addTimeVaryingSystemMatrix(t, systemMatrix);
        return systemMatrix;
    }

    @Override
    public void computeDerivatives(double t, double[] y, double[] yDot) {
        if (Double.isNaN(t)) {
            throw new IllegalStateException("NaN detected during integration.");
        }

        int numTypes = this.parameterization.getNTypes();

        RealMatrix yMatrix = Utils.toMatrix(y, numTypes);
        RealMatrix systemMatrix = this.buildSystemMatrix(t);
        RealMatrix yDotMatrix = systemMatrix.multiply(yMatrix);
        Utils.fillArray(yDotMatrix, yDot);
    }

    @Override
    protected void handleParameterizationIntervalBoundary(double boundaryTime, int oldInterval, int newInterval, double[] state) {
        super.handleParameterizationIntervalBoundary(boundaryTime, oldInterval, newInterval, state);

        // include rho sampling effects

        for (int i = 0; i < this.parameterization.getNTypes(); i++) {
            for (int j = 0; j < this.parameterization.getNTypes(); j++) {
                state[i * this.parameterization.getNTypes() + j] *= (1 - this.parameterization.getRhoValues()[newInterval][i]);
            }
        }
    }

    /**
     * Computes the initial states (preconditioners) for the given strategy and intervals.
     */
    List<InitialState> getInitialStates(String initialMatrixStrategy, List<Interval> intervals) {
        return switch (initialMatrixStrategy) {
            case "random" -> {
                RealMatrix matrix = Utils.getRandomMatrix(
                        this.parameterization.getNTypes(), this.seed
                );

                // we condition on having at most a certain condition number for numeric safety

                double maxConditionNumber = this.parameterization.getNTypes() * 3;
                int seedCorrection = 0;
                while (maxConditionNumber < new SingularValueDecomposition(matrix).getConditionNumber()) {
                    matrix = Utils.getRandomMatrix(
                            this.parameterization.getNTypes(), this.seed + ++seedCorrection
                    );
                }

                RealMatrix inverse = MatrixUtils.inverse(matrix);
                List<InitialState> initialStates = new ArrayList<>();
                for (Interval ignored : intervals) {
                    double[] array = new double[this.parameterization.getNTypes() * this.parameterization.getNTypes()];
                    Utils.fillArray(matrix, array);
                    initialStates.add(new InitialState(array, inverse));
                }

                yield initialStates;
            }
            case "identity" -> {
                RealMatrix matrix = MatrixUtils.createRealIdentityMatrix(this.parameterization.getNTypes());

                double[] array = new double[this.parameterization.getNTypes() * this.parameterization.getNTypes()];
                Utils.fillArray(matrix, array);

                List<InitialState> initialStates = new ArrayList<>();
                for (Interval ignored : intervals) {
                    initialStates.add(new InitialState(array, matrix));
                }

                yield initialStates;
            }
            case "average_inverse" -> intervals.stream().parallel().map((interval) -> {
                double h = interval.end() - interval.start();
                RealMatrix startA = this.buildSystemMatrix(interval.start() + bdmmprime.util.Utils.globalPrecisionThreshold);
                RealMatrix midA = this.buildSystemMatrix((interval.start() + interval.end()) / 2.0);
                RealMatrix threeQuarterA = this.buildSystemMatrix(
                        interval.start() + 3.0 * (interval.end() - interval.start()) / 4.0
                );
                RealMatrix endA =  this.buildSystemMatrix(interval.end() - bdmmprime.util.Utils.globalPrecisionThreshold);

                RealMatrix startInvX = Utils.expm(
                        startA.add(midA.scalarMultiply(4)).add(endA).scalarMultiply(h / 6.0)
                );
                RealMatrix midInvX = Utils.expm(
                        midA.add(threeQuarterA.scalarMultiply(4)).add(endA).scalarMultiply(h / 2.0 / 6.0)
                );
                RealMatrix endInvX = MatrixUtils.createRealIdentityMatrix(this.parameterization.getNTypes());

                RealMatrix averageInvX = startInvX.add(midInvX.scalarMultiply(4)).add(endInvX).scalarMultiply(1.0 / 6.0);

                double[] array = new double[this.parameterization.getNTypes() * this.parameterization.getNTypes()];
                Utils.fillArray(averageInvX, array);

                RealMatrix inverse = MatrixUtils.inverse(averageInvX);

                return new InitialState(array, inverse);
            }).toList();
            default -> throw new RuntimeException(
                    "Error: initial state strategy not known."
            );
        };
    }

    /**
     * Calculates the flow integral using the given intervals.
     *
     * @return the calculated flow.
     */
    @Override
    public IFlow calculateFlowIntegral(
            String initialMatrixStrategy,
            boolean parallelize
    ) {
        this.splitUpIntervals();
        boolean resetInitialStateAtIntervalBoundaries = 1 < this.intervals.size();

        List<InitialState> initialStates = this.getInitialStates(initialMatrixStrategy, this.intervals);

        ContinuousOutputModel[] rawOutputs = this.integrateBackwards(
                initialStates.stream().map(InitialState::initialState).toList(),
                this.intervals,
                resetInitialStateAtIntervalBoundaries,
                parallelize
        );

        return new Flow(
                rawOutputs,
                this.parameterization.getNTypes(),
                initialStates,
                resetInitialStateAtIntervalBoundaries
        );
    }

    /**
     * Splits up the stored intervals if numerical issues are expected. Depending on
     * this.useLoucaPennellIntervals, we use their interval heuristic or our own.
     */
    protected void splitUpIntervals() {
        double logMaxConditionNumber = Math.log(this.maxConditionNumber);

        List<Interval> newIntervals = this.intervals.stream().parallel().map((currentOldInterval) -> {
            double currentIntervalEnd = currentOldInterval.end();
            List<Interval> subIntervals = new ArrayList<>();

            while (true) {
                double minNewIntervalStart = currentOldInterval.start();

                RealMatrix currentEndSystemMatrix = this.buildSystemMatrix(currentIntervalEnd - bdmmprime.util.Utils.globalPrecisionThreshold);
                RealMatrix currentMidSystemMatrix = this.buildSystemMatrix((currentIntervalEnd + minNewIntervalStart) / 2);
                RealMatrix currentStartSystemMatrix = this.buildSystemMatrix(minNewIntervalStart + bdmmprime.util.Utils.globalPrecisionThreshold);

                double maxIntervalSize;
                if (this.useLoucaPennellIntervals) {
                    SingularValueDecomposition decomposition = new SingularValueDecomposition(currentEndSystemMatrix);
                    double maxSingularValue = Arrays.stream(decomposition.getSingularValues()).max().orElseThrow();
                    maxIntervalSize = logMaxConditionNumber / (2.0 * maxSingularValue);
                } else {
                    double endSpread = Utils.getHermitianSpread(currentEndSystemMatrix);
                    double midSpread = Utils.getHermitianSpread(currentMidSystemMatrix);
                    double startSpread = Utils.getHermitianSpread(currentStartSystemMatrix);

                    maxIntervalSize = currentIntervalEnd - minNewIntervalStart;
                    maxIntervalSize = Math.min(maxIntervalSize, logMaxConditionNumber / endSpread);
                    maxIntervalSize = Math.min(maxIntervalSize, logMaxConditionNumber / midSpread);
                    maxIntervalSize = Math.min(maxIntervalSize, logMaxConditionNumber / startSpread);
                }

                double newIntervalStart = Math.max(currentIntervalEnd - maxIntervalSize, currentOldInterval.start());

                // find containing parameterization interval ends
                List<Integer> containingParameterizationIntervalEnds = new ArrayList<>();
                for (int j = 0; j < parameterization.getTotalIntervalCount(); j++) {
                    double parameterizationIntervalEndTime = parameterization.getIntervalEndTimes()[j];
                    if (
                            bdmmprime.util.Utils.lessThanWithPrecision(newIntervalStart, parameterizationIntervalEndTime) &&
                                    bdmmprime.util.Utils.lessThanWithPrecision(parameterizationIntervalEndTime, currentIntervalEnd)
                    ) {
                        containingParameterizationIntervalEnds.add(j);
                    }
                }

                Interval newInterval = new Interval(
                        0, currentOldInterval.parameterizationInterval(), newIntervalStart, currentIntervalEnd
                );
                subIntervals.add(0, newInterval);

                if (bdmmprime.util.Utils.equalWithPrecision(newIntervalStart, currentOldInterval.start())) {
                    // we reached the end of the current old interval
                    break;
                }
                currentIntervalEnd = newIntervalStart;
            }

            return subIntervals;
        }).flatMap(List::stream).collect(Collectors.toList());

        // update interval indices

        for (int i = 0; i < newIntervals.size(); i++) {
            Interval interval = newIntervals.get(i);
            Interval intervalWithCorrectIdx = new Interval(
                    i, interval.parameterizationInterval(), interval.start(), interval.end()
            );
            newIntervals.set(i, intervalWithCorrectIdx);
        }
        this.intervals = newIntervals;
    }
}
