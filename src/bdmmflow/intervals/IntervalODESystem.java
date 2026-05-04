package bdmmflow.intervals;

import bdmmflow.utils.Result;
import bdmmprime.parameterization.Parameterization;
import bdmmprime.util.Utils;
import org.apache.commons.math3.exception.*;
import org.apache.commons.math3.ode.ContinuousOutputModel;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
import org.apache.commons.math3.ode.nonstiff.*;

import java.util.List;
import java.util.Set;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ForkJoinPool;
import java.util.concurrent.TimeUnit;
import java.util.stream.IntStream;
import java.util.stream.Stream;


/**
 * This class allows to represent ODE systems where the time is split up into intervals.
 * <p>
 * There are two reasons why a time is split into intervals:
 * <p>
 * -    The parameterization of the BDMM model can specify parameterization intervals, where the ODE
 * boundary conditions change at the parameterization interval boundaries.
 * A concrete implementation of this class can inherit the handleParameterizationIntervalBoundary method
 * to specify these changes.
 * <p>
 * -    For numerical stability, it can be useful to split up the time even more fine-grained and to restart
 * integration in each interval.
 */
public abstract class IntervalODESystem implements FirstOrderDifferentialEquations {

    protected List<Interval> intervals;
    protected Parameterization parameterization;

    protected double absoluteTolerance;
    protected double relativeTolerance;
    protected double integrationMinStep;
    protected double integrationMaxStep;

    public IntervalODESystem(Parameterization parameterization, List<Interval> intervals, double absoluteTolerance, double relativeTolerance) {
        this.parameterization = parameterization;
        this.intervals = intervals;
        this.integrationMinStep = this.parameterization.getTotalProcessLength() * 1e-15;
        this.integrationMaxStep = this.parameterization.getTotalProcessLength() / 5;
        this.absoluteTolerance = absoluteTolerance;
        this.relativeTolerance = relativeTolerance;
    }

    /**
     * Integrates over the system forward in time. Integration is restarted at the given intervals.
     * <p>
     * The parameterization interval boundaries are handled automatically
     * by calling handleParameterizationIntervalBoundary at the boundaries.
     *
     * @param initialStates             the initial states at the interval starts-
     * @param intervals                 the list of intervals where integration is restarted. Should include the parameterization intervals.
     *                                  Use IntervalUtils.getIntervals to generate these.
     * @param alwaysStartAtInitialState if the integration should be restarted at the initial state at the interval boundaries.
     *                                  this can increase numerical stability.
     * @return the integration result.
     */
    public ContinuousOutputModel[] integrateForwards(List<double[]> initialStates, List<Interval> intervals, boolean alwaysStartAtInitialState, boolean parallelize) {
        ContinuousOutputModel[] outputModels = new ContinuousOutputModel[intervals.size()];

        if (alwaysStartAtInitialState && parallelize) {

            Stream<Result<Object>> executionResults = IntStream.range(0, intervals.size()).parallel().mapToObj(i -> Result.of(() -> {
                Interval interval = intervals.get(i);
                double[] state = initialStates.get(i).clone();

                this.handleParameterizationIntervalBoundaryIfNecessary(interval.start(), state);
                outputModels[interval.interval()] = this.integrate(state, interval.start(), interval.end(), interval);

                return null;
            }));
            Result.throwIfFailure(executionResults);

        } else {

            double[] state = initialStates.get(0).clone();

            for (int i = 0; i < intervals.size(); i++) {
                Interval interval = intervals.get(i);

                if (alwaysStartAtInitialState) {
                    state = initialStates.get(i).clone();
                }

                this.handleParameterizationIntervalBoundaryIfNecessary(interval.start(), state);
                outputModels[interval.interval()] = this.integrate(state, interval.start(), interval.end(), interval);
            }

        }

        return outputModels;
    }

    /**
     * Integrates over the system backwards in time. Integration is restarted at the given intervals.
     * <p>
     * The parameterization interval boundaries are handled automatically
     * by calling handleParameterizationIntervalBoundary at the boundaries.
     *
     * @param initialStates             the initial states at the interval ends.
     * @param intervals                 the list of intervals where integration is restarted. Should include the parameterization intervals.
     *                                  Use IntervalUtils.getIntervals to generate these.
     * @param alwaysStartAtInitialState if the integration should be restarted at the initial state at the interval boundaries.
     *                                  this can increase numerical stability.
     * @return the integration result.
     */
    public ContinuousOutputModel[] integrateBackwards(List<double[]> initialStates, List<Interval> intervals, boolean alwaysStartAtInitialState, boolean parallelize) {
        ContinuousOutputModel[] outputModels = new ContinuousOutputModel[intervals.size()];

        if (alwaysStartAtInitialState && parallelize) {

            Stream<Result<Object>> executionResults = IntStream.range(0, intervals.size()).parallel().mapToObj(i -> Result.of(() -> {
                Interval interval = intervals.get(i);
                double[] state = initialStates.get(i).clone();

                this.handleParameterizationIntervalBoundaryIfNecessary(interval.end(), state);
                outputModels[intervals.size() - interval.interval() - 1] = this.integrate(state, interval.end(), interval.start(), interval);
                return null;
            }));
            Result.throwIfFailure(executionResults);

        } else {

            double[] state = initialStates.get(initialStates.size() - 1).clone();
            for (int i = intervals.size() - 1; i >= 0; i--) {
                Interval interval = intervals.get(i);

                if (alwaysStartAtInitialState) {
                    state = initialStates.get(i).clone();
                }

                this.handleParameterizationIntervalBoundaryIfNecessary(interval.end(), state);
                outputModels[intervals.size() - interval.interval() - 1] = this.integrate(state, interval.end(), interval.start(), interval);
            }

        }

        return outputModels;
    }

    /**
     * Integrate the system along the given interval from start to end using the given initialState.
     */
    protected ContinuousOutputModel integrate(double[] initialState, double start, double end, Interval interval) {
        try {
            ContinuousOutputModel intervalResult = new ContinuousOutputModel();

            DormandPrince853Integrator integrator = new DormandPrince853Integrator(
                    this.integrationMinStep, this.integrationMaxStep, this.absoluteTolerance, this.relativeTolerance
            );
            integrator.addStepHandler(intervalResult);
            integrator.integrate(this, start, initialState, end, initialState);
            integrator.clearStepHandlers();

            return intervalResult;
        } catch (IllegalStateException e) {
            // NaN was found during integration
            // we switch to the slower but more robust DormandPrince54Integrator
            // with lower relative tolerance and try again

            ContinuousOutputModel intervalResult = new ContinuousOutputModel();

            DormandPrince54Integrator integrator = new DormandPrince54Integrator(
                    this.integrationMinStep, this.integrationMaxStep,
                    this.absoluteTolerance, this.relativeTolerance / 100.0
            );
            integrator.addStepHandler(intervalResult);
            integrator.integrate(this, start, initialState, end, initialState);
            integrator.clearStepHandlers();

            return intervalResult;
        }
    }

    /**
     * This method is called on parameterization boundaries. Inherit this method for custom handling of these boundaries.
     */
    protected void handleParameterizationIntervalBoundary(double boundaryTime, int oldInterval, int newInterval, double[] state) {
        // do nothing
    }

    /**
     * This method is called on parameterization boundaries. Inherit this method for custom handling of these boundaries.
     */
    protected void handleParameterizationIntervalBoundary(double boundaryTime, double[] state) {
        this.handleParameterizationIntervalBoundary(boundaryTime, this.getCurrentParameterizationInterval(boundaryTime), this.getCurrentParameterizationInterval(boundaryTime), state);
    }

    /**
     * Calls handleParameterizationIntervalBoundary if time lies on a parameterization interval boundary.
     */
    protected void handleParameterizationIntervalBoundaryIfNecessary(double time, double[] state) {
        if (this.isParameterizationIntervalBoundary(time)) {
            this.handleParameterizationIntervalBoundary(time, state);
        }
    }

    /**
     * Returns the parameterization interval for the given time.
     */
    public int getCurrentParameterizationInterval(double time) {
        return this.parameterization.getIntervalIndex(time);
    }

    /**
     * Checks if the given time lies on a parameterization interval boundary.
     */
    boolean isParameterizationIntervalBoundary(double time) {
        for (int i = 0; i < this.parameterization.getTotalIntervalCount() - 1; i++) {
            double endTime = this.parameterization.getIntervalEndTimes()[i];
            if (Utils.equalWithPrecision(endTime, time)) return true;
        }
        return false;
    }

}
