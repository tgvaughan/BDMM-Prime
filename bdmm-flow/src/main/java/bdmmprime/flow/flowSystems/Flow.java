package bdmmflow.flowSystems;

import bdmmflow.utils.Utils;
import org.apache.commons.math3.linear.*;
import org.apache.commons.math3.ode.ContinuousOutputModel;

import java.util.List;
import java.util.concurrent.ConcurrentHashMap;

/**
 * This class is a lightweight wrapper of the result of the Flow ODE integration. It allows to easily query the
 * flow at every point in time and to use the flow to efficiently integrate over a time span.
 * It supports intervals and also reset of the initial state at each interval start.
 */
public class Flow implements IFlow {
    ContinuousOutputModel[] outputModels;

    List<InitialState> initialStates;
    boolean wasInitialStateResetAtEachInterval;
    int n;

    ConcurrentHashMap<Double, RealMatrix>[] flowCache;

    public Flow(ContinuousOutputModel[] outputModels, int n, List<InitialState> initialStates, boolean wasInitialStateResetAtEachInterval) {
        this.outputModels = outputModels;
        this.n = n;
        this.wasInitialStateResetAtEachInterval = wasInitialStateResetAtEachInterval;
        this.initialStates = initialStates;

        this.flowCache = new ConcurrentHashMap[outputModels.length];
        for (int i = 0; i < outputModels.length; i++) {
            this.flowCache[i] = new ConcurrentHashMap<>();
        }
    }

    /**
     * Allows to integrate over an edge of a tree using the pre-computed flow.
     *
     * @param timeStart the time of the node closer to the root.
     * @param timeEnd   the time of the node closer to the leaves.
     * @param endState  the initial state at the node closer to the leaves.
     * @return the integration result at the time of the node closer to the root.
     */
    @Override
    public IntegrationResult integrateUsingFlow(double timeStart, double timeEnd, double[] endState) {
        int intervalEnd = this.getLeftInterval(timeEnd);
        RealMatrix flowMatrixEnd = this.getFlow(intervalEnd, timeEnd);

        RealVector likelihoodVectorEnd = Utils.toVector(endState);

        DecompositionSolver linearSolver = new QRDecomposition(flowMatrixEnd, 1e-10).getSolver();

        RealVector solution = null;
        try {
            solution = linearSolver.solve(likelihoodVectorEnd);
        } catch (SingularMatrixException e) {
            // we fall back to an SVD least-squares solver
            SingularValueDecomposition svd = new SingularValueDecomposition(flowMatrixEnd);

            if (Double.isInfinite(svd.getConditionNumber())) {
                throw new IllegalStateException("Infinite condition number found.");
            }

            linearSolver = svd.getSolver();
            solution = linearSolver.solve(likelihoodVectorEnd);
        }

        return this.operateFlow(timeStart, intervalEnd, solution);
    }

    /**
     * Operates the flow at a given time on the given vector.
     * This method supports when the flow integration was restarted using the same initial state
     * at the beginning of every interval. In this case, the flow is calculated by accumulatively
     * multiplying the end flows of the intervals between startingAtInterval and time.
     *
     * @param time               the time for which to query the flow from.
     * @param startingAtInterval where to start the accumulation of the flow if initial state resetting
     *                           was used.
     * @param vector             the vector to multiply the flow with.
     * @return the flow at the given time multiplied with the vector.
     */
    public IntegrationResult operateFlow(double time, int startingAtInterval, RealVector vector) {
        int timeInterval = this.getRightInterval(time);

        RealVector accumulatedVector = vector;
        double logScalingFactor = Utils.rescale(accumulatedVector, 0.0);

        for (int i = startingAtInterval; i < timeInterval ; i++) {
            RealMatrix flowEnd = this.getFlow(i, this.outputModels[i].getFinalTime());

            accumulatedVector = flowEnd.operate(accumulatedVector);
            logScalingFactor = Utils.rescale(accumulatedVector, logScalingFactor);

            accumulatedVector = this.initialStates.get(this.initialStates.size() - i - 2).inverse().operate(accumulatedVector);
            logScalingFactor = Utils.rescale(accumulatedVector, logScalingFactor);
        }

        accumulatedVector = this.getFlow(timeInterval, time).operate(accumulatedVector);
        logScalingFactor = Utils.rescale(accumulatedVector, logScalingFactor);

        return new IntegrationResult(accumulatedVector.toArray(), logScalingFactor);
    }

    RealMatrix getFlow(int interval, double time) {
        RealMatrix flow = this.flowCache[interval].get(time);

        if (flow == null) {
            ContinuousOutputModel output = this.outputModels[interval];

            synchronized (output) {
                output.setInterpolatedTime(time);
                flow = Utils.toMatrix(output.getInterpolatedState(), n);
            }

            this.flowCache[interval].put(time, flow);
        }

        return flow;
    }

    /**
     * Returns the interval corresponding to the given time.
     *
     * @param time the time to get the interval for.
     * @return the interval.
     */
    public int getLeftInterval(double time) {
        if (this.outputModels[0].getInitialTime() < time) {
            return 0;
        }

        for (int i = 0; i < this.outputModels.length; i++) {
            ContinuousOutputModel model = this.outputModels[i];

            if (model.getFinalTime() < time && time <= model.getInitialTime()) {
                return i;
            }
        }

        return this.outputModels.length - 1;
    }

    /**
     * Returns the interval corresponding to the given time.
     *
     * @param time the time to get the interval for.
     * @return the interval.
     */
    public int getRightInterval(double time) {
        if (this.outputModels[0].getInitialTime() <= time) {
            return 0;
        }

        for (int i = 0; i < this.outputModels.length; i++) {
            ContinuousOutputModel model = this.outputModels[i];

            if (model.getFinalTime() <= time && time < model.getInitialTime()) {
                return i;
            }
        }

        return this.outputModels.length - 1;
    }
}
