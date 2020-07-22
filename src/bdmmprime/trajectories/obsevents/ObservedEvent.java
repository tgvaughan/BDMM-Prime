package bdmmprime.trajectories.obsevents;

import bdmmprime.parameterization.Parameterization;
import bdmmprime.trajectories.Trajectory;

public abstract class ObservedEvent {

    public double time;
    public int type, multiplicity;
    public int[] lineages; // before the event

    public boolean isFinalEvent() {
        return false;
    }

    public abstract int[] getNextLineageCounts();

    /**
     * Apply the observed event to the chosen trajectory, returning the probability of the observed
     * event given the trajectory.
     *
     * @param param parameterization object
     * @param interval skyline interval
     * @param trajectory trajectory to update
     * @return log probability of observed event given trajectory.
     */
    public abstract double applyToTrajectory(Parameterization param, int interval, Trajectory trajectory);
}
