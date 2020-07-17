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

    public abstract double applyToTrajectory(Parameterization param, int interval, Trajectory trajectory);
}
