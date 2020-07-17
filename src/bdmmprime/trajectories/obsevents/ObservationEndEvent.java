package bdmmprime.trajectories.obsevents;

import bdmmprime.parameterization.Parameterization;
import bdmmprime.trajectories.Trajectory;

public class ObservationEndEvent extends ObservedEvent {

    public ObservationEndEvent(double time) {
        this.time = time;
    }

    @Override
    public int[] getNextLineageCounts() {
        throw new RuntimeException("getNextLineageCounts() not implemented for ObservationEndEvent!");
    }

    @Override
    public double applyToTrajectory(Parameterization param, int interval, Trajectory trajectory) {
        throw new RuntimeException("applyToTrajectory() not implemented for ObservationEndEvent!");
    }

    @Override
    public boolean isFinalEvent() {
        return true;
    }
}
