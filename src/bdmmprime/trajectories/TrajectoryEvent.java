package bdmmprime.trajectories;

public abstract class TrajectoryEvent {
    double time;

    public abstract void updateState(double[] state);

    public abstract void reverseUpdateState(double[] state);
}
