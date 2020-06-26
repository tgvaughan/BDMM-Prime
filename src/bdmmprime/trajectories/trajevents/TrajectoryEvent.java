package bdmmprime.trajectories.trajevents;

public abstract class TrajectoryEvent {
    public double time;
    public int multiplicity;

    public abstract void updateState(double[] state);

    public abstract void reverseUpdateState(double[] state);
}
