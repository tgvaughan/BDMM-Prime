package bdmmprime.trajectories;

public class SamplingWithoutRemoval extends TrajectoryEvent {

    int type;

    public SamplingWithoutRemoval(double time, int type) {
        this.time = time;
        this.type = type;
    }

    @Override
    public void updateState(double[] state) { }

    @Override
    public void reverseUpdateState(double[] state) { }
}
