package bdmmprime.trajectories;

public class BirthEvent extends TrajectoryEvent {

    int type;

    public BirthEvent(double time, int type) {
        this.time = time;
        this.type = type;
    }

    @Override
    public void updateState(double[] state) {
        state[type] += 1;
    }

    @Override
    public void reverseUpdateState(double[] state) {
        state[type] -= 1;
    }
}
