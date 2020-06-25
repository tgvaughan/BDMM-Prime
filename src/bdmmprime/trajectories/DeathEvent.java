package bdmmprime.trajectories;

public class DeathEvent extends TrajectoryEvent {

    int type;

    public DeathEvent(double time, int type) {
        this.time = time;
        this.type = type;
    }

    @Override
    public void updateState(double[] state) {
        state[type] -= 1;
    }

    public void reverseUpdateState(double[] state) {
        state[type] += 1;
    }
}
