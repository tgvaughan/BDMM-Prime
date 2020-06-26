package bdmmprime.trajectories.trajevents;

public class DeathEvent extends TrajectoryEvent {

    int type;

    public DeathEvent(double time, int type, int multiplicity) {
        this.time = time;
        this.type = type;
        this.multiplicity = multiplicity;
    }

    public DeathEvent(double time, int type) {
        this(time, type, 1);
    }

    @Override
    public void updateState(double[] state) {
        state[type] -= 1;
    }

    public void reverseUpdateState(double[] state) {
        state[type] += 1;
    }

    @Override
    public String toString() {
        return "DeathEvent{" +
                "type=" + type +
                '}';
    }
}
