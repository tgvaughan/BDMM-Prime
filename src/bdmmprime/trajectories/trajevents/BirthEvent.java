package bdmmprime.trajectories.trajevents;

public class BirthEvent extends TrajectoryEvent {

    int type;

    public BirthEvent(double time, int type, int multiplicity) {
        this.time = time;
        this.type = type;
        this.multiplicity = multiplicity;
    }

    public BirthEvent(double time, int type) {
        this(time, type, 1);
    }

    @Override
    public void updateState(double[] state) {
        state[type] += 1;
    }

    @Override
    public void reverseUpdateState(double[] state) {
        state[type] -= 1;
    }

    @Override
    public String toString() {
        return "BirthEvent{" +
                "type=" + type +
                '}';
    }
}
