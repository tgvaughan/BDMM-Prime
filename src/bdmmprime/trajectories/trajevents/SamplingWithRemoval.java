package bdmmprime.trajectories.trajevents;

public class SamplingWithRemoval extends TrajectoryEvent {

    int type;

    public SamplingWithRemoval(double time, int type, int multiplicity) {
        this.time = time;
        this.type = type;
        this.multiplicity = multiplicity;
    }

    public SamplingWithRemoval(double time, int type) {
        this(time, type, 1);
    }

    @Override
    public void updateState(double[] state) {
        state[type] -= multiplicity;
    }

    @Override
    public void reverseUpdateState(double[] state) {
        state[type] += multiplicity;
    }

    @Override
    public String toString() {
        return "SamplingWithRemoval{" +
                "type=" + type +
                ", multiplicity=" + multiplicity +
                '}';
    }
}
