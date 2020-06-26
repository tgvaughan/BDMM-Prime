package bdmmprime.trajectories.trajevents;

public class MigrationEvent extends TrajectoryEvent {

    int srcType, destType;

    public MigrationEvent(double time, int srcType, int destType, int multiplicity) {
        this.time = time;
        this.srcType = srcType;
        this.destType = destType;
        this.multiplicity = multiplicity;
    }

    public MigrationEvent(double time, int srcType, int destType) {
        this(time, srcType, destType, 1);
    }

    @Override
    public void updateState(double[] state) {
        state[srcType] -= 1;
        state[destType] += 1;
    }

    @Override
    public void reverseUpdateState(double[] state) {
        state[srcType] += 1;
        state[destType] -= 1;
    }

    @Override
    public String toString() {
        return "MigrationEvent{" +
                "srcType=" + srcType +
                ", destType=" + destType +
                '}';
    }
}
