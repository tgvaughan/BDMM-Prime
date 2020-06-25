package bdmmprime.trajectories;

public class MigrationEvent extends TrajectoryEvent {

    int srcType, destType;

    public MigrationEvent(double time, int srcType, int destType) {
        this.time = time;
        this.srcType = srcType;
        this.destType = destType;
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
}
