package bdmmprime.trajectories;

public class CrossBirthEvent extends TrajectoryEvent {

    int srcType, destType;

    public CrossBirthEvent(double time, int srcType, int destType) {
        this.time = time;
        this.srcType = srcType;
        this.destType = destType;
    }

    @Override
    public void updateState(double[] state) {
        state[destType] += 1;
    }
}
