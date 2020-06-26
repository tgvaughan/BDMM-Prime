package bdmmprime.trajectories.trajevents;

public class CrossBirthEvent extends TrajectoryEvent {

    int srcType, destType;

    public CrossBirthEvent(double time, int srcType, int destType, int multiplicity) {
        this.time = time;
        this.srcType = srcType;
        this.destType = destType;
        this.multiplicity = multiplicity;
    }

    public CrossBirthEvent(double time, int srcType, int destType) {
        this(time, srcType, destType, 1);
    }

    @Override
    public void updateState(double[] state) {
        state[destType] += 1;
    }

    @Override
    public void reverseUpdateState(double[] state) {
        state[destType] -= 1;
    }

    @Override
    public String toString() {
        return "CrossBirthEvent{" +
                "srcType=" + srcType +
                ", destType=" + destType +
                '}';
    }
}
