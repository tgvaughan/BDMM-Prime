package bdmmprime.trajectories.trajevents;

import bdmmprime.trajectories.obsevents.SampledAncestorEvent;

public class SamplingWithoutRemoval extends TrajectoryEvent {

    int type;

    public SamplingWithoutRemoval(double time, int type, int multiplicity) {
        this.time = time;
        this.type = type;
        this.multiplicity = multiplicity;
    }

    public SamplingWithoutRemoval(double time, int type) {
        this(time, type, 1);
    }

    @Override
    public void updateState(double[] state) { }

    @Override
    public void reverseUpdateState(double[] state) { }

    @Override
    public String toString() {
        return "SamplingWithoutRemoval{" +
                "type=" + type +
                ", multiplicity=" + multiplicity +
                '}';
    }
}
