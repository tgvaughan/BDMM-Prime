package bdmmprime.trajectories.obsevents;

public class SampledAncestorEvent extends ObservedEvent {

    public int type;

    public SampledAncestorEvent(double time, int type, int multiplicity) {
        this.time = time;
        this.type = type;
        this.multiplicity = multiplicity;
    }
}
