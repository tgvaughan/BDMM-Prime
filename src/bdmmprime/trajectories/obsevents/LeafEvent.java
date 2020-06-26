package bdmmprime.trajectories.obsevents;

public class LeafEvent extends ObservedEvent {
    public int type;

    public LeafEvent(double time, int type, int multiplicity) {
        this.time = time;
        this.type = type;
        this.multiplicity = multiplicity;
    }
}
